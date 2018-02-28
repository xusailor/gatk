package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import biz.k11i.xgboost.Predictor;
import biz.k11i.xgboost.learner.ObjFunction;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.*;
import org.broadinstitute.hellbender.tools.spark.utils.IntHistogram;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.IOException;
import java.io.InputStream;
import java.util.*;

/**
 * A class that acts as a filter for breakpoint evidence.
 * It passes only that evidence that is part of a putative cluster.
 */
public final class XGBoostEvidenceFilter implements Iterator<BreakpointEvidence> {
    static private final boolean useFastMathExp = true;
    static private final boolean useClusterNumCoherent = false;
    private static final List<String> defaultEvidenceTypeOrder = Arrays.asList(
            "TemplateSizeAnomaly",
            "MateUnmapped", "InterContigPair",
            "SplitRead", "LargeIndel", "WeirdTemplateSize", "SameStrandPair", "OutiesPair"
    );

    private final ReadMetadata readMetadata;
    private final PartitionCrossingChecker partitionCrossingChecker;

    private final Predictor predictor;
    private final Map<String, Integer> evidenceTypeMap;
    private final double coverage;
    private final int minEvidenceMapQ;
    private final double thresholdProbability;
    private final Map<String, String> readGroupToLibraryMap;
    private final Map<String, LibraryStatistics> libraryStatisticsMap;

    private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;

    private Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
    private Iterator<BreakpointEvidence> listItr;

    public XGBoostEvidenceFilter(
            final Iterator<BreakpointEvidence> evidenceItr,
            final ReadMetadata readMetadata,
            final StructuralVariationDiscoveryArgumentCollection.FindBreakpointEvidenceSparkArgumentCollection params,
            final PartitionCrossingChecker partitionCrossingChecker
    ) {
        this.predictor = loadPredictor(params.svEvidenceFilterModelFile);
        this.evidenceTypeMap = loadEvidenceTypeMap(params.svCategoricalVariablesFile);
        this.readMetadata = readMetadata;
        this.partitionCrossingChecker = partitionCrossingChecker;
        this.coverage = readMetadata.getCoverage();
        this.minEvidenceMapQ = params.minEvidenceMapQ;
        this.thresholdProbability = params.svEvidenceFilterThresholdProbability;
        this.readGroupToLibraryMap = readMetadata.getReadGroupToLibraryMap();
        this.libraryStatisticsMap = readMetadata.getAllLibraryStatistics();

        this.evidenceTree = buildTree(evidenceItr);
        this.treeItr = evidenceTree.iterator();
        this.listItr = null;
    }


    private static InputStream getInputStream(final String fileLocation) {
        return fileLocation.startsWith("gatk-resources::") ?
                XGBoostEvidenceFilter.class.getResourceAsStream(fileLocation.substring(16))
                : BucketUtils.openFile(fileLocation);
    }

    public static Predictor loadPredictor(final String modelFileLocation) {
        ObjFunction.useFastMathExp(useFastMathExp);
        try(final InputStream inputStream = getInputStream(modelFileLocation)) {
            return new Predictor(inputStream);
        } catch(IOException e) {
            throw new GATKException(
                    "Unable to load predictor from classifier file " + modelFileLocation + ": " + e.getMessage()
            );
        }
    }

    @VisibleForTesting
    static Map<String, Integer> loadEvidenceTypeMap(final String categoricalVariablesMapFile) {
        final HashMap<String, Integer> evidenceTypeMap = new HashMap<>();
        if(categoricalVariablesMapFile == null || categoricalVariablesMapFile.isEmpty()) {
            // no categorical variables file, just use default evidence order
            for(int index = 0; index < defaultEvidenceTypeOrder.size(); ++index) {
                evidenceTypeMap.put(defaultEvidenceTypeOrder.get(index), index);
            }
            return evidenceTypeMap;
        }
        // override default evidence order with saved order
        try(final InputStream inputStream = getInputStream(categoricalVariablesMapFile)) {
            final JsonNode testDataNode = new ObjectMapper().readTree(inputStream);
            final JsonNode evidenceTypeArrayNode = testDataNode.get("evidence_type");
            if(!evidenceTypeArrayNode.isArray()) {
                throw new IllegalArgumentException("evidenceTypeNode does not encode a valid array");
            }
            int index = 0;
            for(final JsonNode evidenceTypeNode : evidenceTypeArrayNode) {
                evidenceTypeMap.put(evidenceTypeNode.asText(), index);
                ++index;
            }
            return evidenceTypeMap;
        } catch(IOException e) {
            throw new GATKException(
                    "Unable to load evidence-type map from file " + categoricalVariablesMapFile + ": " + e.getMessage()
            );
        }
    }

    @Override
    public boolean hasNext() {
        if ( listItr != null && listItr.hasNext() ) {
            return true;
        }
        listItr = null;
        boolean result = false;
        while ( !result && treeItr.hasNext() ) {
            final SVIntervalTree.Entry<List<BreakpointEvidence>> entry = treeItr.next();
            final SVInterval curInterval = entry.getInterval();
            final List<BreakpointEvidence> evidenceList = entry.getValue();
            if( isValidated(entry.getValue()) || partitionCrossingChecker.onBoundary(curInterval) ) {
                // already validated (no need to mark validated again) or on partition boundary (punt for now)
                result = true;
            } else if( anyPassesFilter(evidenceList) ) {
                evidenceList.forEach(ev -> ev.setValidated(true));
                result = true;
            }

            if ( result ) {
                listItr = entry.getValue().iterator();
            }
        }
        return result;
    }

    @Override
    public BreakpointEvidence next() {
        if ( !hasNext() ) {
            throw new NoSuchElementException("No next element.");
        }
        return listItr.next();
    }

    private static SVIntervalTree<List<BreakpointEvidence>> buildTree( final Iterator<BreakpointEvidence> evidenceItr ) {
        SVIntervalTree<List<BreakpointEvidence>> tree = new SVIntervalTree<>();
        while ( evidenceItr.hasNext() ) {
            final BreakpointEvidence evidence = evidenceItr.next();
            addToTree(tree, evidence.getLocation(), evidence);
        }
        return tree;
    }

    private boolean isValidated( final List<BreakpointEvidence> evList ) {
        for ( final BreakpointEvidence ev : evList ) {
            if ( ev.isValidated() ) return true;
        }
        return false;
    }

    boolean anyPassesFilter(final List<BreakpointEvidence> evidenceList) {
        for(final BreakpointEvidence evidence : evidenceList) {
            final double evidenceGoodProbability = predictor.predictSingle(getFeatures(evidence));
            if(evidenceGoodProbability > thresholdProbability) {
                return true;
            }
        }
        return false;
    }

    @VisibleForTesting
    EvidenceFeatures getFeatures(final BreakpointEvidence evidence) {
        // create new struct for these two, use CigarOperator to update if instanceof ReadEvidence
        final CigarQualityInfo cigarQualityInfo = new CigarQualityInfo(evidence);
        final double evidenceType = evidenceTypeMap.get(evidence.getClass().getSimpleName());
        final double mappingQuality = getMappingQuality(evidence);
        final double templateSize = getTemplateSize(evidence);
        // calculate these similar to BreakpointDensityFilter, but always calculate full totals, never end early.
        final CoverageScaledOverlapInfo overlapInfo
                = useClusterNumCoherent ? getClusterOverlapInfo(evidence) : getIndividualOverlapInfo(evidence);

        return new EvidenceFeatures(new double[]{cigarQualityInfo.basesMatched, cigarQualityInfo.referenceLength,
                                                 evidenceType, mappingQuality, templateSize,
                                                 overlapInfo.numOverlap, overlapInfo.numCoherent});
    }

    private double getMappingQuality(final BreakpointEvidence evidence) {
        return evidence instanceof BreakpointEvidence.ReadEvidence ?
            ((BreakpointEvidence.ReadEvidence) evidence).getMappingQuality() : 0.0;
    }

    private double getTemplateSize(final BreakpointEvidence evidence) {
        if(evidence instanceof BreakpointEvidence.ReadEvidence) {
            // For ReadEvidence, return templateSize as percentile of library's cumulative density function:
            final String readGroup = ((BreakpointEvidence.ReadEvidence) evidence).getReadGroup();
            final String library = readGroupToLibraryMap.get(readGroup);
            final LibraryStatistics libraryStatistics = libraryStatisticsMap.get(library);
            final IntHistogram.CDF templateSizeCDF = libraryStatistics.getCDF();
            /*
            final IntHistogram.CDF templateSizeCDF = libraryStatisticsMap.get(
                    ((BreakpointEvidence.ReadEvidence) evidence).getReadGroup()
            ).getCDF();
            */
            final int cdfBin = Integer.min(Math.abs(((BreakpointEvidence.ReadEvidence) evidence).getTemplateSize()),
                                           templateSizeCDF.size() - 1);
            return templateSizeCDF.getFraction(cdfBin);
        } else if(evidence instanceof BreakpointEvidence.TemplateSizeAnomaly) {
            // For TemplateSizeAnomaly, return readCount scaled by coverage
            return (double)(((BreakpointEvidence.TemplateSizeAnomaly) evidence).getReadCount()) / coverage;
        } else {
            throw new IllegalArgumentException(
                    "templateSize feature only defined for ReadEvidence and TemplateSizeAnomaly"
            );
        }
    }

    private CoverageScaledOverlapInfo getClusterOverlapInfo(final BreakpointEvidence evidence) {
        final SVInterval interval = evidence.getLocation();
        final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> itr = evidenceTree.overlappers(interval);
        PairedStrandedIntervalTree<BreakpointEvidence> targetIntervalTree = new PairedStrandedIntervalTree<>();
        int numOverlap = -1; // correct for fact that evidence will overlap with itself
        int numCoherent = 0;

        while (itr.hasNext()) {
            final List<BreakpointEvidence> evidenceForInterval = itr.next().getValue();
            numOverlap += evidenceForInterval.size();
            // store any overlappers with distal targets in PairStrandedIntervalTree
            for (final BreakpointEvidence overlapper : evidenceForInterval) {
                if (overlapper.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
                    final List<StrandedInterval> distalTargets = overlapper.getDistalTargets(readMetadata, minEvidenceMapQ);
                    for (int i = 0; i < distalTargets.size(); i++) {
                        targetIntervalTree.put(
                                new PairedStrandedIntervals(
                                        new StrandedInterval(overlapper.getLocation(), overlapper.isEvidenceUpstreamOfBreakpoint()),
                                        distalTargets.get(i)),
                                overlapper
                        );
                    }
                }
            }
        }

        // Find maximal numCoherent of any BreakpointEvidence in PairStrandedIntervalTree
        final Iterator<Tuple2<PairedStrandedIntervals, BreakpointEvidence>> targetLinkIterator = targetIntervalTree.iterator();
        while (targetLinkIterator.hasNext()) {
            Tuple2<PairedStrandedIntervals, BreakpointEvidence> next = targetLinkIterator.next();
            // subtract 1 to offset for coherence with self
            final int targetNumCoherent = (int) Utils.stream(targetIntervalTree.overlappers(next._1())).count() - 1;
            if (targetNumCoherent > numCoherent) {
                numCoherent = targetNumCoherent;
            }
        }

        return new CoverageScaledOverlapInfo(numOverlap, numCoherent, coverage);
    }


    private CoverageScaledOverlapInfo getIndividualOverlapInfo(final BreakpointEvidence evidence) {
        final SVInterval interval = evidence.getLocation();
        final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> itr = evidenceTree.overlappers(interval);
        int numOverlap = -1; // correct for fact that evidence will overlap with itself
        int numCoherent;
        if (!evidence.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
            numCoherent = 0;
            while (itr.hasNext()) {
                numOverlap += itr.next().getValue().size();
            }
            return new CoverageScaledOverlapInfo(numOverlap, numCoherent, coverage);
        }

        numCoherent = -1; // correct for fact that evidence will cohere with itself
        final boolean strand = evidence.isEvidenceUpstreamOfBreakpoint();
        List<StrandedInterval> distalTargets = evidence.getDistalTargets(readMetadata, minEvidenceMapQ);
        PairedStrandedIntervalTree<BreakpointEvidence> targetIntervalTree = new PairedStrandedIntervalTree<>();

        while (itr.hasNext()) {
            final List<BreakpointEvidence> evidenceForInterval = itr.next().getValue();
            numOverlap += evidenceForInterval.size();

            for (final BreakpointEvidence overlapper : evidenceForInterval) {
                if (overlapper.isEvidenceUpstreamOfBreakpoint() != strand
                        || !overlapper.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
                    continue;
                }
                final List<StrandedInterval> overlapperDistalTargets
                        = overlapper.getDistalTargets(readMetadata, minEvidenceMapQ);
                for (final StrandedInterval distalTarget : distalTargets) {
                    for (final StrandedInterval overlapperDistalTarget : overlapperDistalTargets) {
                        if (distalTarget.getStrand() == overlapperDistalTarget.getStrand()
                                && distalTarget.getInterval().overlaps(overlapperDistalTarget.getInterval())) {
                            numCoherent += 1;
                        }
                    }
                }
            }
        }

        return new CoverageScaledOverlapInfo(numOverlap, numCoherent, coverage);
    }


    private static <T> void addToTree( final SVIntervalTree<List<T>> tree,
                                   final SVInterval interval,
                                   final T value ) {
        final SVIntervalTree.Entry<List<T>> entry = tree.find(interval);
        if ( entry != null ) {
            entry.getValue().add(value);
        } else {
            final List<T> valueList = new ArrayList<>(1);
            valueList.add(value);
            tree.put(interval, valueList);
        }
    }

    private static class CoverageScaledOverlapInfo {
        final double numOverlap;
        final double numCoherent;

        CoverageScaledOverlapInfo(final int numOverlap, final int numCoherent, final double coverage) {
            this.numOverlap = ((double)numOverlap) / coverage;
            this.numCoherent = ((double)numCoherent) / coverage;
        }
    }

    private static class CigarQualityInfo {
        final int basesMatched;
        final int referenceLength;

        CigarQualityInfo(final BreakpointEvidence evidence) {
            int numMatched = 0;
            int refLength = 0;
            if(evidence instanceof BreakpointEvidence.ReadEvidence) {
                final String cigarString = ((BreakpointEvidence.ReadEvidence) evidence).getCigarString();
                for (final CigarElement element : TextCigarCodec.decode(cigarString).getCigarElements()) {
                    final CigarOperator op = element.getOperator();
                    if (op.consumesReferenceBases()) {
                        refLength += element.getLength();
                        if (op.consumesReadBases()) {
                            numMatched += element.getLength();
                        }
                    }
                }
            }
            basesMatched = numMatched;
            referenceLength = refLength;
        }
    }
}
