package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import biz.k11i.xgboost.Predictor;
import biz.k11i.xgboost.learner.ObjFunction;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
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
import scala.tools.ant.sabbus.Break;

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

    static private final int NUM_OVERLAP = 0;
    static private final int NUM_COHERENT = 1;
    static private final int OVERLAP_MAPPING_QUALITY = 2;
    static private final int COHERENT_MAPPING_QUALITY = 3;

    private final PartitionCrossingChecker partitionCrossingChecker;

    private final Predictor predictor;
    private final Map<String, Integer> evidenceTypeMap;
    private final double coverage;
    private final double thresholdProbability;
    private final Map<String, String> readGroupToLibraryMap;
    private final Map<String, LibraryStatistics> libraryStatisticsMap;

    private final EvidenceOverlapChecker evidenceOverlapChecker;
    private final Map<BreakpointEvidence, int[]> rawFeatureCache;

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
        this.partitionCrossingChecker = partitionCrossingChecker;
        this.coverage = readMetadata.getCoverage();
        this.thresholdProbability = params.svEvidenceFilterThresholdProbability;
        this.readGroupToLibraryMap = readMetadata.getReadGroupToLibraryMap();
        this.libraryStatisticsMap = readMetadata.getAllLibraryStatistics();

        this.evidenceOverlapChecker = new EvidenceOverlapChecker(evidenceItr, readMetadata, params.minEvidenceMapQ);
        this.rawFeatureCache = new HashMap<>();

        this.listItr = null;
        this.treeItr = evidenceOverlapChecker.getTreeItr();
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

    private static InputStream getInputStream(final String fileLocation) {
        return fileLocation.startsWith("gatk-resources::") ?
                XGBoostEvidenceFilter.class.getResourceAsStream(fileLocation.substring(16))
                : BucketUtils.openFile(fileLocation);
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
        final double mappingQuality = (double)getMappingQuality(evidence);
        final double templateSize = getTemplateSize(evidence);
        // calculate these similar to BreakpointDensityFilter, but always calculate full totals, never end early.
        final CoverageScaledOverlapInfo individualOverlapInfo = getIndividualOverlapInfo(evidence);
        final CoverageScaledOverlapInfo clusterOverlapInfo = getClusterOverlapInfo(evidence);

        return new EvidenceFeatures(
            new double[]{
                cigarQualityInfo.basesMatched, cigarQualityInfo.referenceLength, evidenceType, mappingQuality, templateSize,
                individualOverlapInfo.numOverlap, individualOverlapInfo.overlapMappingQuality,
                individualOverlapInfo.meanOverlapMappingQuality, individualOverlapInfo.numCoherent,
                individualOverlapInfo.coherentMappingQuality,
                clusterOverlapInfo.numOverlap, clusterOverlapInfo.overlapMappingQuality,
                clusterOverlapInfo.meanOverlapMappingQuality, clusterOverlapInfo.numCoherent,
                clusterOverlapInfo.coherentMappingQuality
            }
        );
    }

    private int getMappingQuality(final BreakpointEvidence evidence) {
        // Note: return "max" mapping quality for non-ReadEvidence. Reasoning: some features depend on sum or average of
        // read qualities. Non-ReadEvidence isn't *bad* per se, so give it a good score.
        return evidence instanceof BreakpointEvidence.ReadEvidence ?
            ((BreakpointEvidence.ReadEvidence) evidence).getMappingQuality() : 60;
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

    private void cacheFeatures(final BreakpointEvidence evidence) {
        /*
        final int evidenceMappingQuality = getMappingQuality(evidence);
        int numOverlap = -1;
        int overlapMappingQuality = -evidenceMappingQuality;
        for(final Iterator<BreakpointEvidence> overlapper_itr = evidenceOverlapChecker.overlappers(evidence); overlapper_itr.hasNext();) {
            final BreakpointEvidence overlapper = overlapper_itr.next();
            if(overlapper == evidence) {
                continue; // don't count self-overlap
            }
            numOverlap += 1;
            overlapMappingQuality += getMappingQuality(overlapper);
        }

        int numCoherent = 0;
        int coherentMappingQuality = 0;
        for(final Iterator<BreakpointEvidence> coherent_itr = evidenceOverlapChecker.coherent(evidence); coherent_itr.hasNext();) {
            final BreakpointEvidence coherent = coherent_itr.next();
            if(coherent == evidence) {
                continue; // don't count self-coherence
            }
            numCoherent += 1;
            coherentMappingQuality += getMappingQuality(coherent);
        }
        */
        int numOverlap = 0;
        int overlapMappingQuality = 0;
        int numCoherent = 0;
        int coherentMappingQuality = 0;
        for(final Iterator<Tuple2<BreakpointEvidence, Boolean>> overlapper_itr = evidenceOverlapChecker.overlappersWithCoherence(evidence);
                overlapper_itr.hasNext();) {
            final Tuple2<BreakpointEvidence, Boolean> itr_results = overlapper_itr.next();
            final BreakpointEvidence overlapper = itr_results._1;
            if(overlapper == evidence) {
                continue; // don't count self-overlap
            }
            numOverlap += 1;
            final int mappingQuality = getMappingQuality(overlapper);
            overlapMappingQuality += mappingQuality;

            final boolean isCoherent = itr_results._2;
            if(isCoherent) {
                numCoherent += 1;
                coherentMappingQuality += mappingQuality;
            }
        }
        final int[] evidenceFeatureCache = new int[4];
        evidenceFeatureCache[NUM_OVERLAP] = numOverlap;
        evidenceFeatureCache[NUM_COHERENT] = numCoherent;
        evidenceFeatureCache[OVERLAP_MAPPING_QUALITY] = overlapMappingQuality;
        evidenceFeatureCache[COHERENT_MAPPING_QUALITY] = coherentMappingQuality;
        rawFeatureCache.put(evidence, evidenceFeatureCache);
    }

    private CoverageScaledOverlapInfo getIndividualOverlapInfo(final BreakpointEvidence evidence) {
        if(!rawFeatureCache.containsKey(evidence)) {
            cacheFeatures(evidence);
        }
        final int[] evidenceFeatureCache = rawFeatureCache.get(evidence);
        final int numOverlap = evidenceFeatureCache[NUM_OVERLAP];
        final int numCoherent = evidenceFeatureCache[NUM_COHERENT];
        final int overlapMappingQuality = evidenceFeatureCache[OVERLAP_MAPPING_QUALITY];
        final int coherentMappingQuality = evidenceFeatureCache[COHERENT_MAPPING_QUALITY];
        final double meanOverlapMappingQuality = ((double)overlapMappingQuality) / numOverlap;
        return new CoverageScaledOverlapInfo(numOverlap, numCoherent, overlapMappingQuality, coherentMappingQuality,
                meanOverlapMappingQuality, coverage);
    }

    private CoverageScaledOverlapInfo getClusterOverlapInfo(final BreakpointEvidence evidence) {
        int clusterNumOverlap = 0;
        int clusterNumCoherent = 0;
        int clusterOverlapMappingQuality = 0;
        int clusterCoherentMappingQuality = 0;
        double clusterMeanOverlapMappingQuality = 0.0;
        for (final Iterator<BreakpointEvidence> overlapper_itr = evidenceOverlapChecker.overlappers(evidence); overlapper_itr.hasNext(); ) {
            final BreakpointEvidence overlapper = overlapper_itr.next();
            if (overlapper == evidence) {
                continue; // don't count self-overlap in cluster features
            }
            if(!rawFeatureCache.containsKey(overlapper)) {
                cacheFeatures(overlapper);
            }
            final int[] overlapperFeatureCache = rawFeatureCache.get(overlapper);
            clusterNumOverlap = Math.max(clusterNumOverlap, overlapperFeatureCache[NUM_OVERLAP]);
            clusterNumCoherent = Math.max(clusterNumCoherent, overlapperFeatureCache[NUM_COHERENT]);
            clusterOverlapMappingQuality = Math.max(clusterOverlapMappingQuality, overlapperFeatureCache[OVERLAP_MAPPING_QUALITY]);
            clusterCoherentMappingQuality = Math.max(clusterCoherentMappingQuality, overlapperFeatureCache[COHERENT_MAPPING_QUALITY]);

            final double meanOverlapMappingQuality = ((double) overlapperFeatureCache[OVERLAP_MAPPING_QUALITY]) / overlapperFeatureCache[NUM_OVERLAP];
            clusterMeanOverlapMappingQuality = Math.max(clusterMeanOverlapMappingQuality, meanOverlapMappingQuality);
        }

        return new CoverageScaledOverlapInfo(clusterNumOverlap, clusterNumCoherent, clusterOverlapMappingQuality,
                clusterCoherentMappingQuality, clusterMeanOverlapMappingQuality, coverage);
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

    private static class EvidenceOverlapChecker {
        private final SVIntervalTree<List<BreakpointEvidence>> evidenceTree;
        //private final PairedStrandedIntervalTree<BreakpointEvidence> targetIntervalTree;
        private final ReadMetadata readMetadata;
        private final int minEvidenceMapQ;

        EvidenceOverlapChecker(final Iterator<BreakpointEvidence> evidenceItr, final ReadMetadata readMetadata,
                               final int minEvidenceMapQ) {
            this.readMetadata = readMetadata;
            this.minEvidenceMapQ = minEvidenceMapQ;
            evidenceTree = new SVIntervalTree<>();
            //targetIntervalTree = new PairedStrandedIntervalTree<>();
            while ( evidenceItr.hasNext() ) {
                final BreakpointEvidence evidence = evidenceItr.next();
                addToTree(evidenceTree, evidence.getLocation(), evidence);
                /*
                final Boolean strand = evidence.isEvidenceUpstreamOfBreakpoint();
                if (strand != null && evidence.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
                    final StrandedInterval strandedEvidence = new StrandedInterval(evidence.getLocation(), strand);
                    for(final StrandedInterval distalTarget : evidence.getDistalTargets(readMetadata, minEvidenceMapQ)) {
                        targetIntervalTree.put(new PairedStrandedIntervals(strandedEvidence, distalTarget), evidence);
                    }
                }
                */
            }
        }

        Iterator<BreakpointEvidence> overlappers(final BreakpointEvidence evidence) {
            return new OverlapperIterator(evidence, evidenceTree);
        }
        /*
        Iterator<BreakpointEvidence> coherent(final BreakpointEvidence evidence) {
            return new CoherentIterator(evidence, targetIntervalTree, readMetadata, minEvidenceMapQ);
        }
        */
        Iterator<Tuple2<BreakpointEvidence, Boolean>> overlappersWithCoherence(final BreakpointEvidence evidence) {
            return new OverlapAndCoherenceIterator(evidence, evidenceTree, readMetadata, minEvidenceMapQ);
        }

        Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> getTreeItr() {
            return evidenceTree.iterator();
        }

        private static class OverlapperIterator  implements Iterator<BreakpointEvidence> {
            private final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
            private Iterator<BreakpointEvidence> listItr;

            OverlapperIterator(final BreakpointEvidence evidence,
                               final SVIntervalTree<List<BreakpointEvidence>> evidenceTree) {
                treeItr = evidenceTree.overlappers(evidence.getLocation());
                listItr = null;
            }

            @Override
            public boolean hasNext() {
                if ( listItr != null && listItr.hasNext() ) {
                    return true;
                }
                while(treeItr.hasNext()) {
                    listItr = treeItr.next().getValue().iterator();
                    if(listItr.hasNext()) {
                        return true;
                    }
                }
                return false;
            }

            @Override
            public BreakpointEvidence next() {
                if ( !hasNext() ) {
                    throw new NoSuchElementException("No next element.");
                }
                return listItr.next();
            }
        }

        private static class OverlapAndCoherenceIterator implements Iterator<Tuple2<BreakpointEvidence, Boolean>> {
            private final ReadMetadata readMetadata;
            private final int minEvidenceMapQ;
            private final Iterator<SVIntervalTree.Entry<List<BreakpointEvidence>>> treeItr;
            private Iterator<BreakpointEvidence> listItr;
            private final boolean checkCoherence;
            private final Boolean strand;
            private final List<StrandedInterval> distalTargets;

            OverlapAndCoherenceIterator(final BreakpointEvidence evidence,
                                        final SVIntervalTree<List<BreakpointEvidence>> evidenceTree,
                                        final ReadMetadata readMetadata,
                                        final int minEvidenceMapQ) {
                this.readMetadata = readMetadata;
                this.minEvidenceMapQ = minEvidenceMapQ;
                this.treeItr = evidenceTree.overlappers(evidence.getLocation());
                this.strand = evidence.isEvidenceUpstreamOfBreakpoint();
                this.checkCoherence = (strand != null) && evidence.hasDistalTargets(readMetadata, minEvidenceMapQ);
                if(this.checkCoherence) {
                    distalTargets = new ArrayList<>(evidence.getDistalTargets(readMetadata, minEvidenceMapQ));
                } else {
                    distalTargets = null;
                }
            }

            @Override
            public boolean hasNext() {
                if ( listItr != null && listItr.hasNext() ) {
                    return true;
                }
                while(treeItr.hasNext()) {
                    listItr = treeItr.next().getValue().iterator();
                    if(listItr.hasNext()) {
                        return true;
                    }
                }
                return false;
            }

            @Override
            public Tuple2<BreakpointEvidence, Boolean> next() {
                if ( !hasNext() ) {
                    throw new NoSuchElementException("No next element.");
                }
                final BreakpointEvidence overlapper = listItr.next();
                Boolean isCoherent = false;
                if(checkCoherence) {
                    Boolean overlapperStrand = overlapper.isEvidenceUpstreamOfBreakpoint();
                    if(overlapperStrand == strand && overlapper.hasDistalTargets(readMetadata, minEvidenceMapQ)) {
                        for(final StrandedInterval distalTarget : distalTargets) {
                            for(final StrandedInterval overlapperDistalTarget : new ArrayList<>(overlapper.getDistalTargets(readMetadata, minEvidenceMapQ))) {
                                if(distalTarget.getStrand() == overlapperDistalTarget.getStrand()
                                        && distalTarget.getInterval().overlaps(overlapperDistalTarget.getInterval())) {
                                    isCoherent = true; // could break, but everything has 1 distal target, so false optimization
                                }
                            }
                        }
                    }
                }
                return new Tuple2<>(overlapper, isCoherent);
            }
        }

        private static class CoherentIterator  implements Iterator<BreakpointEvidence> {
            private Iterator<Tuple2<PairedStrandedIntervals, BreakpointEvidence>> treeItr;
            private final Iterator<StrandedInterval> distalTargetItr;
            private final StrandedInterval strandedEvidence;
            private final PairedStrandedIntervalTree<BreakpointEvidence> coherentTree;

            CoherentIterator(final BreakpointEvidence evidence,
                             final PairedStrandedIntervalTree<BreakpointEvidence> coherentTree,
                             final ReadMetadata readMetadata,
                             final int minEvidenceMapq) {
                treeItr = null;
                final Boolean strand = evidence.isEvidenceUpstreamOfBreakpoint();
                if (strand != null && evidence.hasDistalTargets(readMetadata, minEvidenceMapq)) {
                    distalTargetItr = evidence.getDistalTargets(readMetadata, minEvidenceMapq).iterator();
                    strandedEvidence = new StrandedInterval(evidence.getLocation(), strand);
                    this.coherentTree = coherentTree;
                } else {
                    distalTargetItr = null;
                    strandedEvidence = null;
                    this.coherentTree = null;
                }
            }

            @Override
            public boolean hasNext() {
                if ( treeItr != null && treeItr.hasNext() ) {
                    return true;
                }
                if(distalTargetItr == null) {
                    return false;
                }
                while(distalTargetItr.hasNext()) {
                    treeItr = Lists.newArrayList(
                            coherentTree.overlappers(
                                    new PairedStrandedIntervals(strandedEvidence, distalTargetItr.next())
                            )
                    ).iterator();
                    /*
                    treeItr = coherentTree.overlappers(
                        new PairedStrandedIntervals(strandedEvidence, distalTargetItr.next())
                    );
                    */
                    if (treeItr.hasNext()) {
                        return true;
                    }
                }
                return false;
            }

            @Override
            public BreakpointEvidence next() {
                if ( !hasNext() ) {
                    throw new NoSuchElementException("No next element.");
                }
                return treeItr.next()._2;
            }
        }
    }


    private static class CoverageScaledOverlapInfo {
        final double numOverlap;
        final double overlapMappingQuality;
        final double meanOverlapMappingQuality;
        final double numCoherent;
        final double coherentMappingQuality;

        CoverageScaledOverlapInfo(final int numOverlap, final int numCoherent, final int overlapMappingQuality,
                                  final int coherentMappingQuality, final double meanOverlapMappingQuality,
                                  final double coverage) {
            this.numOverlap = ((double)numOverlap) / coverage;
            this.overlapMappingQuality = ((double)overlapMappingQuality) / coverage;
            this.numCoherent = ((double)numCoherent) / coverage;
            this.coherentMappingQuality = ((double)coherentMappingQuality) / coverage;
            this.meanOverlapMappingQuality = meanOverlapMappingQuality;
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
