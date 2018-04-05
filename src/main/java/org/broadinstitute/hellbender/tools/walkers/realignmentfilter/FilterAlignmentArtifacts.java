package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.commons.math3.util.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.Trilean;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Filter false positive alignment artifacts from a VCF callset.</p>
 *
 * <p>
 *     FilterAlignmentArtifacts identifies alignment artifacts, that is, apparent variants due to reads being mapped to the wrong genomic locus.
 * </p>
 * <p>
 *     Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome
 *     to confuse the alignment algorithm.  This can occur when the aligner for whatever reason overestimate how uniquely a read
 *     maps, thereby assigning it too high of a mapping quality.  It can also occur through no fault of the aligner due to gaps in
 *     the reference, which can also hide the true position to which a read should map.  By using a good alignment algorithm
 *     (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original
 *     bam alignment) and mapping to the best available reference we can avoid these pitfalls.  The last point is especially important:
 *     one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which
 *     the bam was aligned.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 * <p>
 *     The bam input to this tool should be the reassembly bamout produced by HaplotypeCaller or Mutect2 in the process of generating
 *     the input callset.  The original bam will also work but might fail to filter some indels.  The reference passed with the -R argument
 *     but be the reference to which the input bam was realigned.  This does not need to correspond to the reference of the BWAS-MEM
 *     index image.  The latter should be derived from the best available reference, for example hg38 in humans as of February 2018.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterAlignmentArtifacts \
 *   -R hg19.fasta
 *   -V somatic.vcf.gz \
 *   -I somatic_bamout.bam \
 *   --bwa-mem-index-image hg38.index_image \
 *   -O filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter alignment artifacts from a vcf callset.",
        oneLineSummary = "Filter alignment artifacts from a vcf callset.",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class FilterAlignmentArtifacts extends VariantWalker {
    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    public static final int DEFAULT_INDEL_START_TOLERANCE = 5;
    public static final String INDEL_START_TOLERANCE_LONG_NAME = "indel-start-tolerance";
    @Argument(fullName = INDEL_START_TOLERANCE_LONG_NAME,
            doc="Max distance between indel start of aligned read in the bam and the variant in the vcf", optional=true)
    private int indelStartTolerance = DEFAULT_INDEL_START_TOLERANCE;

    public static final int DEFAULT_FRAGMENT_SIZE = 1000;
    public static final String FRAGMENT_SIZE_LONG_NAME = "fragment-size";
    @Argument(fullName = FRAGMENT_SIZE_LONG_NAME,
            doc="Distance away from variant to look for reads' mates.", optional=true)
    private int fragmentSize = DEFAULT_FRAGMENT_SIZE;

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 20;
    public static final String MIN_MAPPING_QUALITY_LONG_NAME = "min-mapping-quality";
    @Argument(fullName = MIN_MAPPING_QUALITY_LONG_NAME,
            doc="Ignore reads with mapping quality less than this", optional=true)
    private int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;

    public static final int DEFAULT_MAX_FAILED_REALIGNMENTS = 0;
    public static final String MAX_FAILED_REALIGNMENTS_LONG_NAME = "max-failed-realignments";
    @Argument(fullName = MAX_FAILED_REALIGNMENTS_LONG_NAME,
            doc="Maximum number of failed read realignments before a variant is rejected.", optional=true)
    private int maxFailedRealignments = DEFAULT_MAX_FAILED_REALIGNMENTS;

    public static final int DEFAULT_SUFFICIENT_GOOD_REALIGNMENTS = 5;
    public static final String SUFFICIENT_GOOD_REALIGNMENTS_LONG_NAME = "sufficient-good-realignments";
    @Argument(fullName = SUFFICIENT_GOOD_REALIGNMENTS_LONG_NAME,
            doc="Sufficient number of good read realignments to accept a variant.", optional=true)
    private int sufficientGoodRealignments = DEFAULT_SUFFICIENT_GOOD_REALIGNMENTS;


    @ArgumentCollection
    protected RealignmentArgumentCollection realignmentArgumentCollection = new RealignmentArgumentCollection();

    private VariantContextWriter vcfWriter;
    private Realigner realigner;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public void onTraversalStart() {
        realigner = new Realigner(realignmentArgumentCollection);
        vcfWriter = createVCFWriter(new File(outputVcf));

        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME));
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        Trilean passesFilter = Trilean.UNKNOWN;

        if (vc.getNAlleles() == 1) {
            passesFilter = Trilean.TRUE;
        }

        final MutableInt failedRealignmentCount = new MutableInt(0);
        final MutableInt succeededRealignmentCount = new MutableInt(0);

        final Map<GATKRead, GATKRead> mates = realignmentArgumentCollection.dontUseMates ? null
                : readsContext.getReadToMateMap(fragmentSize);

        for (final GATKRead read : readsContext) {
            if (passesFilter != Trilean.UNKNOWN) {
                break;
            } else if (read.getMappingQuality() < minMappingQuality || !supportsVariant(read, vc)) {
                continue;
            }

            // First check the read itself.  If it maps uniquely we don't need to bother with the mate
            final Realigner.RealignmentResult readRealignment = realigner.realign(read);
            if (readRealignment.isGood()) {
                succeededRealignmentCount.increment();
            } else if (mates != null && mates.containsKey(read)) {
                final GATKRead mate = mates.get(read);

                // now check the mate.  If it maps uniquely we consider the entire pair good
                final Realigner.RealignmentResult mateRealignment = realigner.realign(mate);
                if (mateRealignment.isGood()) {
                    succeededRealignmentCount.increment();
                } else {
                    // finally we check whether the pair maps uniquely even though the read and mate individually do not
                    final List<BwaMemAlignment> readRealignments = readRealignment.getRealignments();
                    final List<BwaMemAlignment> mateRealignments = mateRealignment.getRealignments();
                    final List<Pair<BwaMemAlignment, BwaMemAlignment>> plausiblePairs = readRealignments.stream()
                            .flatMap(r -> mateRealignments.stream()
                                    .filter(m -> m.getRefId() == r.getRefId() && Math.abs(m.getRefStart() - r.getRefStart()) < realignmentArgumentCollection.maxReasonableFragmentLength)
                                    .map(m -> new Pair<BwaMemAlignment, BwaMemAlignment>(r,m)))
                            .collect(Collectors.toList());

                    if (plausiblePairs.size() <= 1) {
                        succeededRealignmentCount.increment();
                    } else {
                        plausiblePairs.sort(Comparator.comparingInt(pair -> -pairScore(pair)) );
                        final int scoreDiff = pairScore(plausiblePairs.get(0)) - pairScore(plausiblePairs.get(1));
                        if (scoreDiff >= realignmentArgumentCollection.minAlignerScoreDifference) {
                            succeededRealignmentCount.increment();
                        } else {
                            failedRealignmentCount.increment();
                        }
                    }
                }
            } else {
                failedRealignmentCount.increment();
            }

            if (failedRealignmentCount.intValue() > maxFailedRealignments) {
                passesFilter = Trilean.FALSE;
            } else if (succeededRealignmentCount.intValue() >= sufficientGoodRealignments) {
                passesFilter = Trilean.TRUE;
            }
        }

        // if we haven't decided yet due to too few supporting reads, fail if there are more failures than successes
        passesFilter = passesFilter != Trilean.UNKNOWN ? passesFilter :
                Trilean.of(failedRealignmentCount.intValue() <= succeededRealignmentCount.intValue());

        vcfWriter.add(passesFilter == Trilean.TRUE ? vc : new VariantContextBuilder(vc).filter(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME).make());
    }

    private static int pairScore(final Pair<BwaMemAlignment, BwaMemAlignment> pair) {
        return pair.getFirst().getAlignerScore() + pair.getSecond().getAlignerScore();
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

    private boolean supportsVariant(final GATKRead read, final VariantContext vc) {
        final byte[] readBases = read.getBasesNoCopy();

        final int variantPositionInRead = ReadUtils.getReadCoordinateForReferenceCoordinate(read.getSoftStart(), read.getCigar(),
                vc.getStart(), ReadUtils.ClippingTail.RIGHT_TAIL, true);
        if ( variantPositionInRead == ReadUtils.CLIPPING_GOAL_NOT_REACHED || AlignmentUtils.isInsideDeletion(read.getCigar(), variantPositionInRead)) {
            return false;
        }

        for (final Allele allele : vc.getAlternateAlleles()) {
            final int referenceLength = vc.getReference().length();
            if (allele.length() == referenceLength) {   // SNP or MNP -- check whether read bases match variant allele bases
                if (allele.basesMatch(ArrayUtils.subarray(readBases, variantPositionInRead, Math.min(variantPositionInRead + allele.length(), readBases.length)))) {
                    return true;
                }
            } else {    // indel -- check if the read has the right CIGAR operator near position -- we don't require exact an exact match because indel representation is non-unique
                final boolean isDeletion = allele.length() < referenceLength;
                int readPosition = 0;    // offset by 1 because the variant position is one base before the first indel base
                for (final CigarElement cigarElement : read.getCigarElements()) {
                    if (Math.abs(readPosition - variantPositionInRead) <= indelStartTolerance) {
                        if (isDeletion && mightSupportDeletion(cigarElement) || !isDeletion && mightSupportInsertion(cigarElement)) {
                            return true;
                        }
                    } else {
                        readPosition += cigarElement.getLength();
                    }
                }
            }
        }
        return false;
    }

    private static boolean mightSupportDeletion(final CigarElement cigarElement) {
        final CigarOperator cigarOperator = cigarElement.getOperator();
        return cigarOperator == CigarOperator.D || cigarOperator == CigarOperator.S;
    }

    private final static boolean mightSupportInsertion(final CigarElement cigarElement) {
        final CigarOperator cigarOperator = cigarElement.getOperator();
        return cigarOperator == CigarOperator.I || cigarOperator == CigarOperator.S;
    }
}
