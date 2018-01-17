package org.broadinstitute.hellbender.tools.walkers;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.VariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RMSMappingQuality;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.GeneralPloidyFailOverAFCalculatorProvider;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.util.*;

/**
 * Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
 *
 * <p>
 * This tool is designed to perform joint genotyping on multiple samples pre-called with HaplotypeCaller to produce a
 * multi-sample callset in a highly scalable manner. However it can also be run on a single sample at a time to produce
 * a single-sample callset. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller
 * with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
 *
 * <h3>Input</h3>
 * <p>
 * One or more GVCFs produced by in HaplotypeCaller with the `-ERC GVCF` or `-ERC BP_RESOLUTION` settings, containing
 * the samples to joint-genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A final VCF in which all samples have been jointly genotyped.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Perform joint genotyping on a set of GVCFs enumerated in the command line</h4>
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" GenotypeGVCFs \
 *   -R reference.fasta \
 *   -V input1.g.vcf \
 *   -V input2.g.vcf \
 *   -V input3.g.vcf \
 *   -O output.vcf
 * </pre>
 *
 * <h4>Perform joint genotyping on a set of GVCFs listed in a text file, one per line</h4>
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" GenotypeGVCFs \
 *   -R reference.fasta \
 *   -V input_gvcfs.list \
 *   -O output.vcf
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>Only GVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call GVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
@CommandLineProgramProperties(summary = "Perform joint genotyping on one or more samples pre-called with HaplotypeCaller", oneLineSummary = "Perform joint genotyping on one or more samples pre-called with HaplotypeCaller", programGroup = VariantProgramGroup.class)
@DocumentedFeature
public final class MakeVQSRinput extends VariantWalker {

    public static final String PHASED_HOM_VAR_STRING = "1|1";
    public static final String ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME = "onlyOutputCallsStartingInIntervals";
    private static final String GVCF_BLOCK = "GVCFBlock";

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written", optional=false)
    private File outputFile;

    //@Argument(fullName="includeNonVariantSites", shortName="allSites", doc="Include loci found to be non-variant after genotyping", optional=true)
    //TODO This option is currently not supported.
    private boolean includeNonVariants = false;

    @ArgumentCollection
    private GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    @ArgumentCollection
    private final VariantAnnotationArgumentCollection variantAnnotationArgumentCollection = new VariantAnnotationArgumentCollection(
            Arrays.asList(StandardAnnotation.class.getSimpleName()),
            Collections.emptyList(),
            Collections.emptyList());

    /**
     * This option can only be activated if intervals are specified.
     */
    @Advanced
    @Argument(fullName= ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME,
            doc="Restrict variant output to sites that start within provided intervals",
            optional=true)
    private boolean onlyOutputCallsStartingInIntervals = false;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set
     * when appropriate. Note that dbSNP is not used in any way for the genotyping calculations themselves.
     */
    @ArgumentCollection
    private final DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private GenotypingEngine<?> genotypingEngine;
    // the annotation engine
    private VariantAnnotatorEngine annotationEngine;

    private ReferenceConfidenceVariantContextMerger merger;

    private VariantContextWriter vcfWriter;

    /** these are used when {@link #onlyOutputCallsStartingInIntervals) is true */
    private List<SimpleInterval> intervals;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        final VCFHeader inputVCFHeader = getHeaderForVariants();

        if(onlyOutputCallsStartingInIntervals) {
            if( !hasIntervals()) {
                throw new CommandLineException.MissingArgument("-L or -XL", "Intervals are required if --" + ONLY_OUTPUT_CALLS_STARTING_IN_INTERVALS_FULL_NAME + " was specified.");
            }
        }
        intervals = hasIntervals() ? intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()) :
                Collections.emptyList();

        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples()); //todo should this be getSampleNamesInOrder?

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(variantAnnotationArgumentCollection, dbsnp.dbsnp, Collections.emptyList());

        merger = new ReferenceConfidenceVariantContextMerger();

        setupVCFWriter(inputVCFHeader, samples);
    }

    private static boolean annotationShouldBeSkippedForHomRefSites(VariantAnnotation annotation) {
        return annotation instanceof RankSumTest || annotation instanceof RMSMappingQuality || annotation instanceof AS_RMSMappingQuality;
    }

    private void setupVCFWriter(VCFHeader inputVCFHeader, SampleList samples) {
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // add headers for annotations added by this tool
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));
        headerLines.add(GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        vcfWriter = createVCFWriter(outputFile);

        final Set<String> sampleNameSet = samples.asSetOfSamples();
        final VCFHeader vcfHeader = new VCFHeader(headerLines, new TreeSet<>(sampleNameSet));
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        ref.setWindow(10, 10); //TODO this matches the gatk3 behavior but may be unnecessary
        final VariantContext mergedVC = merger.merge(Collections.singletonList(variant), variant, includeNonVariants ? ref.getBase() : null, true, false);
        if ( !mergedVC.isVariant() || !GenotypeGVCFs.isProperlyPolymorphic(mergedVC) || mergedVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) == 0) {
            return;
        }

        ArrayList<Genotype> calledGenotypes = new ArrayList<>();
        for(final Genotype g : mergedVC.getGenotypes()) {
            final GenotypeBuilder gbuilder = new GenotypeBuilder(g);
            if (g.hasLikelihoods()) {
                GATKVariantContextUtils.makeGenotypeCall(g.getPloidy(), gbuilder, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, g.getLikelihoods().getAsVector(), mergedVC.getAlleles());
                calledGenotypes.add(gbuilder.make());
            }

        }

        //Add in annotations
        final VariantContext reannotated = annotationEngine.annotateContext(mergedVC, features, ref, null, a -> true);
        VariantContextBuilder builder = new VariantContextBuilder(RMSMappingQuality.getInstance().finalizeRawMQ(reannotated));
        builder.genotypes(calledGenotypes);


        if (!reannotated.hasAttribute(GATKVCFConstants.RAW_QUAL_APPROX_KEY))
            logger.warn("Variant is missing the QUALapprox key -- if this tool was run with GenomicsDB input, check the vidmap.json annotation info");
        final double QUALapprox = reannotated.getAttributeAsDouble(GATKVCFConstants.RAW_QUAL_APPROX_KEY, 0.0);
        if(QUALapprox < genotypeArgs.STANDARD_CONFIDENCE_FOR_CALLING - 10*Math.log10(genotypeArgs.snpHeterozygosity))  //we don't apply the prior to the QUAL approx, so do it here
            return;
        final int variantDP = reannotated.getAttributeAsInt(GATKVCFConstants.VARIANT_DEPTH_KEY, 0);

        double QD = QUALapprox / (double)variantDP;

        builder.attribute(GATKVCFConstants.QUAL_BY_DEPTH_KEY, QD);

        VariantContext result = builder.make();

        SimpleInterval variantStart = new SimpleInterval(result.getContig(), result.getStart(), result.getStart());
        if (!onlyOutputCallsStartingInIntervals || intervals.stream().anyMatch(interval -> interval.contains(variantStart))) {
            vcfWriter.add(result);
            //vcfWriter.add(new VariantContextBuilder(combinedVC).noGenotypes().make());
        }

    }

    @VisibleForTesting
    static boolean isSpanningDeletion(final Allele allele){
        return allele.equals(Allele.SPAN_DEL) || allele.equals(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED);
    }

    /**
     * Cleans up genotype-level annotations that need to be updated.
     * 1. move MIN_DP to DP if present
     * 2. propagate DP to AD if not present
     * 3. remove SB if present
     * 4. change the PGT value from "0|1" to "1|1" for homozygous variant genotypes
     * 5. move GQ to RGQ if the site is monomorphic
     *
     * @param vc            the VariantContext with the Genotypes to fix
     * @param createRefGTs  if true we will also create proper hom ref genotypes since we assume the site is monomorphic
     * @return a new set of Genotypes
     */
    @VisibleForTesting
    static List<Genotype> cleanupGenotypeAnnotations(final VariantContext vc, final boolean createRefGTs) {
        final GenotypesContext oldGTs = vc.getGenotypes();
        final List<Genotype> recoveredGs = new ArrayList<>(oldGTs.size());
        for ( final Genotype oldGT : oldGTs ) {
            final Map<String, Object> attrs = new HashMap<>(oldGT.getExtendedAttributes());

            final GenotypeBuilder builder = new GenotypeBuilder(oldGT);
            int depth = oldGT.hasDP() ? oldGT.getDP() : 0;

            // move the MIN_DP to DP
            if ( oldGT.hasExtendedAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY) ) {
                depth = parseInt(oldGT.getAnyAttribute(GATKVCFConstants.MIN_DP_FORMAT_KEY));
                builder.DP(depth);
                attrs.remove(GATKVCFConstants.MIN_DP_FORMAT_KEY);
            }

            attrs.remove(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);

            // update PGT for hom vars
            if ( oldGT.isHomVar() && oldGT.hasExtendedAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ) {
                attrs.put(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, PHASED_HOM_VAR_STRING);
            }

            // create AD if it's not there
            if ( !oldGT.hasAD() && vc.isVariant() ) {
                final int[] AD = new int[vc.getNAlleles()];
                AD[0] = depth;
                builder.AD(AD);
            }

            if ( createRefGTs ) {
                // move the GQ to RGQ
                if (oldGT.hasGQ()) {
                    builder.noGQ();
                    attrs.put(GATKVCFConstants.REFERENCE_GENOTYPE_QUALITY, oldGT.getGQ());
                }

                //keep 0 depth samples and 0 GQ samples as no-call
                if (depth > 0 && oldGT.hasGQ() && oldGT.getGQ() > 0) {
                    final List<Allele> refAlleles = Collections.nCopies(oldGT.getPloidy(), vc.getReference());
                    builder.alleles(refAlleles);
                }

                // also, the PLs are technically no longer usable
                builder.noPL();
            }

            recoveredGs.add(builder.noAttributes().attributes(attrs).make());
        }
        return recoveredGs;
    }

    private static int parseInt(Object attribute){
        if( attribute instanceof String) {
            return Integer.parseInt((String)attribute);
        } else if ( attribute instanceof Number){
            return ((Number) attribute).intValue();
        } else {
            throw new IllegalArgumentException("Expected a Number or a String but found something else.");
        }
    }

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        final UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = new GenotypeCalculationArgumentCollection(genotypeArgs);
        return uac;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null) {
            vcfWriter.close();
        }
    }
}

