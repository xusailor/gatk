package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.annotations.Test;

import java.io.File;
import java.util.stream.StreamSupport;

import static org.testng.Assert.*;

public class FilterAlignmentArtifactsIntegrationTest extends CommandLineProgramTest {

    private static final String DREAM_BAMS_DIR = largeFileTestDir + "mutect/dream_synthetic_bams/";
    private static final String DREAM_VCFS_DIR = toolsTestDir + "mutect/dream/vcfs/";
    private static final String DREAM_MASKS_DIR = toolsTestDir + "mutect/dream/masks/";

    @Test
    public void testDreamTruthData() {
        final File tumorBam = new File(DREAM_BAMS_DIR, "tumor_1.bam");
        final File truthVcf = new File(DREAM_VCFS_DIR, "sample_1.vcf");
        final File mask = new File(DREAM_MASKS_DIR, "mask1.list");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-V", truthVcf.getAbsolutePath(),
                "-L", "20",
                "--bwa-mem-index-image", "/Users/davidben/Desktop/bwa_mem_hg_38/Homo_sapiens_assembly38.index_bundle",
                "-XL", mask.getAbsolutePath(),
                "-O", filteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        int j = 4;
    }

    @Test
    public void testFalseNegatives() {
        final File tumorBam = new File("/Volumes/cga_tcga-gsc/benchmark/data/realignments/synthetic.challenge.set4.tumor/synthetic.challenge.set4.tumor.bam");
        final File falseNegatives = new File("/Volumes/davidben/mutect/validations/test-realignment/dream/dream-4-false-negatives-from-realigner.vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-V", falseNegatives.getAbsolutePath(),
                "-L", "1", "-L", "2", "-L", "3", "-L", "4",
                "--bwa-mem-index-image", "/Users/davidben/Desktop/bwa_mem_hg_38/Homo_sapiens_assembly38.index_bundle",
                "-O", filteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        final long numFailingFilter = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                .filter(vc -> vc.getFilters().contains(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME))
                .count();
        final long numPassingFilter = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                .filter(vc -> !vc.getFilters().contains(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME))
                .count();
        int j = 10;
    }

    @Test
    public void testTrueNegatives() {
        final File tumorBam = new File("/Volumes/cga_tcga-gsc/benchmark/data/realignments/synthetic.challenge.set4.tumor/synthetic.challenge.set4.tumor.bam");
        final File trueNegatives = new File("/Volumes/davidben/mutect/validations/test-realignment/dream/dream-4-true-negatives-from-realigner.vcf");
        final File filteredVcf = createTempFile("filtered", ".vcf");

        final String[] args = {
                "-I", tumorBam.getAbsolutePath(),
                "-V", trueNegatives.getAbsolutePath(),
                "--bwa-mem-index-image", "/Users/davidben/Desktop/bwa_mem_hg_38/Homo_sapiens_assembly38.index_bundle",
                "-O", filteredVcf.getAbsolutePath()
        };

        runCommandLine(args);

        final long numFailingFilter = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                .filter(vc -> vc.getFilters().contains(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME))
                .count();
        final long numPassingFilter = StreamSupport.stream(new FeatureDataSource<VariantContext>(filteredVcf).spliterator(), false)
                .filter(vc -> !vc.getFilters().contains(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME))
                .count();
        int j = 10;
    }
}