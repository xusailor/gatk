package org.broadinstitute.hellbender.tools.walkers;

import com.intel.genomicsdb.GenomicsDBFeatureReader;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBConstants;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.testng.Assert.*;

/**
 * Created by gauthier on 3/9/18.
 */
public class MakeVQSRinputIntegrationTest extends CommandLineProgramTest {
    private static final VCFHeader VCF_HEADER = VariantContextTestUtils.getCompleteHeader();
    private static final String HG_00096 = largeFileTestDir + "gvcfs/HG00096.g.vcf.gz";
    private static final String HG_00268 = largeFileTestDir + "gvcfs/HG00268.g.vcf.gz";



    private static String getQueryJsonForGenomicsDB(String vidMappingFile, String callsetMappingFile, String tiledbWorkspace,
                                                    String referenceGenome, boolean produceGTField) throws IOException {
        //Produce temporary JSON query config file
        String indentString = "    ";
        String queryJSON = "{\n";
        queryJSON += indentString + "\"scan_full\": true,\n";
        queryJSON += indentString + "\"workspace\": \""+tiledbWorkspace+"\",\n";
        queryJSON += indentString + "\"array\": \""+GenomicsDBConstants.DEFAULT_ARRAY_NAME+"\",\n";
        queryJSON += indentString + "\"vid_mapping_file\": \""+vidMappingFile+"\",\n";
        queryJSON += indentString + "\"callset_mapping_file\": \""+callsetMappingFile+"\",\n";
        queryJSON += indentString + "\"produce_GT_field\": true,\n";
        queryJSON += indentString + "\"reference_genome\": \""+referenceGenome+"\"";
        queryJSON += "\n}\n";
        File tmpQueryJSONFile = File.createTempFile("queryJSON", ".json");
        tmpQueryJSONFile.deleteOnExit();
        FileWriter fptr = new FileWriter(tmpQueryJSONFile);
        fptr.write(queryJSON);
        fptr.close();
        return tmpQueryJSONFile.getAbsolutePath();
    }
    //Produce temporary JSON query config file

    private static GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> getGenomicsDBFeatureReader(final String workspace, final String reference, boolean produceGTField) throws IOException {
        if (produceGTField) {
            return new GenomicsDBFeatureReader<>(
                    "",
                    getQueryJsonForGenomicsDB(new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).getAbsolutePath(),
                            new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).getAbsolutePath(),
                            workspace,
                            reference,
                            produceGTField),
                    new BCF2Codec());
        } else {
            return new GenomicsDBFeatureReader<>(
                    new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).getAbsolutePath(),
                    new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).getAbsolutePath(),
                    workspace,
                    GenomicsDBConstants.DEFAULT_ARRAY_NAME,
                    reference,
                    new File(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME).getAbsolutePath(),
                    new BCF2Codec());
        }
    }

    private void writeToGenomicsDB(final List<String> vcfInputs, final SimpleInterval interval, final String workspace,
                                   final int batchSize, final Boolean useBufferSize, final int bufferSizePerSample, int threads) {
        final ArgumentsBuilder args = new ArgumentsBuilder();
        args.addArgument(GenomicsDBImport.WORKSPACE_ARG_LONG_NAME, workspace);
        args.addArgument("L", IntervalUtils.locatableToString(interval));
        vcfInputs.forEach(vcf -> args.addArgument("V", vcf));
        args.addArgument("batch-size", String.valueOf(batchSize));
        args.addArgument(GenomicsDBImport.VCF_INITIALIZER_THREADS_LONG_NAME, String.valueOf(threads));
        if (useBufferSize)
            args.addArgument("genomicsdb-vcf-buffer-size", String.valueOf(bufferSizePerSample));

        runCommandLine(args);
    }

    private static void checkJSONFilesAreWritten(final String workspace) {
        Assert.assertTrue(new File(workspace, GenomicsDBConstants.DEFAULT_VIDMAP_FILE_NAME).exists());
        Assert.assertTrue(new File(workspace, GenomicsDBConstants.DEFAULT_CALLSETMAP_FILE_NAME).exists());
        Assert.assertTrue(new File(workspace, GenomicsDBConstants.DEFAULT_VCFHEADER_FILE_NAME).exists());
    }

    private static void checkGenomicsDBAgainstExpected(final String workspace, final SimpleInterval interval, final String expectedCombinedVCF, final String referenceFile, final boolean testAll) throws IOException {
        final GenomicsDBFeatureReader<VariantContext, PositionalBufferedStream> genomicsDBFeatureReader =
                getGenomicsDBFeatureReader(workspace, referenceFile, !testAll);

        final AbstractFeatureReader<VariantContext, LineIterator> combinedVCFReader =
                AbstractFeatureReader.getFeatureReader(expectedCombinedVCF, new VCFCodec(), true);

        try (CloseableTribbleIterator<VariantContext> actualVcs =
                     genomicsDBFeatureReader.query(interval.getContig(), interval.getStart(), interval.getEnd());

             CloseableTribbleIterator<VariantContext> expectedVcs =
                     combinedVCFReader.query(interval.getContig(), interval.getStart(), interval.getEnd())) {

            BaseTest.assertCondition(actualVcs, expectedVcs, (a, e) -> {
                // Test that the VCs match
                if (testAll) {
                    // To correct a discrepancy between genotypeGVCFs which outputs empty genotypes as "./." and GenomicsDB
                    // which returns them as "." we simply remap the empty ones to be consistent for comparison
                    List<Genotype> genotypes = a.getGenotypes().stream()
                            .map(g -> g.getGenotypeString().equals(".")?new GenotypeBuilder(g).alleles(GATKVariantContextUtils.noCallAlleles(2)).make():g)
                            .collect(Collectors.toList());
                    a = new VariantContextBuilder(a).genotypes(genotypes).make();
                    VariantContextTestUtils.assertVariantContextsAreEqualAlleleOrderIndependent(a, e, Collections.emptyList(), VCF_HEADER);

                    // Test only that the genotypes match
                } else {
                    List<Genotype> genotypes = e.getGenotypes().stream()
                            .map(g -> g.getGenotypeString().equals(".")?new GenotypeBuilder(g).alleles(Collections.emptyList()).make():g)
                            .collect(Collectors.toList());
                    e = new VariantContextBuilder(e).genotypes(genotypes).make();
                    VariantContextTestUtils.assertVariantContextsHaveSameGenotypes(a, e);
                }
            });
        }
    }

    @DataProvider(name="VCFdata")
    public Object[][] getVCFdata() {
        return new Object[][] {
                new Object[]{1, },
                new Object[]{2, }
        };
    }


    @Test (dataProvider = "VCFdata")
    public void testUsingGenomicsDB() throws IOException {
        final String workspace = createTempDir("genomicsdb-tests-").getAbsolutePath() + "/workspace";

        writeToGenomicsDB(vcfInputs, interval, workspace, 0, false, 0, 1);
        checkJSONFilesAreWritten(workspace);
        File expectedCombinedVCF = runCombineGVCFs(vcfInputs, interval, b38_reference_20_21, CombineGVCFArgs);
        checkGenomicsDBAgainstExpected(workspace, interval, expectedCombinedVCF.getAbsolutePath(), b38_reference_20_21, true);
    }

}