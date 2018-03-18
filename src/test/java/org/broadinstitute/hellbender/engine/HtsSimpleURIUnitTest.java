package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class HtsSimpleURIUnitTest {

    @DataProvider(name = "validURIs")
    public Object[][] getValidURIs() {
        return new Object[][] {
                {"file:///path/to/file.bam"},
                {"file:///path.bam"},
                {"gs://abucket/bucket",},
                {"gs://hellbender/test/resources/"},
                {"gs://path.bam"},
                {"hdfs://namenode/path/to/file"},
                {"jimfs:afile"},

                // has a "#" in the path part, which is valid URI syntax, but without encoding will be
                // treated as a fragment delimiter; if URLEncoded, it will no longer be a valid (file/path) name
                {"/project/gvcf-pcr/23232_1#1/1.g.vcf.gz"},

                // TODO: should these be treated as having a "file" scheme ?
                {"/path/to/file.bam"},
                {"path/to/file.bam"},
                {"localFile"},

                {FeatureDataSource.GENOMIC_DB_URI_SCHEME + "adb"},
                {CloudStorageFileSystem.URI_SCHEME + "//abucket/bucket"},
                {CloudStorageFileSystem.GCS_VIEW + "//abucket/bucket"},
        };
    }

    @DataProvider(name = "invalidURIs")
    public Object[][] getInvalidURIs() {
        return new Object[][] {
                {"file://^"},
                {"file://"},
        };
    }

    @DataProvider(name = "validPaths")
    public Object[][] getValidPaths() {
        return new Object[][]{
                {"/project/gvcf-pcr/23232_1#1/1.g.vcf.gz"},
                {"file:///path/to/file.bam"},
                {"file:///path.bam"},
                {"gs://hellbender/test/resources/"},
                {"gs://path.bam"},
                {"/path/to/file.bam"},
                {"path/to/file.bam"},
                {"localFile"},
                {CloudStorageFileSystem.URI_SCHEME + "//abucket/bucket"},
                {CloudStorageFileSystem.GCS_VIEW + "//abucket/bucket"},
                {FeatureDataSource.GENOMIC_DB_URI_SCHEME + "adb"},
        };
    }

    @DataProvider(name = "invalidPaths")
    public Object[][] getInvalidPaths() {
        return new Object[][]{
                // valid URIs that are not valid as a path

                // unknown namenode
                {"hdfs://nonexistentnamenode/path/to/file"},

                // 2 slashes (vs. 3) causes "path" to be interpreted as an invalid authority name
                {"file://path/to/file.bam"},

                // TODO: this should result in an unknown provider exception, but the current implementation as taken
                // TODO: from IOUtils accepts any unknown provider and creates a new FileSystem for it ?
                //{"alice://foobar"},
        };
    }

    @Test(dataProvider = "validURIs")
    public void testValidURIs(final String validURIString) {
        final HtsURI htsURI = new HtsSimpleURI(validURIString);
        Assert.assertNotNull(htsURI);
        Assert.assertEquals(htsURI.getURI().toString(), validURIString);
    }

    @Test(dataProvider = "invalidURIs", expectedExceptions = IllegalArgumentException.class)
    public void testInvalidURIs(final String invalidURIString) {
        new HtsSimpleURI(invalidURIString);
    }

    @Test(dataProvider = "validPaths")
    public void testIsPathValidPaths(final String invalidPathString) {
        final HtsURI htsURI = new HtsSimpleURI(invalidPathString);
        Assert.assertTrue(htsURI.isPath());
        Assert.assertEquals(htsURI.getURI().toString(), htsURI.getURIString());
    }

    @Test(dataProvider = "invalidPaths")
    public void testIsPathInvalidPaths(final String invalidPathString) {
        final HtsURI htsURI = new HtsSimpleURI(invalidPathString);
        Assert.assertFalse(htsURI.isPath());
    }

    @Test(dataProvider = "validPaths")
    public void testToPathValidPaths(final String validPathString) {
        final HtsURI htsURI = new HtsSimpleURI(validPathString);
        Assert.assertTrue(htsURI.isPath());
        Assert.assertNotNull(htsURI.toPath());
    }

    @Test(dataProvider = "invalidPaths", expectedExceptions = IllegalArgumentException.class)
    public void testToPathInvalidPaths(final String invalidPathString) {
        final HtsURI htsURI = new HtsSimpleURI(invalidPathString);
        Assert.assertFalse(htsURI.isPath());
        final String reason = htsURI.getToPathFailureReason();
        //System.out.println(reason);
        htsURI.toPath();
    }

    @Test(dataProvider = "validPaths")
    public void testEncodedValidPaths(final String invalidPathString) {
        final HtsURI htsURI = new HtsSimpleURI(invalidPathString);
        Assert.assertTrue(htsURI.isPath());
        Assert.assertNotNull(htsURI.toPath());
    }


    // IOUtils tests for comparison
    @Test(dataProvider = "validPaths")
    public void testIOUtilsValidPaths(final String validPathString) {
        Assert.assertNotNull(IOUtils.getPath(validPathString));
    }

    @Test(dataProvider = "invalidPaths", expectedExceptions = IllegalArgumentException.class)
    public void testIOUtilsInvalidPaths(final String invalidPathString) {
        Assert.assertNotNull(IOUtils.getPath(invalidPathString));
    }

}
