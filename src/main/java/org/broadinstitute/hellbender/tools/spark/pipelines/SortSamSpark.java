package org.broadinstitute.hellbender.tools.spark.pipelines;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import scala.Tuple2;

import java.util.Collections;
import java.util.List;

@DocumentedFeature
@CommandLineProgramProperties(summary = "Sorts the input SAM/BAM/CRAM",
        oneLineSummary = "SortSam on Spark (works on SAM/BAM/CRAM)",
        programGroup = ReadDataManipulationProgramGroup.class)
@BetaFeature
public final class SortSamSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    public static final String SORT_ORDER_LONG_NAME = "sort-order";
    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc="the output file path", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String outputFile;

    @Argument(doc="the order to sort the file into", fullName = SORT_ORDER_LONG_NAME, optional = true)
    private SAMFileHeader.SortOrder sortOrder = SAMFileHeader.SortOrder.coordinate;

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    protected void onStartup() {
        if( sortOrder.getComparatorInstance() == null){
            throw new UserException.BadInput("Cannot sort a file in " + sortOrder + " order.  That ordering doesnt define a valid comparator.  "
                                                     + "Please choose a valid sort order");
        }
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> reads = getReads();
        int numReducers = getRecommendedNumReducers();
        logger.info("Using %s reducers", numReducers);

        final SAMFileHeader readsHeader = getHeaderForReads();
        ReadCoordinateComparator comparator = new ReadCoordinateComparator(readsHeader);
        JavaRDD<GATKRead> sortedReads;
        if (shardedOutput) {
            sortedReads = reads
                    .mapToPair(read -> new Tuple2<>(read, null))
                    .sortByKey(comparator, true, numReducers)
                    .keys();
        } else {
            sortedReads = reads; // sorting is done by writeReads below
        }
        readsHeader.setSortOrder(sortOrder);
        writeReads(ctx, outputFile, sortedReads);
    }
}
