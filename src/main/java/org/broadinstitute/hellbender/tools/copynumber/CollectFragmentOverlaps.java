package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.MappingQualityReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.HDF5SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.MetadataUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Collects fragment counts at specified intervals.  The count for each interval is calculated by counting
 * the number of fragments that overlap the interval.  The start and end positions of fragments are
 * inferred from read information; only properly paired, first-of-pair reads that pass the specified read filters are
 * used.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SAM format read data
 *     </li>
 *     <li>
 *         Intervals at which counts will be collected.
 *         The argument {@code interval-merging-rule} must be set to {@link IntervalMergingRule#OVERLAPPING_ONLY}
 *         and all other common arguments for interval padding or merging must be set to their defaults.
 *     </li>
 *     <li>
 *         Output file format.  This can be used to select TSV or HDF5 output.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Counts file.
 *         By default, the tool produces HDF5 format results. This can be changed with the {@code format} option
 *         to TSV format. Using HDF5 files with {@link CreateReadCountPanelOfNormals}
 *         can decrease runtime, by reducing time spent on IO, so this is the default output format.
 *         The HDF5 format contains information in the paths defined in {@link HDF5SimpleCountCollection}. HDF5 files may be viewed using
 *         <a href="https://support.hdfgroup.org/products/java/hdfview/">hdfview</a> or loaded in python using
 *         <a href="http://www.pytables.org/">PyTables</a> or <a href="http://www.h5py.org/">h5py</a>.
 *         The TSV format has a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in
 *         {@link SimpleCountCollection.SimpleCountTableColumn}, and the corresponding entry rows.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk CollectFragmentCounts \
 *          -I sample.bam \
 *          -L intervals.interval_list \
 *          --interval-merging-rule OVERLAPPING_ONLY \
 *          -O sample.counts.hdf5
 * </pre>
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Collects fragment counts at specified intervals",
        oneLineSummary = "Collects fragment counts at specified intervals",
        programGroup = CoverageAnalysisProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class CollectFragmentOverlaps extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(CollectFragmentOverlaps.class);

    private static final int DEFAULT_MINIMUM_MAPPING_QUALITY = 30;

    enum Format {
        TSV, HDF5
    }

    public static final String FORMAT_LONG_NAME = "format";
    public static final String MAXIMUM_FRAGMENT_LENGTH_LONG_NAME = "maximum-fragment-length";

    @Argument(
            doc = "Output file for fragment counts.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputCountsFile = null;

    @Argument(
            doc = "Output file format.",
            fullName = FORMAT_LONG_NAME,
            optional = true
    )
    private Format format = Format.HDF5;

    @Argument(
            doc = "Maximum fragment length.  Fragments with lengths greater than this will not be counted.",
            fullName = MAXIMUM_FRAGMENT_LENGTH_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int maximumFragmentLength = 10000;

    /**
     * Metadata contained in the BAM file.
     */
    private SampleLocatableMetadata metadata;

    /**
     * Overlap detector used to determine when fragments overlap with input intervals.
     */
    private OverlapDetector<SimpleInterval> intervalOverlapDetector;

    private Multiset<SimpleInterval> intervalMultiset = HashMultiset.create();

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>(super.getDefaultReadFilters());
        filters.add(ReadFilterLibrary.MAPPED);
        filters.add(ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT);
        filters.add(ReadFilterLibrary.NOT_DUPLICATE);
        filters.add(ReadFilterLibrary.FIRST_OF_PAIR); // this will make sure we don't double count
        filters.add(ReadFilterLibrary.PROPERLY_PAIRED);
        filters.add(new MappingQualityReadFilter(DEFAULT_MINIMUM_MAPPING_QUALITY));
        // this will only keep reads in pairs that are properly oriented and mapped on same chromosome
        // and lie within a few standard deviations from the mean of fragment size distributions
        return filters;
    }

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {
        metadata = MetadataUtils.fromHeader(getHeaderForReads(), Metadata.Type.SAMPLE_LOCATABLE);
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        //this check is currently redundant, since the master dictionary is taken from the reads;
        //however, if any other dictionary is added in the future, such a check should be performed
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(metadata.getSequenceDictionary(), sequenceDictionary)) {
            logger.warn("Sequence dictionary in BAM does not match the master sequence dictionary.");
        }

        CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);

        logger.info("Initializing and validating intervals...");
        final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);
        intervalOverlapDetector = OverlapDetector.create(intervals);

        //verify again that intervals do not overlap
        Utils.validateArg(intervals.stream().noneMatch(i -> intervalOverlapDetector.getOverlaps(i).size() > 1),
                "Input intervals may not be overlapping.");

        logger.info("Collecting fragment counts...");
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {
        //TODO collect information on reads that do not have a properly paired mate

        try {
            final SimpleInterval fragment = ReadOrientation.getFragment(read);
            logger.debug(fragment);
            if (fragment.getLengthOnReference() > maximumFragmentLength) {
                return;
            }
            final Set<SimpleInterval> overlappingIntervals = intervalOverlapDetector.getOverlaps(fragment);
            logger.debug(overlappingIntervals);

            //if fragment doesn't overlap any of the provided intervals, do nothing
            if (overlappingIntervals.isEmpty()) {
                return;
            }
            intervalMultiset.addAll(overlappingIntervals);
        } catch (final IllegalArgumentException e) {
            logger.warn(String.format("Exception encountered when calculating fragment count, skipping read: %s", read));
        }
    }

    @Override
    public Object onTraversalSuccess() {
        logger.info("Writing fragment counts to " + outputCountsFile);
        final SimpleCountCollection fragmentCounts = new SimpleCountCollection(
                metadata,
                intervalOverlapDetector.getAll().stream()
                        .map(i -> new SimpleCount(i, intervalMultiset.count(i)))
                        .collect(Collectors.toList()));

        if (format == Format.HDF5) {
            fragmentCounts.writeHDF5(outputCountsFile);
        } else {
            fragmentCounts.write(outputCountsFile);
        }

        return "SUCCESS";
    }

    /**
     * Helper class to calculate fragment corresponding to a properly paired read
     */
    @VisibleForTesting
    protected enum ReadOrientation {

        /**
         * Read was located on forward strand
         */
        FORWARD(read -> {
            logger.debug("FORWARD");
            return new SimpleInterval(
                    read.getContig(),
                    read.getUnclippedStart(),
                    read.getUnclippedStart() + read.getFragmentLength());
        }),

        /**
         * Read was located on reverse strand
         */
        REVERSE(read -> {
            logger.debug("REVERSE");
            return new SimpleInterval(
                    read.getContig(),
                    read.getUnclippedEnd()  + read.getFragmentLength(),
                    read.getUnclippedEnd());
        });

        private final Function<GATKRead, SimpleInterval> readToFragmentMapper;

        ReadOrientation(final Function<GATKRead, SimpleInterval> readToFragmentMapper) {
            this.readToFragmentMapper = readToFragmentMapper;
        }

        /**
         * Get a function that maps the read to a fragment.
         */
        protected Function<GATKRead, SimpleInterval> getReadToFragmentMapper() {
            return readToFragmentMapper;
        }

        /**
         * Get {@link ReadOrientation} instance corresponding to the orientation of the read.
         */
        protected static ReadOrientation getReadOrientation(final GATKRead read) {
            return read.getFragmentLength() > 0 ? FORWARD : REVERSE;
        }

        /**
         * Compute the fragment corresponding to the read.
         */
        protected static SimpleInterval getFragment(final GATKRead read) {
            return getReadOrientation(read).getReadToFragmentMapper().apply(read);
        }
    }
}
