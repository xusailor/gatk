package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.BwaMemIndexImageCreator;

public class RealignmentArgumentCollection {
    public static final int DEFAULT_MIN_SEED_LENGTH = 20;
    public static final double DEFAULT_DROP_RATIO = 0.5;
    public static final double DEFAULT_SEED_SPLIT_FACTOR = 1.2;
    public static final int DEFAULT_MAX_REASONABLE_FRAGMENT_LENGTH = 10000;

    /**
     * BWA-mem index image created by {@link BwaMemIndexImageCreator}
     */
    @Argument(fullName = "bwa-mem-index-image", shortName = "index", doc = "BWA-mem index image")
    public String bwaMemIndexImage;

    /**
     * Turn off the default mate-aware realignment
     */
    @Argument(fullName = "dont-use-mates", doc = "Realign individual reads without using their mates", optional = true)
    public boolean dontUseMates = false;

    /**
     * Maximum fragment length to be considered a reasonable pair alignment
     */
    @Argument(fullName = "max-reasonable-fragment-length", doc = "MAximum fragment length to be considered a reasonable pair alignment", optional = true)
    public int maxReasonableFragmentLength = DEFAULT_MAX_REASONABLE_FRAGMENT_LENGTH;

    /**
     * BWA-mem minimum seed length
     */
    @Argument(fullName = "minimum-seed-length", shortName = "min-seed-length", optional = true, doc = "Minimum number of matching bases to seed a MEM")
    public int minSeedLength = DEFAULT_MIN_SEED_LENGTH;

    /**
     * BWA-mem extension drop ratio
     */
    @Argument(fullName = "drop-ratio", shortName = "drop-ratio", doc = "Fraction of best MEM extension score below which other extensions are dropped")
    public double dropRatio = DEFAULT_DROP_RATIO;

    /**
     * BWA-mem seed split factor
     */
    @Argument(fullName = "seed-split-factor", shortName = "split-factor", doc = "MEMs longer than the minimum seed length times this factor are split and re-seeded.")
    public double splitFactor = DEFAULT_SEED_SPLIT_FACTOR;

}