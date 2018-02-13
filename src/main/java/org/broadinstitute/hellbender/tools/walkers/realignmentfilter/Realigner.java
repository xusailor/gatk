package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;


import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class Realigner {
    private final BwaMemAligner aligner;
    private final int maxReasonableFragmentLength;

    public Realigner(final RealignmentArgumentCollection rfac) {
        final BwaMemIndex index = new BwaMemIndex(rfac.bwaMemIndexImage);
        maxReasonableFragmentLength = rfac.maxReasonableFragmentLength;
        aligner = new BwaMemAligner(index);
        aligner.setMinSeedLengthOption(rfac.minSeedLength);
        aligner.setDropRatioOption((float) rfac.dropRatio);
        aligner.setSplitFactorOption((float) rfac.splitFactor);
        aligner.setFlagOption(BwaMemAligner.MEM_F_ALL);
    }

    public RealignmentResult realign(final GATKRead read) {
        final List<BwaMemAlignment> alignments = aligner.alignSeqs(Arrays.asList(read), GATKRead::getBasesNoCopy).get(0);
        return checkAlignments(alignments);
    }

    public RealignmentResult checkAlignments(final List<BwaMemAlignment> alignments) {
        if (alignments.isEmpty()) {
            return new RealignmentResult(false, Collections.emptyList());
        } else if (alignments.get(0).getRefId() < 0) {
            return new RealignmentResult(false, alignments);
        }


        // TODO: THIS IS COMPLETELY WRONG
        if (alignments.size() > 1) {
            //TODO: check difference between best and second best aligner score
            return new RealignmentResult(false, alignments);
        }


        // TODO: perhaps check number of mismatches in second best alignment?

        // TODO: we need to check that contig is the same or equivalent up to hg38 alt contig
        // TODO: do this by looking at index.getReferenceContigNames() and bamHeader.getSequenceDictionary().getSequences()
        // TODO: in IDE and seeing what the correspondence could be
        // TODO: put in check that start position is within eg 10 Mb of original mapping

        return new RealignmentResult(true, alignments);
    }

    public static class RealignmentResult {
        private final boolean mapsToSupposedLocation;
        private final List<BwaMemAlignment> realignments;

        public RealignmentResult(boolean mapsToSupposedLocation, List<BwaMemAlignment> realignments) {
            this.mapsToSupposedLocation = mapsToSupposedLocation;
            this.realignments = realignments;
        }

        public boolean isGood() { return mapsToSupposedLocation;  }

        public List<BwaMemAlignment> getRealignments() { return realignments; }
    }
}