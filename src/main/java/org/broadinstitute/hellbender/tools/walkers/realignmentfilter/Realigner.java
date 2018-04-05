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
    private final int minAlignerScoreDifference;
    private final int numberOfRegularContigs;

    public Realigner(final RealignmentArgumentCollection rfac) {
        final BwaMemIndex index = new BwaMemIndex(rfac.bwaMemIndexImage);
        maxReasonableFragmentLength = rfac.maxReasonableFragmentLength;
        minAlignerScoreDifference = rfac.minAlignerScoreDifference;
        numberOfRegularContigs = rfac.numRegularContigs;
        aligner = new BwaMemAligner(index);
        aligner.setMinSeedLengthOption(rfac.minSeedLength);
        aligner.setDropRatioOption((float) rfac.dropRatio);
        aligner.setSplitFactorOption((float) rfac.splitFactor);
        aligner.setFlagOption(BwaMemAligner.MEM_F_ALL);
    }

    public RealignmentResult realign(final GATKRead read) {
        final List<BwaMemAlignment> alignments = aligner.alignSeqs(Arrays.asList(read), GATKRead::getBasesNoCopy).get(0);

        final List<BwaMemAlignment> nonAltAlignments = alignments.size() == 1 ? alignments :
                alignments.stream().filter(a -> a.getRefId() < numberOfRegularContigs).collect(Collectors.toList());
        return checkAlignments(nonAltAlignments);
    }

    public RealignmentResult checkAlignments(final List<BwaMemAlignment> alignments) {
        if (alignments.isEmpty()) {
            return new RealignmentResult(false, Collections.emptyList());
        } else if (alignments.get(0).getRefId() < 0) {
            return new RealignmentResult(false, alignments);
        } else if (alignments.size() == 1) {
            return new RealignmentResult(true, alignments);
        } else {
            final int scoreDifference = alignments.get(0).getAlignerScore() - alignments.get(1).getAlignerScore();
            //TODO: check difference between best and second best aligner score
            return new RealignmentResult(scoreDifference >= minAlignerScoreDifference, alignments);
        }

        // TODO: we need to check that contig is the same or equivalent up to hg38 alt contig
        // TODO: do this by looking at index.getReferenceContigNames() and bamHeader.getSequenceDictionary().getSequences()
        // TODO: in IDE and seeing what the correspondence could be
        // TODO: put in check that start position is within eg 10 Mb of original mapping
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