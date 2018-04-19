package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.test.VariantContextTestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import spire.math.All;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.testng.Assert.*;

public class AssemblyBasedCallerGenotypingEngineUnitTest extends GATKBaseTest {

    @DataProvider(name = "getVcsAtThisLocation")
    public Object[][] getVcsAtThisLocationData() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{new ArrayList<>(), 1000, new ArrayList<>(), new ArrayList<>()});

        final Haplotype snpHaplotype = new Haplotype("ACTGGTCAACTGGTCAACTGGTCAACTGGTCA".getBytes());
        final List<Allele> snpAlleles = Arrays.asList(Allele.create("A", true), Allele.create("G"));
        final VariantContextBuilder snpVCBuilder = new VariantContextBuilder("a", "20", 1000, 1000, snpAlleles);
        final VariantContext snpVc = snpVCBuilder.make();
        snpHaplotype.setEventMap(new EventMap(Arrays.asList(snpVc)));

        // this one matches the snp haplotype above (to test duplicate removal)
        final Haplotype snpHaplotypeDuplicate = new Haplotype("ACTGGTCAACTGGTCAACTGGTCAACTGGACA".getBytes());
        final List<Allele> snpAlleles2 = Arrays.asList(Allele.create("A", true), Allele.create("G"));
        final VariantContextBuilder svpVC2Builder = new VariantContextBuilder("a", "20", 1000, 1000, snpAlleles2);
        final VariantContext snpVc2 = svpVC2Builder.make();
        final List<Allele> snpAlleles3 = Arrays.asList(Allele.create("T", true), Allele.create("A"));
        final VariantContextBuilder snpVC3Builder = new VariantContextBuilder("a", "20", 1020, 1020, snpAlleles3);
        final VariantContext snpVc3 = snpVC3Builder.make();
        snpHaplotypeDuplicate.setEventMap(new EventMap(Arrays.asList(snpVc2, snpVc3)));


        final Haplotype deletionHaplotype = new Haplotype("ACTGGTCAGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAlleles = Arrays.asList(Allele.create("ACTGGTCAACT", true), Allele.create("A"));
        final VariantContextBuilder deletionVCBuilder = new VariantContextBuilder("a", "20", 995, 1005, deletionAlleles);
        final VariantContext deletionVc = deletionVCBuilder.make();
        deletionHaplotype.setEventMap(new EventMap(Arrays.asList(deletionVc)));

        // matches the deletion alleles above but at a different position (to catch an edge case in duplicate removal)
        final Haplotype deletionHaplotypeFalseDuplicate = new Haplotype("ACTGGTCAGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAllelesFalseDuplicate = Arrays.asList(Allele.create("ACTGGTCAACT", true), Allele.create("A"));
        final VariantContextBuilder deletionFalseDuplicateBuilder = new VariantContextBuilder("a", "20", 998, 1008, deletionAllelesFalseDuplicate);
        final VariantContext deletionVcFalseDuplicate = deletionFalseDuplicateBuilder.make();
        deletionHaplotypeFalseDuplicate.setEventMap(new EventMap(Arrays.asList(deletionVcFalseDuplicate)));

        // doesn't overlap 1000
        final Haplotype deletionHaplotypeNoSpan = new Haplotype("CAACTGGTCAACTGGTCAACTGGTCAACTGGTCAACTGGTCA".getBytes());
        final List<Allele> deletionAllelesNoSpan = Arrays.asList(Allele.create("GTCAA", true), Allele.create("G"));
        final VariantContextBuilder deletionVcNoSpanBuilder = new VariantContextBuilder("a", "20", 990, 994, deletionAllelesNoSpan);
        final VariantContext deletionVcNoSpan = deletionVcNoSpanBuilder.make();
        deletionHaplotypeNoSpan.setEventMap(new EventMap(Arrays.asList(deletionVcNoSpan)));

        tests.add(new Object[]{Arrays.asList(snpHaplotype), 1000, new ArrayList<>(), Arrays.asList(snpVc)});
        tests.add(new Object[]{Arrays.asList(snpHaplotype, snpHaplotypeDuplicate), 1000, new ArrayList<>(), Arrays.asList(snpVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype), 995, new ArrayList<>(), Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype), 1000, new ArrayList<>(), Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype, deletionHaplotypeNoSpan), 1000, new ArrayList<>(), Arrays.asList(deletionVc)});
        tests.add(new Object[]{Arrays.asList(deletionHaplotype, deletionHaplotypeFalseDuplicate, deletionHaplotypeNoSpan), 1000, new ArrayList<>(), Arrays.asList(deletionVc, deletionVcFalseDuplicate)});

        tests.add(new Object[]{Arrays.asList(deletionHaplotype, snpHaplotype), 1000, new ArrayList<>(), Arrays.asList(deletionVc, snpVc)});

        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(snpVc), Arrays.asList(snpVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 995, Arrays.asList(deletionVc), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc, snpVc),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), snpVCBuilder.source("Comp1Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc, deletionVcNoSpan), Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make())});
        tests.add(new Object[]{new ArrayList<>(), 1000, Arrays.asList(deletionVc, deletionVcFalseDuplicate, deletionVcNoSpan),
                Arrays.asList(deletionVCBuilder.source("Comp0Allele0").make(), deletionFalseDuplicateBuilder.source("Comp1Allele0").make())});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "getVcsAtThisLocation")
    public void testGetVCsAtThisLocation(final List<Haplotype> haplotypes,
                                         final int loc,
                                         final List<VariantContext> activeAllelesToGenotype,
                                         final List<VariantContext> expectedVcsAtThisLocation) {

        final List<VariantContext> vcsAtThisPosition = AssemblyBasedCallerGenotypingEngine.getVCsAtThisLocation(haplotypes, loc, activeAllelesToGenotype);
        Assert.assertEquals(vcsAtThisPosition.size(), expectedVcsAtThisLocation.size());
        for (int i = 0; i < expectedVcsAtThisLocation.size(); i++) {
            VariantContextTestUtils.assertVariantContextsAreEqual(vcsAtThisPosition.get(i), expectedVcsAtThisLocation.get(i), new ArrayList<>());
            Assert.assertEquals(vcsAtThisPosition.get(i).getSource(), expectedVcsAtThisLocation.get(i).getSource());
        }
    }

}