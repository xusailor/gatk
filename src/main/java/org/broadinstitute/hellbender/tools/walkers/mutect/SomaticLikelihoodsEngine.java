package org.broadinstitute.hellbender.tools.walkers.mutect;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Created by David Benjamin on 3/9/17.
 */
public class SomaticLikelihoodsEngine {

    public static final double CONVERGENCE_THRESHOLD = 0.001;

    /**
     * Given a likelihoods matrix, calculate the parameters of the Dirichlet posterior distribution on their allele
     * fractions, which define a discrete distribution.
     * @param log10Likelihoods matrix of alleles x reads
     * @param priorPseudocounts
     */
    public static double[] alleleFractionsPosterior(final RealMatrix log10Likelihoods, final double[] priorPseudocounts) {
        final int numberOfAlleles = log10Likelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == priorPseudocounts.length, "Must have one pseudocount per allele.");

        double[] dirichletPosterior = new IndexRange(0, numberOfAlleles).mapToDouble(n -> 1.0);  // initialize flat posterior
        boolean converged = false;

        while(!converged) {
            // alleleCounts = \sum_r \bar{z}_r, where \bar{z}_r is an a-dimensional vector of the expectation of z_r with respect to q(f)
            final double[] alleleCounts = getEffectiveCounts(log10Likelihoods, dirichletPosterior);
            final double[] newDirichletPosterior = MathArrays.ebeAdd(alleleCounts, priorPseudocounts);
            converged = MathArrays.distance1(dirichletPosterior, newDirichletPosterior) < CONVERGENCE_THRESHOLD;
            dirichletPosterior = newDirichletPosterior;
        }

        return dirichletPosterior;
    }

    //same with flat prior
    public static double[] alleleFractionsPosterior(final RealMatrix log10Likelihoods) {
        final double[] flatPrior = new IndexRange(0, log10Likelihoods.getRowDimension()).mapToDouble(n -> 1);
        return alleleFractionsPosterior(log10Likelihoods, flatPrior);
    }


    /**
     * Given data log likelihoods and a Dirichlet prior for a categorical distribution, obtain the array of total
     * responsibilities for each category
     * @param log10Likelihoods
     * @param dirichletPrior
     * @return
     */
    @VisibleForTesting
    protected static double[] getEffectiveCounts(RealMatrix log10Likelihoods, double[] dirichletPrior) {
        final double[] effectiveLog10Weights = new Dirichlet(dirichletPrior).effectiveLog10MultinomialWeights();
        return MathUtils.sumArrayFunction(0, log10Likelihoods.getColumnDimension(),
                read -> MathUtils.posteriors(effectiveLog10Weights, log10Likelihoods.getColumn(read)));
    }


    /**
     * @param log10Likelihoods matrix of alleles x reads
     * @param priorPseudocounts
     */
    public static double log10Evidence(final RealMatrix log10Likelihoods, final double[] priorPseudocounts) {
        final int numberOfAlleles = log10Likelihoods.getRowDimension();
        Utils.validateArg(numberOfAlleles == priorPseudocounts.length, "Must have one pseudocount per allele.");
        final double[] alleleFractionsPosterior = alleleFractionsPosterior(log10Likelihoods, priorPseudocounts);
        final double priorContribution = log10DirichletNormalization(priorPseudocounts);
        final double posteriorContribution = -log10DirichletNormalization(alleleFractionsPosterior);

        final double[] log10AlleleFractions = new Dirichlet(alleleFractionsPosterior).effectiveLog10MultinomialWeights();

        final double likelihoodsAndEntropyContribution = new IndexRange(0, log10Likelihoods.getColumnDimension()).sum(r -> {
            final double[] log10LikelihoodsForRead = log10Likelihoods.getColumn(r);
            final double[] responsibilities = MathUtils.posteriors(log10AlleleFractions, log10LikelihoodsForRead);
            final double likelihoodsContribution = MathUtils.sum(MathArrays.ebeMultiply(log10LikelihoodsForRead, responsibilities));
            final double entropyContribution = Arrays.stream(responsibilities).map(SomaticLikelihoodsEngine::xLog10x).sum();
            return likelihoodsContribution - entropyContribution;

        });

        return priorContribution + posteriorContribution + likelihoodsAndEntropyContribution;
    }

    // same as above using the default flat prior
    public static double log10Evidence(final RealMatrix log10Likelihoods) {
        final double[] flatPrior = new IndexRange(0, log10Likelihoods.getRowDimension()).mapToDouble(n -> 1);
        return log10Evidence(log10Likelihoods, flatPrior);
    }

    private static double xLog10x(final double x) {
        return x < 1e-8 ? 0 : x * Math.log10(x);
    }

    public static double log10DirichletNormalization(final double[] dirichletParams) {
        final double logNumerator = Gamma.logGamma(MathUtils.sum(dirichletParams));
        final double logDenominator = MathUtils.sum(MathUtils.applyToArray(dirichletParams, Gamma::logGamma));
        return MathUtils.logToLog10(logNumerator - logDenominator);
    }

    // compute the likelihoods that each allele is contained at some allele fraction in the sample
    public static <A extends Allele> PerAlleleCollection<Double> somaticLog10Odds(final LikelihoodMatrix<A> log10Matrix) {
        final double log10EvidenceWithAllAlleles = log10Matrix.numberOfReads() == 0 ? 0 :
                log10Evidence(SomaticGenotypingEngine.getAsRealMatrix(log10Matrix));

        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int refIndex = SomaticGenotypingEngine.getRefIndex(log10Matrix);
        IntStream.range(0, log10Matrix.numberOfAlleles()).filter(a -> a != refIndex).forEach(a -> {
            final A allele = log10Matrix.getAllele(a);
            final LikelihoodMatrix<A> log10MatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(log10Matrix, allele);
            final double log10EvidenceWithoutThisAllele = log10MatrixWithoutThisAllele.numberOfReads() == 0 ? 0 :
                    log10Evidence(SomaticGenotypingEngine.getAsRealMatrix(log10MatrixWithoutThisAllele));
            lods.setAlt(allele, log10EvidenceWithAllAlleles - log10EvidenceWithoutThisAllele);
        });
        return lods;
    }
}
