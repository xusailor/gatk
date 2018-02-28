package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import biz.k11i.xgboost.util.FVec;
import com.google.common.annotations.VisibleForTesting;

public class EvidenceFeatures  implements FVec {
    private static final long serialVersionUID = 1L;
    private final double[] values;

    EvidenceFeatures(final int numFeatures) {this.values = new double[numFeatures];}

    EvidenceFeatures(final double[] values) {this.values = values;}

    @Override
    public double fvalue(final int index){return values[index];}

    public int length() {return this.values.length;}

    public void set_value(final int index, final double value){values[index] = value;}

    @VisibleForTesting
    public double[] getValues() {return this.values;}
}
