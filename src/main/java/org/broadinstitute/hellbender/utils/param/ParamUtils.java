package org.broadinstitute.hellbender.utils.param;

import org.apache.commons.math3.exception.NotFiniteNumberException;
import org.apache.commons.math3.util.MathUtils;

/**
 * This class should eventually be merged into Utils, which is in hellbender, and then this class should be deleted.
 *
 * Created by lichtens on 8/25/15.
 */
public class ParamUtils {

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     *
     * <p>Note that min can be greater than max, but that will guarantee that an exception is thrown.</p>
     *
     * @param val value to check
     * @param min minimum value for val
     * @param max maximum value for val
     * @param message the text message that would be pass to the exception thrown when {@code o == null}.
     * @return the same value
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static double inRange(final double val, final double min, final double max, final String message) {
        if ((val < min) || (val > max)){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is within range and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param min minimum value for val
     * @param max maximum value for val
     * @param message the text message that would be pass to the exception thrown when val gt min or val lt max.
     * @return the same value
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static long inRange(final long val, final double min, final double max, final String message) {
        if ((val < min) || (val > max)){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static long isPositiveOrZero(final long val, final String message) {
        if (val < 0){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is positive or zero and returns the same value or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static double isPositiveOrZero(final double val, final String message) {
        if (val < 0){
            throw new IllegalArgumentException(message);
        }
        return val;
    }

    /**
     * Checks that the  input is not infinity nor NaN or throws an {@link IllegalArgumentException}
     * @param val value to check
     * @param message the text message that would be pass to the exception thrown
     * @return the same value
     * @throws IllegalArgumentException if a {@code o == null}
     */
    public static double isFinite(final double val, final String message) {
        try {
            MathUtils.checkFinite(val);
        } catch (final NotFiniteNumberException ne) {
            throw new IllegalArgumentException(message);
        }
        return val;
    }
}
