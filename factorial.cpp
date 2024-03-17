#include "factorial.h"

/**
 * Calculate the factorial of a given integer.
 * 
 * param n The integer to calculate the factorial for
 * 
 * return The factorial of n
 */
double factorial(int n) {
    int result = 1;

    // Multiply every integer from n to 1
    for (int i = n; i > 0; i--) {
        result *= i;
    }
    return result;
}

/**
 *  Calculate the binomial coefficient of two given integers.
 * 
 * param m The first integer
 * param n The second integer
 * 
 * return The binomial coefficient of m and n
 */
double binomialCoef(int m, int n) {
    // Binommial coefficients are only defined for m >= n
    if (m < n) {
        return 0;
    }
    return factorial(m) / (factorial(n) * factorial(m - n));
}


/**
 * @brief Calculate the double factorial of a given integer.
 * 
 * @param n The integer to calculate the double factorial for
 * 
 * @return The double factorial of n
 */
double doubleFactorial(int n) {
    int result = 1;

    // Multiply every other integer from n to 1
    for (int i = n; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

