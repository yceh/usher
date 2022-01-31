/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   stat.h
 * CREATED ON:  23 August 2011, 16:45
 * AUTHOR:      Maciej F. Boni
 * 
 * DESCRIPTION: Implement statistical methods.
 * 
 * HISTORY:     Version    Date         Description
 *              1.0        2011-08-23   created
 * 
 * VERSION:     1.0
 * LAST EDIT:   23 August 2011
 * 
 * NOTE:        All the methods here are implemented with
 *              "long double" type for more precision.
 */

#ifndef STAT_H
#define STAT_H

namespace stat {

    namespace siegmund {

        long double continuousApprox(int m, int n, int k);

        long double discreteApprox(int m, int n, int k);

    }

    namespace correction {

        /**
         * Dunn-Sidak statistical correction for multiple comparisons
         * @param pVal              The P-value of the current sample.
         * @param numTotalSamples   Total number of samples.
         * @return  The statistically corrected P-value.
         * @note    If the given P-value is too small, the result will be
         *          the same as using Bonferroni correction.
         */
        long double dunnSidak(long double pVal, long double numTotalSamples);

        /**
         * Bonferroni statistical correction for multiple comparisons
         * @param pVal              The P-value of the current sample.
         * @param numTotalSamples   Total number of samples.
         * @return  The statistically corrected P-value.
         */
        long double bonferroni(long double pVal, long double numTotalSamples);
    }

}

#endif	/* STAT_H */
