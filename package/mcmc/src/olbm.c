
/*
*
* mcmc and MCMC package for R
* Copyright (c) 2005 Charles J. Geyer
*
* All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, and/or sell copies of the
* Software, and to permit persons to whom the Software is furnished to do so,
* provided that the above copyright notice(s) and this permission notice appear
* in all copies of the Software and that both the above copyright notice(s) and
* this permission notice appear in supporting documentation.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY RIGHTS.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS NOTICE
* BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL DAMAGES,
* OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
* WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
* ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*
* Except as contained in this notice, the name of a copyright holder shall
* not be used in advertising or otherwise to promote the sale, use or other
* dealings in this Software without prior written authorization of the
* copyright holder.
*/

#include <R.h>

/* overlapping batch means for vector time series
*
*  input:
*
*    x       time series, n x p matrix, n time points and p components
*    len     batch length
*
*  output:
*
*    mean    sample mean, a p vector
*    var     estimated variance of sample mean, a p x p matrix
*
*/

#define X(I,J)    	x[(I) + n * (J)]
#define VAR(I,J)    var[(I) + p * (J)]

void olbm(double *x, long *nin, long *pin, long *lin, double *mean,
    double *var, long *nocalcin)
{
    int n = nin[0];
    int p = pin[0];
    int len = lin[0];
    double nbatch = n - len + 1;
    int nocalc = nocalcin[0];
    double *work = (double *) R_alloc(p, sizeof(double));

    int i, j, k, l;

    if (len > n)
    	error("len > n\n");

    if (! nocalc)
    	for (i=0; i<p; i++) {
    		double sum = 0.0;
    		for (k=0; k<n; k++)
    			sum += X(k,i);
    		mean[i] = sum / n;
    	}

    /* easier to work with len * means, change means to that until
     * further notice
     */
    for (i=0; i<p; i++)
    	mean[i] *= len;

    for (i=0; i<p; i++) {
    	work[i] = 0.0;
    	for (k=0; k<len; k++)
    		work[i] += X(k,i);
    	for (j=i; j>=0; j--)
    		VAR(i,j) = (work[i] - mean[i]) * (work[j] - mean[j]);
    }

    for (k=0, l=len; l<n; k++, l++)
    	for (i=0; i<p; i++) {
    		work[i] -= X(k,i);
    		work[i] += X(l,i);
    		for (j=i; j>=0; j--)
    			VAR(i,j) += (work[i] - mean[i]) *
    				(work[j] - mean[j]);
    	}

    /* fix up means and variances, divide out factors of len and len^2 */
    for (i=0; i<p; i++)
    	mean[i] /= len;

    for (i=0; i<p; i++)
    	for (j=0; j<=i; j++) {
    		VAR(i,j) /= nbatch * n * len;
    		if (j < i) VAR(j,i) = VAR(i,j);
    	}

}

