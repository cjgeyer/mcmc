
#include <R.h>
#include "mcmc.h"

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

void olbm(double *x, int *nin, int *pin, int *lin, double *mean,
    double *var, int *nocalcin)
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

