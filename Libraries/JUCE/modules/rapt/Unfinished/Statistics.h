#ifndef RAPT_STATISTICS_H
#define RAPT_STATISTICS_H

// wrap into class

/** Computes the sample cross-correlation sequence of 2 sequences x, y of length N and stores the
result in r which is also of length N. You may use this function also to compute an
autocorrelation sequence by passing the same pointer for x and y. The array r for the result must
be distinct from x and y. The sample cross-correlation is a biased estimate of the true
cross-corrleation. The bias is such that a value for lag k is underestimated by the factor
(N-k)/N. This puts lower weight on higher lags which is often desirable because these values are
less reliably estimated due to a smaller number of averaged terms. The function automatically
chooses between direct and FFT-based computation according to the sequence length N. It can
optionally remove the bias, i.e. call rsRemoveCorrelationBias after the computation. */
template<class T>
void rsCrossCorrelation(T x[], T y[], int N, T r[], bool removeBias = false);

/** Direct computation of cross-correlation - more efficient for short sequences.  */
template<class T>
void rsCrossCorrelationDirect(T x[], T y[], int N, T r[]);

/** FFT based computation of cross-correlation - more efficient for long sequences. */
template<class T>
void rsCrossCorrelationFFT(T x[], T y[], int N, T r[]);

/** FFT based computation of autocorrelation. For FFT-based computuation, the autocorrelation
saves some computations compared to the cross-correlation, that's why we have a separate
function.  */
template<class T>
void rsAutoCorrelationFFT(T x[], int N, T r[]);

/** Given a biased sample (cross- or auto-) correlation sequence of length N in x, this function
removes the bias of the estimate and stores the result in r (which may point to the same array
as x). @see crossCorrelation
*/
template<class T>
void rsRemoveCorrelationBias(T x[], int N, T r[]);

/** Given two arrays x and y of lengths Nx and Ny, this function computes the cross-correlation
between the arrays, defined as E{x*y} / sqrt(E{x^2}*E{y^2}), where E{...} denotes the expectation
value (which is approximated by the arithmetic mean in the actual computation). If the lengths
Nx, Ny are not equal, you can imagine that the effective length is taken to be that of the longer
array and the shorter array is conceptually zero-padded to match the length of the longer
(there's no literal zero-padding done, just the summation loop is shorter). */
template<class T>
T rsCrossCorrelation(T *x, int Nx, T *y, int Ny);

/** Given two arrays x and y of lengths Nx and Ny, this function computes the cross-correlation
between the arrays where the shorter one of the arrays is first stretched to the length of the
longer one by means of linear interpolation. This might be a good measure of shape similarity
of two arrays regardless of their overall lengths (i.e. a length-invariant shape-similarity
measure). An array with a certain shape and another array with a compressed or stretched version
of the same shape should have a "stretched" correlation value of 1 (up to error due to linear
interpolation). */
template<class T>
T rsStretchedCrossCorrelation(T *x, int Nx, T *y, int Ny);

#endif
