#ifndef RS_TRANSFORMS_H
#define RS_TRANSFORMS_H

namespace RSLib
{

  /**
  FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)

  Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the time domain data in
  fftBuffer[0...2*fftFrameSize-1]. The FFT array takes and returns the cosine and sine parts in an
  interleaved manner, ie. fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. fftFrameSize
  must be a power of 2. It expects a complex input signal (see footnote 2), ie. when working with
  'common' audio signals our input signal has to be passed as {in[0],0.,in[1],0.,in[2],0.,...} asf.
  In that case, the transform of the frequencies of interest is in fftBuffer[0...fftFrameSize]. */
  template<class T>
  void smbFft(T *fftBuffer, long fftFrameSize, long sign);

  /** A direct implementation of the DFT with the slow O(N^2) complexitiy scaling. Its main utility 
  is to check the outputs of the FFT algorithms for debugging.
  \todo: check this function.
  */
  template<class T>
  void rsDFT(rsComplex<T> *buffer, int N);

  /** A general purpose FFT-routine for complex inputs. */
  template<class T>
  void rsFFT(rsComplex<T> *buffer, int N);

  /** Computes FFT-magnitudes and (optionally) also phases of a given real signal of length N. The
  resulting magnitude and phase arrays are both of length N/2 due to the symmetry of the Fourier
  transform of a real signal. */
  template<class T>
  void rsMagnitudeAndPhase(T *signal, int N, T *magnitudes, T *phases = NULL);

  /** Inverse transformation of rsFFT. */
  template<class T>
  void rsIFFT(rsComplex<T> *buffer, int N);

  /** A radix-2 FFT-routine. */
  template<class T>
  void rsRadix2FFT(rsComplex<T> *buffer, int N);

  /** An FFT-routine for arbitrary input sizes that uses the Bluestein algorithm. */
  //template<class T>
  //void rsBluesteinFFT(rsComplex<T> *buffer, int N);

}

#endif
