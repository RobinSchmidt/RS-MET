
//-------------------------------------------------------------------------------------------------
// setup:

void ConvolverPartitioned::setImpulseResponse(double *newImpulseResponse, int newLength)
{
  if( newLength < 0 ) { RAPT::rsError("Length must be >= 0"); return; }

  M = newLength;
  directConvolver.setImpulseResponse(newImpulseResponse, RAPT::rsMin(newLength,
    directConvolutionLength));

  int accu             = directConvolutionLength;
  int currentLength    = directConvolutionLength;
  int numFftConvolvers = 0;
  while(newLength > accu)
  {
    numFftConvolvers += 1;
    accu             += currentLength;
    currentLength    *= 2;
  }
  fftConvolvers.resize(numFftConvolvers);

  int currentStart = directConvolutionLength;
  currentLength    = directConvolutionLength;
  for(int c=0; c<numFftConvolvers; c++)
  {
    if(c == numFftConvolvers-1)
    {
      // Last block might be shorter than currentLength, so we pass a zero-padded version:
      double* finalBlock  = new double[currentLength];
      int     finalLength = newLength-currentStart;    // length of non-zero part
      int k;
      for(k=0; k<finalLength; k++)
        finalBlock[k] = newImpulseResponse[currentStart+k];
      for(k=finalLength; k<currentLength; k++)
        finalBlock[k] = 0.0;
      fftConvolvers[c].setImpulseResponse(finalBlock, currentLength);
      delete[] finalBlock;
      // ToDo: try to avoid the memory allocation, maybe setImpulseResponse should take an 
      // optional parameter "zeroPadding"
    }
    else
      fftConvolvers[c].setImpulseResponse(&(newImpulseResponse[currentStart]), currentLength);

    currentStart  += currentLength;
    currentLength *= 2;
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void ConvolverPartitioned::clearImpulseResponse()
{
  directConvolver.clearImpulseResponse();
  for(size_t c = 0; c < fftConvolvers.size(); c++)
    fftConvolvers[c].clearImpulseResponse();
}

void ConvolverPartitioned::clearInputBuffers()
{
  directConvolver.clearInputBuffer();
  for(size_t c = 0; c < fftConvolvers.size(); c++)
    fftConvolvers[c].clearInputBuffer();
}


/*

ToDo:
-templatize and move to rapt
-experiment with using integer arithmetic to avoid roundoff noise - this requires using NTT instead
 of complex FFT...maybe that's faster 
 -is modular arithemtic is faster than complex? -> benchmark! certainly, if the modulus is apower 
  of 2 such the we may use bitmasking for the modulo operation...but does that work for NTT?
  ...or maybe we do not need to take the remainder after each multiplcation but only one at the
  end (or whenever there's an risk of overflow)
 -use the lower 23 bits of 32 bit integers, if the input comes from float32
-May avoid the bit-reversed ordering

Ideas:
-To balance the load of the FFT computations, "de-atomize" the FFT routine and compute at each 
 sample a few operations of all running FFTs, rather then just accumulating samples at each sample and 
 triggering full one-go FFTs at some samples
-To this end, introduce a class rsFourierTransformerBalanced that accepts samples and the input, 
 one at a time, does some FFT operations and returns. At the next sample, it does some more and 
 returns. After N samples have been accepted, the full FFT of the previous block should be available
 -it may be beneficial to replace the for-loops with while-loops - then we can directly jump into the 
  funtion and continue where we left off, for any value of h,k,j,jf,jl,Wjk ...which, i think, are
  needed for the status of the FFT, i.e. these should be passed as (struct of) reference parameters


Questions:
-using uinform block-sizes, input spectra can be re-used by just passing the input spectrum to the
 next block...but can this also be made to work, when the next block has twice the length...isn't 
 the DIF-FFT splitting the signal in the time-domain into two halves - that would mean, one half is 
 already computed and we only need to compute the other half

Resources:
http://www.cs.ust.hk/mjg_lib/bibs/DPSu/DPSu.Files/Ga95.PDF
https://www.researchgate.net/publication/280979094_Partitioned_convolution_algorithms_for_real-time_auralization

*/