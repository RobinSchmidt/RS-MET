



/*

Background/Ideas:

make a class rsBlep (or rsBlepper, BlepApplicator).or rsStepBandLimiter
-double getSample(double in) produces an output sample at a time as usual
-void addStep(double time, double amplitude) may be called immediately before getSample to announce 
 that a step has occurred in between the previous call to getSample and the one to follow, "time"
 gives the fractional sample instant at which the step occurred (which determines at which phase we
 have to sample the corrector to be added) and "amplitude" gives the size of the step (which 
 determines the multiplication factor by which the corrector has to be scaled)
-the class should manage the required delay (i.e. delay the output by whatever number of samples is 
 needed)
-the user should be able to set the number of samples for the blep
-implementation:
 -use a delayline on N samples
  -setStep adds the corrector signal into the delayline (i.e. adds the scaled and shifted 
   step-corrector into the whole delayline
  -getSample adds to the incoming input a sample from the delayline at current index i, clears the 
   delay-element at the current index i and increments the index - the corrector signal as index i
   has been consumed by getSample
-maybe extend it to allow for supressing jumps in the derivative as well:
 -addDerivativeStep(double time, double amplitude)
 -see: http://dafx16.vutbr.cz/dafxpapers/18-DAFx-16_paper_33-PN.pdf
    http://research.spa.aalto.fi/publications/papers/dafx16-blamp/
-to fill the belp/blamp/blarabola/blolynomial buffer (the delayline), we need a function to 
 evaluate the an approximation of the residual (which is the difference between a heaviside-step 
 and a bandlimited step) at arbitrary times
pseudocode (likely wrong):
void addStep(double time, double amplitude)
{
  for(int i = 0; i < L; i++)  // L is the length of the delayline
    blepBuffer[wrap(bufIndex + i)] += amplitude  * blepResidual(time + i);
}
double getSample(double in)
{
  double y = delayBuffer[bufIndex] + blepBuffer[bufIndex];
  blepBuffer[bufIndex] = 0;            // clear blep at this position - it has been consumed
  delayBuffer[wrap(bufIndex+L)] = in;  // write input into delayline
  bufIndex = wrap(bufIndex + 1);
}
-maybe the blepResidual function could use a function oject rsBlepApproximator
-approximate the blep by writing a sampled sinc-function into an array (let the user select an 
 integer that says, how many samples should be taken per sinc-zero-crossing: 
 1: sample at zero-crossings, 2: at zero-crossings and maxima, etc.
 -a 2nd array is filled with the sinc-derivative values at the same instants (this derivative 
  should be evaluated analytically)
 -the sinc values together with the derivative values can be used for a hermite approximation
 -this piecewise polynomial approximation can be integrated to obtain the blep and integrated twice
  to obtain the blamp (alternatively, the blep can be calculated analytically)
 -maybe apply window to the sinc before differentiating/integrating ..or obtain the resisuals first
  from the unwidowed versions and the window the residual - try both, use whichever gives better 
  results
 -maybe, may compute the blep and blamp residual simultaneously (by evaluating a polynomial and its
  derivative simultaneously)


*/