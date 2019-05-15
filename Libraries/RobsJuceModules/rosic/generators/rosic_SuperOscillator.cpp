//#include "rosic_SuperOscillator.h"
//using namespace rosic;

SuperOscillator::SuperOscillator()
{
	numVoices     = 1;
	detune        = 0.0;
 freqSpacing   = 0;
 freqsAreDirty = true;
 phaseSpread   = 2.0/3.0;
	ampScale      = 1.0 / sqrt((double)numVoices);

	for(int i=0; i<maxVoices; i++)
	{
		increments1[i]  = 0.0;
		increments2[i]  = 0.0;
		phaseIndices[i] = 0.0;
	}

 detuneRatio = 0.5*(sqrt(5.0)-1.0);
 detuneDamp  = 1.0;

 sampleCount = 0;
}

SuperOscillator::~SuperOscillator()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void SuperOscillator::setSampleRate(double newSampleRate)
{
	if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

	sampleRateRec   = 1.0 / sampleRate;

 basicIncrement  = tableLengthDbl*freq*sampleRateRec;

	increments1[0] = basicIncrement*pulseFactor1;
	increments2[0] = basicIncrement*pulseFactor2;

 freqsAreDirty = true;
}

void SuperOscillator::setStartPhase(double newStartPhase)
{
 if( (newStartPhase>=0.0) && (newStartPhase<=360.0) )
  startIndex = (newStartPhase/360.0)*tableLengthDbl;
}

void SuperOscillator::setNumVoices(int newNumVoices)
{
	if( RAPT::rsIsOdd(newNumVoices) ) //check, if it is an odd number - only those are allowed
		numVoices = newNumVoices;
	ampScale   = 1.0 / sqrt((double)numVoices);
 freqsAreDirty = true;
}

void SuperOscillator::setFreqSpacing(int newFreqSpacing)
{
	if( newFreqSpacing >= 0 && newFreqSpacing <= 127)
  freqSpacing = newFreqSpacing;
 freqsAreDirty = true;
}

void SuperOscillator::setDetuneRatio(double newDetuneRatio)
{
 if( (newDetuneRatio>=0.0) && (newDetuneRatio<=1.0) )
  detuneRatio = newDetuneRatio;
 freqsAreDirty = true;
}

void SuperOscillator::setDetunePhaseSpread(double newDetunePhaseSpread)
{
 phaseSpread = newDetunePhaseSpread;
}

//-------------------------------------------------------------------------------------------------
// event processing:

void SuperOscillator::resetPhase()
{
	for(long i=0; i<numVoices; i++)
  phaseIndices[i] = startIndex + phaseSpread*(i*tableLengthDbl/numVoices);
 sampleCount = 0;
}

/*
Ideas:


-try other frequency ratios - maybe we should try to make them as irrational as possible, see here:
  https://www.youtube.com/watch?v=CaasbfdJdJg
  -maybe use other metallic ratios https://en.wikipedia.org/wiki/Metallic_mean
   ..but the golden ratio (maybe scaled by some simple rational ratio) is maximally irrational 
   number applies only to a pair of oscs - we want maximally irrational ratios between all pairs
   ...can we derive a formula (maybe based on continued fractions?) that somhow maximizes the 
   irrationality of all pairs in a set of numbers? if a/b = golden-ratio * simple-ratio, we have a
   single maximally irrational ratio - but in the case of 3 oscs, we want a/c and b/c also be 
   maximally irrational
   lat a = 1, b = 1.618.. - what should c be, such that a/c and b/c are also both(!) maximally
   irrational?
  -could it be that the metallic ratios are also mutually maximally irrational? -> figure out!
   ...what's the ratio between two metallic ratios, say the golden and the silver ratio? how does
   the continued fraction of goldRatio/silverRatio look like? ...how do we get the continued 
   fraction of a fraction of two continued fraction? is the some formula or algo for that? how 
   does the arithmetic of continued fractions work anyway? how are the added, multiplied and 
   divided? are there recursion formulas for the coefficients ci as functions of ai, bi when
   c = a/b?
   see here:
   https://math.stackexchange.com/questions/76036/arithmetic-of-continued-fractions-does-it-exist
   https://hack.org/mc/writings/hackerdom/hakmem/cf.html
   https://crypto.stanford.edu/pbc/notes/contfrac/compute.html
   https://github.com/mjdominus/cf
   https://github.com/mjdominus/cf/blob/master/cf_arith.c
  -metallic ratio rn for integer n is given by (n + sqrt(n^2+4))/2 - so the ratio of two
   metallic ratios for integers m, n is: rn/rm = (n + sqrt(n^2+4)) / (m + sqrt(m^2+4))
   ...maybe obtain the continued fraction expansion of these numerically and see, if we get small
   coefficients
  -or maybe just use continued fraction expansions that end with a string of 1s but may have some
   non-unity coeffs at the start
   [1,1,1,1,1,...], [2,1,1,1,1,...], [1,2,1,1,1,...], ...
   https://en.wikipedia.org/wiki/Continued_fraction
  -or how about sqrt of integers? try with sage things like: continued_fraction(sqrt(2)/sqrt(6))
   ...they seem to indeed have simple and sometimes small coeffs
   continued_fraction(sqrt(2)/sqrt(18)) is finite ...why? shouldn't it be finite only when the
   terms in the squra-roots are perfect squares (like 4, 9, 16, ..)? ..no! they seem to be finite,
   when it's twice a perfect square - why?
  -or: start with the two outer frequencies and compute at which time instant they happen to be 
   phase-aligned again and what their phase is and select the third frequency is such a way that it
   happens to be maximally out of phase with the other two at this instant - then compute, when the
   3 oscs are in phase sync again and select the 4th freq to be out of phase with them...and so on


 -maybe let client code set the ratios by passing a ratio-array

-to counteract the beating artifacts, we may somehow introduce a form of chaos:
 -take the output of the saw-stack
 -downshift the center frequency to 0 Hz by single-sideband modulation (the complex modulator sine 
  should  have the same frequency as the stack's center frequency)
 -modulate the frequencies of the saws with a signal obtained from that downshifted output (maybe
  apply a lowpass and/or bandpass, maybe clip it, maybe use 1/(1+x^2) waveshaper, etc.)
 -to do the SSB properly, we need to oversample by 2 - but we may use a simple low-quality 
  up/downsampler - the remaining aliasing may actually help for the chaoticity - or maybe don't 
  oversample at all just let the ssb signal alias as it may ...or maybe just use ring-modulation 
  with a sine...or maybe with the saw itself....whatever...all we want is to obtain a signal that 
  contains some sort of LFOish low-freq content - any sort of ringmodulation of the super-osc 
  output with a signal at the same frequency should produce some content around 0 Hz
-maybe bandpass the supersaw output at the center-freq and use the bandpassed version for ringmod:
 supersaw -> bandpass -> ringmod -> lowpass (should give erratic DC) -> frequency modulation
 -maybe we should modulate the increments in such a way that tends to spread the phases apart, when
  the envelope is high and drag them together when the envelope is low
 -or instead of the envelope, we use some other measure of total phase alignment
 

*/
