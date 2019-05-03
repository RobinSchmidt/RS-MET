#include "OnePoleFilter.h"

//----------------------------------------------------------------------------
// construction/destruction:

OnePoleFilter::OnePoleFilter()
{
	setSampleRate(44100.0);  // sampleRate = 44100 Hz by default
	setMode      (0);        // bypass by default
	setCutoff    (20000.0);  // cutoff = 20000 Hz by default

	//calcCoeffs() is called by setMode and setCutoff

	resetBuffers();          //reset memorized samples to zero
}

OnePoleFilter::~OnePoleFilter()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void OnePoleFilter::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;
 sampleRateRec = 1.0 / sampleRate;

 calcCoeffs();
 return;
}

void OnePoleFilter::setMode(int newMode)
{
 mode = newMode; //0:bypass, 1:Low Pass, 2:High Pass
 calcCoeffs();
}

void OnePoleFilter::setCutoff(double newCutoff)
{
 if( (newCutoff>0.0) && (newCutoff<=20000.0) )
  cutoff = newCutoff;
	else
		cutoff = 20000.0;

 calcCoeffs();
 return;
}

void OnePoleFilter::setCoeffs(double newB0, double newB1, double newA1)
{
 b0 = newB0;
 b1 = newB1;
 a1 = newA1;
}

//----------------------------------------------------------------------------
//others:

void OnePoleFilter::calcCoeffs()
{
 //intermediate variable for calculation (x is the amount of decay
 //between adjacent samples):
 double x = exp( -2.0 * PI * cutoff * sampleRateRec);

 switch(mode)
 {
  case 1:  //calculate low-pass coefficients:
  {
   b0 = 1-x;
   b1 = 0.0;
   a1 = x;
  }break;

  case 2:  //calculate high-pass coefficients:
  {
   b0 =  0.5*(1+x);
   b1 = -0.5*(1+x);
   a1 = x;
  }break;

  default://bypass:
  {
   b0 = 1.0;
   b1 = 0.0;
   a1 = 0.0;
  }break;
 }//end of switch

}

void OnePoleFilter::resetBuffers()
{
 x_1 = 0.0;
 y_1 = 0.0;
}
