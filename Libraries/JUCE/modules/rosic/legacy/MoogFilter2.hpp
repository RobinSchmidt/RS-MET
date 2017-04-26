/* 

MoogFilter2.h: interface for the MoogFilter2 class.

© Braindoc 2006 (www.braindoc.de)

This class implements a nonlinear model of the Moog-Ladder filter based on the paper
by Antti Huovilainen.

*/

#if !defined(MoogFilter2_h_Included)
#define MoogFilter2_h_Included

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "VstTools.h"
#include "AntiAliasFilter.h"


//const sample teeBeeFactor  = 33.0/18.0;

class MoogFilter2
{
public:
 //construction/destruction:
 MoogFilter2();
 virtual ~MoogFilter2();

 static const long maxOverSampling = 64;  //filter operates on a higher sample-rate than the system 
                                          //for closer approximation of a delay free feedback-loop

 //parameter settings:
 void   setSampleRate(sample SampleRate);
 void   setOverSamp  (long   OverSamp);
 void   setOutputTap (long   OutputTap);  //decide after which filter stage the output is take
                                          //range (1-4)
	void   setDrive     (sample Drive);      //filter input drive (in dB), input clips at ? 



 //parameters that are supposed to be updated at sample-rate (inlined access-functions):
 __forceinline void setCutoff  (sample Cutoff);
 __forceinline void setFeedback(sample Feedback);   //amount of inverted feedback from the output to
                                                    //the input for the resonant peak
 __forceinline void calcCoeffs ();                  //causes the filter to re-calculate the coeffiecients

 //audio processing:
 __forceinline sample getSample  (sample InSamp);

 //others:
 void   resetBuffers();

protected:
 //embedded objects:
 AntiAliasFilter antiAliasFlt;

	sampleAl sampleRate;
 sampleAl sampleRateRec;            //reciprocal of the sampleRate

 //parameter variables
 sampleAl cutoff;
 sampleAl feedback;
	sampleAl drive;                    //drive as amplitude value

 sampleAl V_t;
 sampleAl alpha, beta;

 long   outputTap;                //the number of the filter stage after which the output is taken

 sampleAl prevInSamp;               //previous input sample which is needed for interpolation
                                  //(->oversampling); the filter interpolates with a straight line
                                  //between the previous and the current input sample  

 sampleAl prevOutSamp;              //previous output sample which is fed back to the input

 long   overSampling;             //amount of oversampling (low oversampling rates cause the resonance
                                  //to degrade on high cutoff frequencies, because the phase-shift departs
                                  //from 180 degrees at high frequencies)

 sampleAl overSamplingRec;          //reciprocal

 sampleAl buffer[maxOverSampling]; //holds the straight line between the previous and the 
                                 //current input sample

 //buffering:
 sampleAl W[3];
 sampleAl W_1[3];
 sampleAl y[4];                   //the y[n] samples for the 4 stages (the current outputs)
 sampleAl y_1[4];                 //the y[n-1] samples for the 4 filter stages


 //the outputs of the 4 filter stages (taps):
 //sample tapOut[4];

	//protected methods:
 __forceinline sample calcOverSampledOutSamp(sample InSamp);

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
__forceinline void MoogFilter2::setCutoff(sample Cutoff)
{
 //restrict cutoff frequency to the range between 20 and 20000 Hz:
 if(Cutoff <= 20.0)
  cutoff = 20.0;
 else if(Cutoff >= 20000.0)
  cutoff = 20000.0;
 else cutoff = Cutoff;
}

__forceinline void MoogFilter2::setFeedback(sample Feedback)
{
 if(Feedback>=0)
  feedback = Feedback;
}

__forceinline void MoogFilter2::calcCoeffs()
{
 static sample g;
 g     = 1 - exp(-2*PI*cutoff*sampleRateRec*overSamplingRec);
 alpha = 2*V_t*g;
}

__forceinline sample MoogFilter2::getSample(sample In)
{
 static sampleAl outSamp, intSlope;
 static long   i;

 //fill the buffer with the linear interpolated segment between the previous and
 //the current input sample:

 //calculate the interpolation slope:
 intSlope = (In-prevInSamp)*overSamplingRec;

 //fill the buffer with the oversampled segment:
 for(i=0; i<overSampling; i++)
  buffer[i] = prevInSamp + intSlope*i;

 //filter the buffer:
 for(i=0; i<overSampling; i++)
  outSamp = calcOverSampledOutSamp(buffer[i]);

 //store the current input sample for the next call:
 prevInSamp  = In;
 //prevOutSamp = outSamp; //is already done in calcOverSampledOutSamp()

 return outSamp;
}


__forceinline sample MoogFilter2::calcOverSampledOutSamp(sample In)
{
 y[0] = y_1[0] + alpha*(tanh(beta*(drive*In-feedback*y_1[3]))-W_1[0]);
 W[0] = tanh(beta*y[0]);
 y[1] = y_1[1] + alpha*(W[0]-W_1[1]);
 W[1] = tanh(beta*y[1]);
 y[2] = y_1[2] + alpha*(W[1]-W_1[2]);
 W[2] = tanh(beta*y[2]);
 y[3] = y_1[3] + alpha*(W[2]-tanh(beta*y_1[3]));

 //buffering:
 y_1[0] = y[0];
 y_1[1] = y[1];
 y_1[2] = y[2];
 y_1[3] = y[3];
 W_1[0] = W[0];
 W_1[1] = W[1];
 W_1[2] = W[2];

 return antiAliasFlt.getSample(y[outputTap]);
}


#endif // !defined(MoogFilter2_h_Included)
