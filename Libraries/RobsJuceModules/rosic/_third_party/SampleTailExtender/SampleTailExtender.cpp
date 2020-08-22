//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#define _USE_MATH_DEFINES
#include "SampleTailExtender.h"
#include <assert.h>
#include <cfloat>
#include "FFT.h"
#include <algorithm>

//=================================================================================
SampleTailExtender::SampleTailExtender()
{
    audioAnalysisFrame.resize (audioFrameSize);
}

//=================================================================================
std::vector<double> SampleTailExtender::extendSample (std::vector<double> inputSignal, double sampleRate, int notePitchClass)
{
    //// please use audio signals that are at least 2 seconds in length
    //assert (inputSignal.size() > (int) (sampleRate * 2));
    
    int synthesisStartPointInSamples = (int) (synthesisStartPointInSeconds * sampleRate);
    synthesisStartPointInSamples = std::min ((int)inputSignal.size() - 1, synthesisStartPointInSamples);
    
    int splicePoint = synthesisStartPointInSamples - crossfadeLength;
    
    // fill up the audio analysis frame
    int n = 0;
    for (int i = splicePoint - (audioFrameSize / 2); i < splicePoint + (audioFrameSize / 2); i++)
    {
        audioAnalysisFrame[n] = inputSignal[i];
        n++;
    }
    
    // calcualte the FFT
    FFT fft (audioFrameSize);
    fft.performFFT (audioAnalysisFrame);
    const std::vector<double>& magnitudeSpectrum = fft.getMagnitude();
    const std::vector<double>& phaseSpectrum = fft.getPhase();
    
    // calculate the indices of all harmonics
    std::vector<int> harmonics = getIndicesOfHarmonics (magnitudeSpectrum, sampleRate);
    
    // analyse the major harmonics of the signal for their energy profiles and establish
    // the amount of beating
    harmonicAnalyser.analyseHarmonics (inputSignal, sampleRate, splicePoint, audioFrameSize, notePitchClass);
    
    // calculate decay envelope
    double splicePointPeak = calculateSplicePointPeakEnergy (inputSignal, splicePoint);
    std::vector<double> decayEnvelope = makeDecayEnvelope (splicePointPeak, sampleRate, decayRate);
    
    int numSynthesisSamples = (int)decayEnvelope.size();
    std::vector<double> synthesisedSignal (numSynthesisSamples, 0);
    
    for (auto& h : harmonics)
    {
        double mag = getMagnitudeForSinusoid (magnitudeSpectrum, h);
        double phase = getPhaseForSinusoid (phaseSpectrum, h);
        double magInDecibels = 20. * log10 (mag);
        
        if (magInDecibels > -100.)
        {
            double frequency = getTunedFrequency (magnitudeSpectrum, h, sampleRate);
            std::vector<double> sinusoid = synthesiseSinusoid (numSynthesisSamples, 1.0, phase, frequency, sampleRate);
            
            // if there is harmonic beating to be applied
            if (h < (int)harmonicAnalyser.beatingAmount.size() && harmonicAnalyser.beatingAmount[h] > 0.0)
            {
                applyHarmonicBeatingToSinusoid (sinusoid, sampleRate, beatingStrength * harmonicAnalyser.beatingAmount[h], harmonicAnalyser.beatingFrequency[h]);
            }
            
            for (size_t i = 0; i < synthesisedSignal.size(); i++)
            {
                synthesisedSignal[i] += sinusoid[i] * (decayEnvelope[i] * mag);
            }
        }
    }
    
    // clip the signal when it starts to produce bit noise because it is too quiet
    setSignalToZeroBelowThreshold (synthesisedSignal, cutoffThresholdInDecibels);
    
    // synthesise new signal
    std::vector<double> outputSignal (synthesisStartPointInSamples, 0);
    
    for (size_t i = 0; i < outputSignal.size(); i++)
        outputSignal[i] = inputSignal[i];

    // do crossfade
    for (int i = 0; i < crossfadeLength; i++)
    {
        double alpha = (double)i / (double)crossfadeLength;
        double fadeInEnvelope = cos ( (1. - alpha) * M_PI) * 0.5 + 0.5;
        double fadeOutEnvelope = cos (alpha * M_PI) * 0.5 + 0.5;

		if (splicePoint + i > outputSignal.size()-1 || synthesisedSignal.size() == 0 || i > synthesisedSignal.size()-1)
			break;

        outputSignal[splicePoint + i] = fadeOutEnvelope * outputSignal[splicePoint + i] + fadeInEnvelope * synthesisedSignal[i];
    }

    // add synthesised signal to the end
    for (size_t i = crossfadeLength; i < synthesisedSignal.size(); i++)
        outputSignal.push_back (synthesisedSignal[i]);
    
    return outputSignal;
}


//=================================================================================
void SampleTailExtender::setSynthesisStartPointInSeconds (double seconds)
{
    synthesisStartPointInSeconds = seconds;
}

//=================================================================================
void SampleTailExtender::setDecayRate (double rate)
{
    decayRate = rate;
}

//=================================================================================
void SampleTailExtender::setBeatingStrength (double strength)
{
    beatingStrength = strength;
}

//=================================================================================
void SampleTailExtender::setCutoffThresholdInDecibels (double threshold)
{
    cutoffThresholdInDecibels = threshold;
}

//=================================================================================
std::vector<int> SampleTailExtender::getIndicesOfHarmonics (const std::vector<double>& magnitudeSpectrum, double sampleRate)
{
    std::vector<int> harmonics;
    
    double maxFrequency = 10000.;
    int maxBin = (int) (maxFrequency / (float (sampleRate) / float (audioFrameSize) ));
    
    for (int i = 2; i < maxBin; i++)
    {
        
        bool condition1 = magnitudeSpectrum[i] > magnitudeSpectrum[i - 1];
        bool condition2 = magnitudeSpectrum[i] > magnitudeSpectrum[i - 2];
        bool condition3 = magnitudeSpectrum[i] > magnitudeSpectrum[i + 1];
        bool condition4 = magnitudeSpectrum[i] > magnitudeSpectrum[i + 2];
        
        if (condition1 && condition2 && condition3 && condition4)
            harmonics.push_back (i);
    }
    
    return harmonics;
}

//=================================================================================
std::vector<double> SampleTailExtender::makeDecayEnvelope (double observedPeakEnergy, double sampleRate, double decaySpeed)
{
    std::vector<double> envelope;
    
    double observedDecibels = 20. * log10 (observedPeakEnergy);
    
    // this puts decay speed at the kind of value
    // it needs to be (if we allow the input value to
    // be closer to a normalised value)
    decaySpeed = decaySpeed * 3.84;
    
    double envelopeValueInDecibels = observedDecibels;
    double envelopeChange = decaySpeed / sampleRate;
    
    while (envelopeValueInDecibels > -100.)
    {
        double envelopeValue = pow (10., envelopeValueInDecibels / 20.);
        envelopeValue /= (observedPeakEnergy + FLT_MIN);
        envelope.push_back (envelopeValue);
        
        envelopeValueInDecibels -= envelopeChange;
    }
    
    return envelope;
}

//=================================================================================
double SampleTailExtender::calculateSplicePointPeakEnergy (std::vector<double>& audioSignal, int splicePoint)
{
    double maxVal = 0.0;
    for (int i = splicePoint - 2048; i < splicePoint + 2048; i++)
    {
        if (fabs (audioSignal[i]) > maxVal)
            maxVal = fabs (audioSignal[i]);
    }
    
    return maxVal;
}

//=================================================================================
double SampleTailExtender::getMagnitudeForSinusoid (const std::vector<double>& magnitude, int index)
{
    return 2. * (magnitude[index] / (float) (audioFrameSize / 2));
}

//=================================================================================
double SampleTailExtender::getPhaseForSinusoid (const std::vector<double>& phase, int index)
{
    double correctedPhase;
    
    if (index % 2 == 0)
        correctedPhase = phase[index];
    else
        correctedPhase = phase[index] - M_PI;
    
    while (correctedPhase < (-1. * M_PI))
        correctedPhase += (2 * M_PI);
    
    while (correctedPhase >= M_PI)
        correctedPhase -= (2 * M_PI);
    
    return correctedPhase;
}

//=================================================================================
double SampleTailExtender::getTunedFrequency (const std::vector<double>& magnitude, int index, double sampleRate)
{
    // parabolic interpolation
    double A1 = 20. * log10 (magnitude[index - 1]);
    double A2 = 20. * log10 (magnitude[index]);
    double A3 = 20. * log10 (magnitude[index + 1]);
    
    double v = (A1 - A3) / (2. * (A1 - 2. * A2 + A3));
    
    double tunedBinIndex = (double)index + v;
    double tunedFrequency = tunedBinIndex / (double)audioFrameSize * sampleRate;
    return tunedFrequency;
}

//=================================================================================
std::vector<double> SampleTailExtender::synthesiseSinusoid (int numSamples, double magnitude, double phase, double frequency, double sampleRate)
{
    std::vector<double> sinusoid(numSamples);
    
    for(int i = 0; i < numSamples; i += 1)
    {
        double p = phase + (i / sampleRate) * frequency * 2. * M_PI;
        
        // r = small random amount adjusting phase that helps to reduce the output sounding synthetic
        double r = (double)(rand() % 1000) / 1000.;
        r -= 0.5;
        r /= 1000.;
        
        p += r;
        
        double sample = magnitude * cos (p);
        sinusoid[i] = sample;
    }
    
    return sinusoid;
}

//=================================================================================
void SampleTailExtender::setSignalToZeroBelowThreshold (std::vector<double>& signal, double thresholdInDecibels)
{
	if (signal.size() == 0)
		return;

    int thresholdPoint = -1;
    int blockSize = 4096;
    
    for (size_t i = 0; i < signal.size() - blockSize; i += blockSize)
    {
        double peakEnergy = 0.;
        for (int j = 0; j < blockSize; j++)
        {
            if (fabs (signal[i + j]) > peakEnergy)
                peakEnergy = fabs (signal[i + j]);
        }
        
        double dB = 20. * log10 (peakEnergy);
        
        if (dB < thresholdInDecibels)
        {
            thresholdPoint = i;
            break;
        }
    }
    
    if (thresholdPoint != -1)
    {
        for (size_t i = thresholdPoint; i < signal.size(); i++)
        {
            int j = i - thresholdPoint;
            
            if (j < 10000) // fade out for the first 100000 samples
            {
                double g = (10000. - float (j)) / 10000.;
                signal[i] *= g;
            }
            else
            {
                signal[i] = 0.0;
            }
        }
    }
}

//=================================================================================
void SampleTailExtender::applyHarmonicBeatingToSinusoid (std::vector<double>& sinusoid, double sampleRate, double beatingAmount, double beatingFrequency)
{
    double numSeconds = static_cast<double> (sinusoid.size()) / sampleRate;
    double minValueForBeatingEnvelope = 1. - beatingAmount;
    double rangeOfBeatingEnvelope = 1. - minValueForBeatingEnvelope;
    
    for (size_t i = 0; i < sinusoid.size(); i++)
    {
        double x = static_cast<double> (i) / static_cast<double> (sinusoid.size());
        double y = cos (2. * M_PI * x * beatingFrequency * numSeconds);
        y = (y + 1.) / 2.;
        y = y * rangeOfBeatingEnvelope;
        y += minValueForBeatingEnvelope;
        
        sinusoid[i] *= y;
    }
}
