//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#include "HarmonicAnalyser.h"
#include "FFT.h"
#include <iostream>

//=================================================================================
HarmonicAnalyser::HarmonicAnalyser()
{
    
}

//=================================================================================
void HarmonicAnalyser::analyseHarmonics (std::vector<double>& audioSignal, double sampleRate, int splicePoint, int audioFrameSize, int notePitchClass)
{
    FFT fft (audioFrameSize);
    int hopSize = audioFrameSize / 4;
    int i = 0;
    
    // get the indices of major harmonics
    harmonics = getIndicesOfMajorHarmonics (audioSignal, splicePoint, audioFrameSize);
    
    resetHarmonicAmplitudeProfiles();
    harmonicAmplitudeProfiles.resize (harmonics.size());
    
    // step through the signal in frames, storing the energy of each major harmonic
    // to build up an energy profile
    while (i <= splicePoint)
    {
        std::vector<double> audioFrame (audioSignal.begin() + i, audioSignal.begin() + i + audioFrameSize);
        fft.performFFT (audioFrame);
        
        const std::vector<double>& magnitude = fft.getMagnitude();
        
        // store the harmonic energy for the given frame
        for (size_t j = 0; j < harmonics.size(); j++)
            harmonicAmplitudeProfiles[j].push_back (magnitude[harmonics[j]]);
        
        i += hopSize;
    }
    
    // calculate amount and frequency of harmonic beating
    analyseBeating (audioFrameSize, hopSize, sampleRate, notePitchClass);
}

//=================================================================================
std::vector<int> HarmonicAnalyser::getIndicesOfMajorHarmonics (std::vector<double>& audioSignal, int splicePoint, int audioFrameSize)
{
    std::vector<double> audioFrame (audioSignal.begin() + splicePoint - (audioFrameSize / 2), audioSignal.begin() + splicePoint + (audioFrameSize / 2));
    
    FFT fft (audioFrameSize);
    fft.performFFT (audioFrame);

    std::vector<double> magnitudeSpectrum = fft.getMagnitude();
    
    for (size_t i = 0; i < magnitudeSpectrum.size(); i++)
        magnitudeSpectrum[i] = magnitudeSpectrum[i] / ((double)audioFrameSize / 2.);
    
    double decibelThreshold = -60.;
    double threshold = pow (10., decibelThreshold / 20.);
    
    std::vector<int> harmonics;
    
    for (size_t i = 2; i < (magnitudeSpectrum.size() / 2) - 2; i++)
    {
        bool condition1 = magnitudeSpectrum[i] > magnitudeSpectrum[i - 1];
        bool condition2 = magnitudeSpectrum[i] > magnitudeSpectrum[i - 2];
        bool condition3 = magnitudeSpectrum[i] > magnitudeSpectrum[i + 1];
        bool condition4 = magnitudeSpectrum[i] > magnitudeSpectrum[i + 2];
        bool condition5 = magnitudeSpectrum[i] > threshold;
        
        if (condition1 && condition2 && condition3 && condition4 && condition5)
            harmonics.push_back (i);
    }
    
    return harmonics;
}

//=================================================================================
void HarmonicAnalyser::resetHarmonicAmplitudeProfiles()
{
    for (size_t i = 0; i < harmonicAmplitudeProfiles.size(); i++)
        harmonicAmplitudeProfiles[i].clear();
    
    harmonicAmplitudeProfiles.clear();
}

//=================================================================================
void HarmonicAnalyser::analyseBeating (int audioFrameSize, int hopSize, double sampleRate, int notePitchClass)
{
    beatingAmount.resize (audioFrameSize);
    std::fill (beatingAmount.begin(), beatingAmount.end(), 0.);
    
    beatingFrequency.resize (audioFrameSize);
    std::fill (beatingFrequency.begin(), beatingFrequency.end(), 0.);
    
    for (size_t h = 0; h < harmonics.size(); h++)
    {
        double frequencyOfBeating;
        double amountOfBeating;
        
        // calculate the pitch class of the harmonic
        double f = harmonics[h] * sampleRate / static_cast<double>(audioFrameSize);
        int midiNote = int (round (log (f / 440.) / log (2.) * 12. + 69));
        int harmonicPitchClass = midiNote % 12;
        
        int positiveChangeCount = 0;
        
        // calculate the difference and count how many instances have positive change
        for (size_t i = 1; i < harmonicAmplitudeProfiles[h].size(); i++)
        {
            double difference = harmonicAmplitudeProfiles[h][i] - harmonicAmplitudeProfiles[h][i - 1];
            
            if (difference >= 0)
                positiveChangeCount++;
        }
        
        // the amount of beating is estimated as the number of positive changes as a ratio of half the harmonic profile length
        // (i.e. we expect the harmonic energy to decrease, if it increases in lots of instances, it is probably beating)
        amountOfBeating = static_cast<double>(positiveChangeCount) / static_cast<double>(harmonicAmplitudeProfiles[h].size() / 2);
    
        
        std::vector<double> flattenedHarmonicAmplitudeCurve;
        
        // apply some flattening to the curve before measuring peaks (we don't want the downwards trend to make some peaks unrecognised)
        for (size_t i = 0; i < harmonicAmplitudeProfiles[h].size(); i++)
            flattenedHarmonicAmplitudeCurve.push_back (harmonicAmplitudeProfiles[h][i] * (1. + 0.01 * static_cast<double> (i)));
        
        std::vector<int> peakArray;
        peakArray.push_back (0); // consider first sample a peak
        
        for (size_t i = 2; i < flattenedHarmonicAmplitudeCurve.size() - 2; i++)
        {
            bool condition1 = flattenedHarmonicAmplitudeCurve[i] > flattenedHarmonicAmplitudeCurve[i - 1];
            bool condition2 = flattenedHarmonicAmplitudeCurve[i] > flattenedHarmonicAmplitudeCurve[i - 2];
            bool condition3 = flattenedHarmonicAmplitudeCurve[i] > flattenedHarmonicAmplitudeCurve[i + 1];
            bool condition4 = flattenedHarmonicAmplitudeCurve[i] > flattenedHarmonicAmplitudeCurve[i + 2];
            
            if (condition1 && condition2 && condition3 && condition4)
                peakArray.push_back (i);
        }
        
        //  if we have some peaks, attempt to estimate the "frequency"
        // (the beating doesn't seem to have a very regular frequency from what I have
        // observed - it varies, sometimes drastically, but this should give a rough estimate)
        if (peakArray.size() > 1)
        {
            double sumVal = 0.;
            
            for (size_t k = 1; k < peakArray.size(); k++)
                sumVal += (peakArray[k] - peakArray[k - 1]);
            
            double meanDifference = sumVal / static_cast<double> (peakArray.size());
            frequencyOfBeating = 1. / (meanDifference * (static_cast<double> (hopSize) / sampleRate));
        }
        else // if we don't have any peaks, just randomly estimate the frequency, there probably isn't beating anyway here...
        {
            double randomValue = static_cast<double> (rand() % 1000) / 1000.;
            frequencyOfBeating = 1.5 * randomValue + 0.5; // between 0.5 and 2 Hz
        }
        
        // ensure frequency is in a sensible range
        frequencyOfBeating = fmax (0.5, frequencyOfBeating);
        frequencyOfBeating = fmin (2.0, frequencyOfBeating);
        
        // stop harmonic from beating if it might be the fundamental
        if (h < 2 && harmonicPitchClass == notePitchClass)
            amountOfBeating = 0.0;
        
        beatingAmount[harmonics[h]] = amountOfBeating;
        beatingFrequency[harmonics[h]] = frequencyOfBeating;
    }
    
}
