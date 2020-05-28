
//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#ifndef SampleTailExtender_h
#define SampleTailExtender_h

#include <iostream>
#include <vector>
#include "HarmonicAnalyser.h"

//=================================================================================
/** Main class for extending sample tails. */
class SampleTailExtender
{
public:
    
    //=================================================================================
    /** Constructor */
    SampleTailExtender();
    
    //=================================================================================
    /** Extends a given sample, returning a new vector containing the extended sample.
     * @param inputSignal the input audio signal
     * @param sampleRate the audio sample rate
     * @param notePitchClass the note pitch class where C = 0, C# = 1, ... B = 11. Please
     * see SampleTailExtensionHelperFunctions.h for some helpful functions to derive this from 
     * the file name
     * @returns a vector containing the extended sample 
     */
    std::vector<double> extendSample (std::vector<double> inputSignal, double sampleRate, int notePitchClass);
    
    //=================================================================================
    /** Sets the time in seconds at which the synthesised sample will begin. If the sample is shorter than 
     * the start time, the last possible point in the sample will be used */
    void setSynthesisStartPointInSeconds (double seconds);
    
    /** Sets the rate of decay of the sample which controls how fast it will die away. The normal value is 1.0, so 
     * adjust up and down from there.
     */
    void setDecayRate (double rate);
    
    /** This controls the strength of any beating applied (which is automatically estimated, but this parameter lets
     * you control whether it is more or less severely applied. The default value is 1.0
     */
    void setBeatingStrength (double strength);
    
    /** This is the threshold in decibels at which the tail should be set to zero to avoid bit reduction noise. 
     * The default value is -100 decibels
     */
    void setCutoffThresholdInDecibels (double threshold);

    std::vector<int> getIndicesOfHarmonics(const std::vector<double>& magnitudeSpectrum, double sampleRate);

private:
    
    //=================================================================================
    std::vector<double> makeDecayEnvelope (double observedPeakEnergy, double sampleRate, double decaySpeed);
    double calculateSplicePointPeakEnergy (std::vector<double>& audioSignal, int splicePoint);
    double getMagnitudeForSinusoid (const std::vector<double>& magnitude, int index);
    double getPhaseForSinusoid (const std::vector<double>& phase, int index);
    double getTunedFrequency (const std::vector<double>& magnitude, int index, double sampleRate);
    std::vector<double> synthesiseSinusoid (int numSamples, double magnitude, double phase, double frequency, double sampleRate);
    void setSignalToZeroBelowThreshold (std::vector<double>& signal, double thresholdInDecibels);
    void applyHarmonicBeatingToSinusoid (std::vector<double>& sinusoid, double sampleRate, double beatingAmount, double beatingFrequency);
    
    //=================================================================================
    const int audioFrameSize {8192};
    const int crossfadeLength {8192};
    double synthesisStartPointInSeconds {3.1};
    double cutoffThresholdInDecibels {-100.};
    double decayRate {1.0};
    double beatingStrength {1.0};
    std::vector<double> audioAnalysisFrame;
    HarmonicAnalyser harmonicAnalyser;
};

#endif /* SampleTailExtender_h */
