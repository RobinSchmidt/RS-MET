//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#ifndef HarmonicAnalyser_h
#define HarmonicAnalyser_h

#include <vector>

//=================================================================================
/** A class for analysing the amplitude trajectory of individual harmonics. The amount of
 * harmonic beating and the frequency of that beating are estimated. It should be noted
 * that these are approximate estimates taken from a very small amount of data.
 */
class HarmonicAnalyser
{
public:
    
    //=================================================================================
    /** Constructor */
    HarmonicAnalyser();
    
    //=================================================================================
    /** Analyses the harmonics in the input audio signal up to the splice point. The resulting amounts of beating and frequency
     * of beating are storeed in beatingAmount and beatingFrequency
     */
    void analyseHarmonics (std::vector<double>& audioSignal, double sampleRate, int splicePoint, int audioFrameSize, int notePitchClass);
    
    /** Calculates the indices of major harmonics */
    std::vector<int> getIndicesOfMajorHarmonics (std::vector<double>& audioSignal, int splicePoint, int audioFrameSize);
    
    //=================================================================================
    /** The amount of beating per bin of the magnitude spectrum. This will be mostly zero as most bins are not major harmonics */
    std::vector<double> beatingAmount;
    
    /** The frequency of beating per bin of the magnitude spectrum. This will be mostly zero as most bins are not major harmonics */
    std::vector<double> beatingFrequency;
    
private:
    
    void resetHarmonicAmplitudeProfiles();
    void analyseBeating (int audioFrameSize, int hopSize, double sampleRate, int notePitchClass);
    
    std::vector<int> harmonics;
    std::vector<std::vector<double>> harmonicAmplitudeProfiles;
};

#endif /* HarmonicAnalyser_h */
