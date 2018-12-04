//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#ifndef FFT_h
#define FFT_h

#include "libs/kiss_fft130/kiss_fft.h"
#include <vector>

//=================================================================================
/** A wrapper around Kiss FFT for easy FFT computation */
class FFT
{
public:
    
    //=================================================================================
    /** Constructor - pass in the audio frame size you will use for FFTs */
    FFT (int frameSize);
    ~FFT();
    
    //=================================================================================
    /** Performs the FFT on an audio frame */
    void performFFT (std::vector<double> audioFrame);
    
    /** Returns the magnitude spectrum */
    const std::vector<double>& getMagnitude();
    
    /** Returns the phase spectrum */
    const std::vector<double>& getPhase();
    
private:
    
    void calculateHanningWindow();
    
    int audioFrameSize;
    std::vector<double> hanningWindow;
    std::vector<double> real;
    std::vector<double> imag;
    std::vector<double> magnitude;
    std::vector<double> phase;
    
    kiss_fft_cfg cfg;
    kiss_fft_cpx* fftIn;
    kiss_fft_cpx* fftOut;
};

#endif /* FFT_h */
