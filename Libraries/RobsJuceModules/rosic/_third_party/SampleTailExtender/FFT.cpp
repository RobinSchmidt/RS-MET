//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#define _USE_MATH_DEFINES
#include "FFT.h"
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

//=================================================================================
FFT::FFT (int frameSize)
{
    // store audio frame size
    audioFrameSize = frameSize;

    // setup FFT
    fftIn = new kiss_fft_cpx[audioFrameSize];
    fftOut = new kiss_fft_cpx[audioFrameSize];
    cfg = kiss_fft_alloc (audioFrameSize, 0, 0, 0);

    // resize vectors
    real.resize (audioFrameSize);
    imag.resize (audioFrameSize);
    magnitude.resize ((audioFrameSize / 2) + 1);
    phase.resize ((audioFrameSize / 2) + 1);

    // generate hanning window
    calculateHanningWindow();
}

//=================================================================================
FFT::~FFT()
{
    free (cfg);
    delete[] fftIn;
    delete[] fftOut;
}

//=================================================================================
void FFT::performFFT (std::vector<double> audioFrame)
{
    // copy samples to FFT input array
    for (int i = 0; i < audioFrameSize; i++)
    {
        fftIn[i].r = (float)(audioFrame[i] * hanningWindow[i]);
        fftIn[i].i = 0.0;
    }

    // perform the FFT
    kiss_fft (cfg, fftIn, fftOut);

    // retrieve samples from FFT output array
    for (int i = 0; i < audioFrameSize; i++)
    {
        real[i] = (double)fftOut[i].r;
        imag[i] = (double)fftOut[i].i;
    }

    // calculate magnitude and phase
    for (size_t i = 0; i < magnitude.size(); i++)
    {
        magnitude[i] = sqrt ((real[i] * real[i]) + (imag[i] * imag[i]));
        phase[i] = atan2 (imag[i], real[i]);
    }
}

//=================================================================================
const std::vector<double>& FFT::getMagnitude()
{
    return magnitude;
}

//=================================================================================
const std::vector<double>& FFT::getPhase()
{
    return phase;
}

//=================================================================================
void FFT::calculateHanningWindow()
{
    hanningWindow.resize (audioFrameSize);

    double divisor = (double) (audioFrameSize - 1);

    for (int i = 0; i < audioFrameSize; i++)
        hanningWindow[i] = 0.5 * (1 - cos (2. * M_PI * (i / divisor)));
}
