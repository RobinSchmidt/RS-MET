#ifndef RS_PARTIALEXTRACTIONEXPERIMENTS_H
#define RS_PARTIALEXTRACTIONEXPERIMENTS_H


/** Plots the impulse- or magnitude responses of the multipass bidirectional (bandpass or lowpass) 
filter for different numbers of passes. */
void biDirectionalFilter();


void beatingSines();
void envelopeDeBeating(); // try to remove beating from an extracted envelope


/** Creates a sine wave with an amplitude envelope and tries to retrieve that envelope and recreate
the sine from the envelope (with possibly different frequency and startphase). */
void sineRecreation();

void sineRecreationBandpassNoise();

void sineWithPhaseCatchUp();


/** Creates 3 partials with attack/decay envelope and tries to recover the middle one by means of 
the rsBiDirectionalFilter::applyConstPeakBandpassBwInHz function. */
void partialExtractionTriple();





void partialExtractionBell();

/** Creates 3 partials with attack/decay envelope and tries to recover the middle one via a biquad
filter. */
void partialExtractionViaBiquadTriple();

void partialExtractionSample();
// rename functions...

#endif