#ifndef RAPT_SAMPLEMANIPULATION_H
#define RAPT_SAMPLEMANIPULATION_H

/** Applies a smooth (sinusoidal) fade-out between samples "start" and "end". */
template<class T>
void rsFadeOut(T* buffer, int start, int end);

#endif
