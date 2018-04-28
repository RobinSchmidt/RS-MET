#ifndef RS_SAMPLEMANIPULATION_H
#define RS_SAMPLEMANIPULATION_H

namespace RSLib
{

  /** Applies a smooth (sinusoidal) fade-out between samples "start" and "end". */
  void RSLib_API rsFadeOut(double *buffer, int start, int end);

}

#endif
