#ifndef rojue_AudioSampleBufferFunctions_h
#define rojue_AudioSampleBufferFunctions_h

#include "rojue_ImmediatePlaybackAudioSource.h"

namespace rojue
{

  /** This file contains various functions (among them signal processing functions) to be applied 
  to an AudioSampleBuffer. */

  /** Copies the contents of one AudioSampleBuffer (source) into another AudioSampleBuffer 
  (destination), changing the also size (numChannels, numSamples) of the destination buffer if 
  necesarry. */
  void copyAudioSampleBufferContents(AudioSampleBuffer source, AudioSampleBuffer destination);

  /** Resamples the input audio-buffer. The returned buffer will have to be deleted by the 
  caller (function not yet tested). */
  AudioSampleBuffer* resampleBuffer(AudioSampleBuffer* inBuffer, double inSampleRate, 
    double outSampleRate);

  /** Saves the audio-data contained in bufferToSave into an audio-file given by fileToSaveTo. The 
  file format will be inferred from the extension of the file. The last parameter can be used to
  pop up a warning box (with OK and cancel options), if the file already exists - if false, it 
  will be overwritten without warning. */
  void saveAudioSampleBufferToFile(AudioSampleBuffer* bufferToSave, File fileToSaveTo, 
    double sampleRateToUse, int bitsPerSample, bool warnIfFileExists);

}

#endif 