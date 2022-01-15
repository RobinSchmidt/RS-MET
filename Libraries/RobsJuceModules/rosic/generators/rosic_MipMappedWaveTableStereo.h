#ifndef rosic_MipMappedWaveTableStereo_h
#define rosic_MipMappedWaveTableStereo_h


namespace rosic
{

/** This is a class for storing a stereo waveform, mip-mapping it for anti-aliasing and retrieve
anti-aliased values from it. We do not decimate the versions of the table with less 
high-frequency content. All levels are at the same internal resolution, just with less and less
spectral content.

References:
  https://en.wikipedia.org/wiki/Mipmap

todo: 
-factor out:
 -WaveTableManipulator (does all these time and spectral manipulations)
 -WaveTableMipMapper (renders the set of mip-mapped versions)
-make a class MultiWaveTable that can crossfade between an arbitrary number of "prototype"
 waveforms (each of them with the same manipulations applied and each of them mip-mapped)
-maybe make some time-domain manipulations realtime modulatable, i.e. don't render them into 
 the wavetable(s)  */

class MipMappedWaveTableStereo
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Constructor. */
  MipMappedWaveTableStereo();

  /** Destructor. */
  ~MipMappedWaveTableStereo();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Copies the values in newWaveform into the internal buffers and renders various bandlimited
  versions via FFT/iFFT. The newWaveform array is expected to be an array of two pointers to
  doubles each of which representing one channel. That is: newWaveform[0][0] must be the first
  sample of the left channel data and newWaveform[1][0] must be the first sample of the right
  channel. Similarly, newWaveform[0][lengthInSamples-1] and newWaveform[1][lengthInSamples-1]
  must be the last samples for the two channels. The "lengthInSamples" can be arbitrary. The
  data will be interpolated to the internal table-sizes anyway by means of forward and inverse
  FFTs, where for the forward transform a Bluestein algorithm will be used to accomodate for
  the arbitrary length. The internal length is a power of two, so a radix-2 algorithm can be
  used for the iFFTs. */
  void setWaveform(double** newWaveform, int lengthInSamples);

  /** Same as above but accepts single precision buffers. */
  void setWaveform(float** newWaveform, int lengthInSamples);

  /** Sets the name (i.e. the relative path) of the current sample) - this function is only
  for conveniently handling preset management in a plugIn-context */
  void setSampleName(char* newSampleName);

  /** Time-reverses the waveform. */
  void setTimeReverse(bool shouldBeReversed);

  /** Inverts the polarity of the waveform. */
  void setPolarityInversion(bool shouldBeInversed);

  /** Sets a time warping coeffcient 'a' which maps the phase pointer according to the formula:
  ..... - this is a mapping of the domain 0...1 into the range 0...1
  which is the identity mapping for a=0, introduces concavity for a<0 and convexity for a>0.
  The coefficient should be in the range ]-1...+1[ (not including the extremes). */
  void setFullWavePhaseWarp(double newWarpCoefficient);

  /** Sets the coefficient 'b' in the formula: .... for the time warping of the phase. */
  void setHalfWavePhaseWarp(double newWarpCoefficient);

  /** Sets the harmonic number to which the comb filter is tuned. */
  void setCombHarmonic(double newCombHarmonic);

  /** Sets an offset between two instances of the waveform that are added together to form some
  kind of comb-filtering (in degrees). */
  void setCombOffset(double newCombOffset);

  /** Sets the amount of comb-filtering in percent -100...100. */
  void setCombAmount(double newCombAmount);

  /** Switches into a mode in which in the channels of the input table are converted into sum
  and difference before creating the mipmap. In this mode, the first channel output will
  represent the sum and second will represent the difference.  */
  void setChannelSumAndDifferenceMode(bool shouldUseSumAndDifference);

  /** Adjusts the ratio between even and odd harmonics. */
  void setEvenOddRatio(double newRatio);

  /** Applies a sharp (brickwall) lowpass-filter to the harmonics of the waveform. */
  void setHighestHarmonicToKeep(int newHighestHarmonicToKeep);

  /** Applies a sharp (brickwall) highpass-filter to the harmonics of the waveform. */
  void setLowestHarmonicToKeep(int newLowestHarmonicToKeep);

  /** Sets the spectral contrast - this contrast value acts as an exponent on the magnitude of
  each harmonic, thus, a value of 1 is neutral, values > 1 emphasize differences in magnitude,
  values < 1 even out differences in magnitude, a value of 0 whitens the magnitude spectrum
  and negative values invert the spectrum. Contrast is applied before all the other spectral
  modifiers (such as slope, even/odd ratio, etc.). */
  void setSpectralContrast(double newContrast);

  /** Sets a slope for the spectrum in dB/oct for the tableset to be generated - negative values
  will cause the tables to get a lowpassish character, positive value lead to a highpassish
  character. */
  void setSpectralSlope(double newSlope);

  /** Sets a phase offset between even and odd harmonics */
  void setEvenOddPhaseShift(double newPhaseShift);

  /** Sets a scaling factor for the harmonis phases. */
  void setPhaseScale(double newPhaseScale);

  /** Sets a constant shifting offset for the phases of the harmonics. */
  void setPhaseShift(double newPhaseShift);

  /** Sets the phase-shift of the the harmonics between left and right channel in degrees. The
  left channel's harmoincs will be shifted half this value into the negative direction and the
  right channel's harmonics by the same amount into the positive direction. */
  void setStereoPhaseShift(double newPhaseShift);

  /** Sets an additional phase-shift between even and odd harmonics, half of which is applied
  with negative sign to the left channel and with positive sign to the right channel. */
  void setEvenOddStereoPhaseShift(double newPhaseShift);

  /** Fills the "tableSet"-variable with all zeros. */
  void fillWithAllZeros();

  /** Switches the automatic re-rendering of the mip-map on and off - if on, the mip map will be
  automatically re-rendered whenever a parameter changes which will invalidate the current
  mip-map. You may want to switch this off, when changeing several parameters at a time, such as
  in recalling a preset, to avoid the table to be re-rendered multiple times. */
  void setAutomaticMipMapReRendering(bool shouldAutomaticallyReRender);

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the length of the table. */
  int getTableLength() const { return tableLength; }

  /** Returns the number of stored levels of detail. */
  int getNumLevels() const { return numTables; }

  /** Copies the content of the given channel and mip-mapping level into the given buffer where 
  the buffer must have a length of at least getTableLength(), channel must be 0 or 1 and level 
  must be < getNumLevels(). Useful for visualizing the different levels of the mip-map. */
  void copyDataTo(double* buffer, int channel, int level) /*const*/;

  /** Returns the name (i.e. the relative path) of the sample as a zero-terminated string. */
  char* getSampleName() const { return sampleName; }

  /** Informs, if the waveform is time-reversed. */
  bool isTimeReversed() const { return timeReversed; }

  /** Informs, if the waveform is polarity inverted. */
  bool isPolarityInverted() const { return polarityInverted; }

  /** Return the 'a' coefficient in the phase warping formula. @see: setFullWavePhaseWarp() */
  double getFullWavePhaseWarp() const { return fullWavePhaseWarp; }

  /** Return the 'b' coefficient in the phase warping formula. @see: setHalfWavePhaseWarp() */
  double getHalfWavePhaseWarp() const { return halfWavePhaseWarp; }

  /** Returns the harmonic to which the comb is tuned. */
  double getCombHarmonic() const { return combHarmonic; }

  /** Returns the offset of the added duplicate of the waveform (in degrees). */
  double getCombOffset() const { return combOffset; }

  /** Returns the amplitude of the added duplicate of the waveform (in percent). */
  double getCombAmount() const { return combAmount; }

  /** Returns a constant pointer to a pointer to the prototype wavetable. Concpetuall, it's a 
  double** where you can't modify the contents of the arrays. Use getPrototypeNumSamples() to 
  retrieve the lengths of the arrays pointed to. */
  const double* const * getPrototypeWaveform() const { return prototypeWave; }
  // ToDo: verify that this is the correct kind of constness. Maybe we should have a 3rd const in 
  // between double and *? And/or maybe the 1st const is not needed? What about thread-safety?

  /** Returns the number of samples in the prototype waveform. */
  int getPrototypeNumSamples() /*const*/;

  /** Returns the ratio between even and odd harmonics. */
  double getEvenOddRatio() const { return evenOddRatio; }

  /** Returns the highest harmonic that will be kept in the waveform. */
  int getHighestHarmonicToKeep() const { return highestHarmonicToKeep; }

  /** Returns the lowest harmonic that will be kept in the waveform. */
  int getLowestHarmonicToKeep() const { return lowestHarmonicToKeep; }

  /** Returns the spectral contrast @see: setSpectralContrast. */
  double getSpectralContrast() const { return spectralContrast; }

  /** Returns the spectral slope @see: setSpectralSlope. */
  double getSpectralSlope() const { return spectralSlope; }

  /** Returns the phase shift between even and odd harmonics @see: setEvenOddPhaseShift.  */
  double getEvenOddPhaseShift() const { return evenOddPhaseShift; }

  /** Returns the or the odd harmonics. */
  //double getOddPhaseShift();

  /** Returns the scaling factor fo the phases. */
  double getPhaseScale() const { return phaseScale; }

  /** Returns the phase shifting offset. */
  double getPhaseShift() const { return phaseShift; }

  /** Returns the phase shift between left and right channel @see: setStereoPhaseShift. */
  double getStereoPhaseShift() const { return stereoPhaseShift; }

  /** Returns the additional phase-shift between even and odd harmonics, half of which is
  applied with negative sign to the left channel and with positive sign to the right channel.
  @see: setEvenOddStereoPhaseShift. */
  double getEvenOddStereoPhaseShift() const { return evenOddStereoPhaseShift; }

  /** Informs, whether or not the mip-map is automatically re-rendered whenever a parameter
  changes which invalidates the current mip-map.  */
  bool isMipMapAutoReRenderingActive() const { return autoReRenderMipMap; }

  //-----------------------------------------------------------------------------------------------
  // \name Processing:

  /** Returns the value at 'position' of table 'tableIndex' with interpolation. */
  INLINE double getValue(int channel, int tableIndex, double position);

  /** Returns the value at position = 'positionInt' + 'positionFrac' of table 'tableIndex'
  with interpolation. */
  INLINE double getValue(int channel, int tableIndex, int positionInt, double positionFrac);

  /** Returns the value at position = 'positionInt' + 'positionFrac' of table 'tableIndex'
  with interpolation without safety checks - it is therefore faster but less secure. */
  INLINE double getValueFast(int channel, int tableIndex, int positionInt, double positionFrac);

  //-----------------------------------------------------------------------------------------------
  // \name Misc

  /** Generates a mip-mapped multisample from the prototype tables, where each of the
  successive tables contains one half of the spectrum of the previous one. */
  void renderMipMap();

protected:

  /** Returns a value of the prototype wave at some (possibly non-integer) phase index measured
  in samples (range: 0...prototypeWaveNumSamples). It uses linear interpolation */
  double getPrototypeValueAt(int channel, double phaseIndex);

  /** Warps a phase index measured in samples according to the warping map (input (and output)
  range: 0...prototypeWaveNumSamples). */
  double warpPhaseIndex(double unwarpedIndex);

  /** Number of channels - not to '2' as magic number in the actual code may facilitate later
  genralizations to more than two channels. */
  static const intA numChannels = 2;

  /** The Oscillator class uses a one table-per octave multisampling to avoid aliasing.
  With a table-size of 8192 and a sample-sample rate of 44100, the 12th table will have a
  fundamental frequency (the frequency where the increment is 1) of 11025 which is good for
  the highest frequency. */
  static const intA numTables = 10;    // rename to numLevels
  //static const intA numTables = 12;

  /** Length of the lookup-table. The actual length of the allocated memory is 4 samples longer,
  to store additional samples for the interpolator (which are the same values as at the beginning
  of the buffer) */
  static const intA tableLength = 2048;
  //static const intA tableLength = 8192;

  /** The multisample for anti-aliased waveform generation. The 4 additional values are equal
  to the first 4 values in the table for easier interpolation. The first index is for the
  table-number - index 0 accesses the first version which has full bandwidth, index 1 accesses
  the second version which is bandlimited to Nyquist/2, 2->Nyquist/4, 3->Nyquist/8, etc. */
  float tableSet[numTables][tableLength+4][numChannels];

  char* sampleName;

  double* prototypeWave[2];
  int     prototypeWaveNumSamples;

  bool   channelSumAndDifferenceMode, timeReversed, polarityInverted;
  double fullWavePhaseWarp, halfWavePhaseWarp;
  double combHarmonic;
  double combOffset, combAmount;
  double spectralContrast, evenOddRatio, spectralSlope;
  int    lowestHarmonicToKeep, highestHarmonicToKeep;
  double phaseScale, phaseShift, evenOddPhaseShift, stereoPhaseShift, evenOddStereoPhaseShift;

  bool autoReRenderMipMap;

  // embedded objects:
  FourierTransformerBluestein forwardTransformer;
  FourierTransformerRadix2    inverseTransformer;
  //Interpolator                interpolator;

  MutexLock mutex;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE double MipMappedWaveTableStereo::getValue(int channel, int tableIndex, double position)
{
  double positionFloor = floor(position);
  double positionFrac  = position - positionFloor;
  int    positionInt   = (int)positionFloor;

  return getValue(channel, tableIndex, positionInt, positionFrac);
}

INLINE double MipMappedWaveTableStereo::getValue(int channel, int tableIndex, int positionInt,
  double positionFrac)
{
  // ensure, that the channel and table index are in the valid range:
  if(channel < 0)
    channel = 0;
  else if(channel >= numChannels)
    channel = numChannels;

  if(tableIndex < 0)
    tableIndex = 0;
  else if(tableIndex >= numTables)
    tableIndex = numTables-1;

  // lookup value in the table with interpolation and return it:
  return getValueFast(channel, tableIndex, positionInt, positionFrac);
}

INLINE double MipMappedWaveTableStereo::getValueFast(int c, int t, int pi, double pf)
{
  //return (1.0-pf)*tableSet[t][pi][c] + pf*tableSet[t][pi+1][c];
  return tableSet[t][pi][c] + pf*(tableSet[t][pi+1][c]-tableSet[t][pi][c]);
}

} // end namespace rosic

#endif // rosic_MipMappedWaveTableStereo_h
