# Sample Tail Extension (C++)

This library extends samples with a synthesised tail, given the first few seconds of the sample. 

### Dependencies

The algorithm uses Kiss FFT for fourier transform calculations. This is included in the `libs` directory.

### Usage


	#include "SampleTailExtender.h"
	
	// create the sample tail extender object
	SampleTailExtender se;
	
	// create an audio signal as a vector and fill with samples
	std::vector<double> inputSignal; 
	
	// store the sample rate as a double
	double sampleRate = 44100.;
	
	// store the pitch class as an integer, where C = 0, C# = 1, etc
	int pitchClass = 0; 
	
	// extend the sample
	std::vector<double> outputSignal = se.extendSample (inputSignal, sampleRate, pitchClass);
	
### Parameters

You can set parameters for the algorithm using the following functions:

- setSynthesisStartPointInSeconds()
- setDecayRate()
- setBeatingStrength()
- setCutoffThresholdInDecibels()

Please see `SampleTailExtender.h` which contains documentation of what each of these are used for and which values to use.
	
### Estimating Pitch Class From The Filename

If you have a file name for the sample in the following format...

`C#1.wav`

or 

`G2.wav`

then you can use two functions I've included in a header called `SampleTailExtensionHelperFunctions.h`, e.g:

	#include "SampleTailExtensionHelperFunctions.h"

	std::string filePath = "/path/to/your/audio/C#2.wav";
	std::string fileName = getFileNameFromPath (filePath);
	int pitchClass = getPitchClassFromFileName (fileName);