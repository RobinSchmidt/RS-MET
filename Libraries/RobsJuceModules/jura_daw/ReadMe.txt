If i will at any point write a DAW application, classes for that should go here.



Idea: make a RAPT-DAW

Main concepts from the user's point of view:

Clip:
-represents a short snippet of audio or control data that can be moved around
 in the arrangement
-can have any number of subclips (recursively)
-stores a time (in beats) for where it is located in the arrangement
-can convert that time in beats into a time in samples by using a referenced
 object of a class TimeConverter
-has subclasses AudioClip, ControlClip, etc.
-AudioClips represent a region of audio from a buffer in the SamplePool
-maybe can also have an effects chain, maybe here we can even use offline
 processors
-overlapping clips on the same track will just get added

Track:
-represents one channel of audio
-can have any number of subtracks (recursively)
-clips can be positioned on it
-has an effects chain associated with it (represented by a class 
 RealtimeProcessingChain)

Arrangement:
-represents a full arrangement with any number of tracks
-can be wrapped up into a clip and subsequently be used in a higher-level 
 arrangement
-is a subclass of TimeConverter, such that clips can reference and use their 
 embedding arrangement to convert their stored time in beats into actual time
 in samples (that makes it possible to make global tempo changes on the 
 arrangement level where the clips can still maintain their musical time in beats)

AudioPool:
-stores all the audio-data that is used in the project (in clips or sampling instruments)
-audio clips hold a reference to an AudioBuffer that sits in the AudioPool
-user can grab samples from the pool (in a kind of browser) and place them into
 the arrangement
-the browser can also show a disk directory in another subsection where the user can 
 also grab files from - if they are not already in the pool, they will be automatically
 added (and optionally a local copy in the project folder is created)
 
AudioFileBuffer 
-subclass of juce::AudioBuffer<float>
-used inside the AudioPool for accessing the audio data
-manages direct-from-disk streaming, if necessary (it should be an option of the project
 settings - maybe a preload/buffering size in  samples which can also bet set to 
 "Everything" of "Full" - or maybe that option could be per buffer?)





Main concepts from the programmers point of view:

Software Architechture:
-use the PAC pattern (presentation/abstraction/control)

Threading Model:
-1 main thread that handles user input (GUI-stuff)
-a number of audio threads equal to the number of CPU cores, the processing of
 the individual tracks will be distributed among the threads (and thereby on 
 the CPU cores), we will need some kind of ThreadScheduler class for this
-maybe a drawing thread that is used to draw the waveforms and stuff on the GUI 
 it's not a good idea to do that in the main thread because it can make the GUI
 unresponsive when there's a lot to draw (the waveforms in the clips, etc.)
-maybe a thread for streaming/preloading audio data
-use a lock-free queue for synchronization -> the main/GUI thread produces data 
 that is consumed by the audio threads


maybe call it PlayDaw
