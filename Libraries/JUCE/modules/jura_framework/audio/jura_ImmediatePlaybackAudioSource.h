#ifndef jura_ImmediatePlaybackAudioSource_h
#define jura_ImmediatePlaybackAudioSource_h

/** This class can be used to play an AudioSampleBuffer or an AudioFileBuffer immediately. */

class ImmediatePlaybackAudioSource : public AudioSource
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ImmediatePlaybackAudioSource();

  /** Destructor. */
  virtual ~ImmediatePlaybackAudioSource();

  //---------------------------------------------------------------------------------------------
  // overrides for subclasses AudioSource:

  /** Tells the source to prepare for playing. */
  virtual void prepareToPlay(int samplesPerBlockExpected, double sampleRate);

  /** Allows the source to release anything it no longer needs after playback has stopped. */
  virtual void releaseResources();

  /** Called repeatedly to fetch subsequent blocks of audio data. */
  virtual void getNextAudioBlock(const AudioSourceChannelInfo &bufferToFill);

  //---------------------------------------------------------------------------------------------
  // others:

  /** Passes an AudioSampleBuffer for immediate playback. */
  virtual void startPlayback(AudioSampleBuffer* newBufferToPlay);

  /** Pauses the playback for later continuation (maybe also used to releasr the pause state by
  passing false. */
  virtual void pausePlayback(bool shouldBePaused = true);

  /** Resumes playbabck from the position where it currently is. */
  virtual void resumePlayback();

protected:

  AudioSampleBuffer audioSampleBuffer;
  int               playPosition;
  CriticalSection   lock;
  bool              isPaused;

  juce_UseDebuggingNewOperator;
};

#endif  