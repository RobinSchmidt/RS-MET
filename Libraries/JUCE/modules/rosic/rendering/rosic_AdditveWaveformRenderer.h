#ifndef rosic_WaveformGenerator_h
#define rosic_WaveformGenerator_h

// rosic-indcludes:
#include "../transforms/rosic_FourierTransformerRadix2.h"
#include "../math/rosic_PrimeNumbers.h"

namespace rosic
{









  /**

  This is a class for generating a single-cycle-waveform by means of various algorithms.

  */

  class WaveformGenerator
  {

  public:

    /** Enumeration of the modes that can be used to create the waveforms. */
    enum modes
    {
      STANDARD_WAVEFORM,
      AUDIO_FILE,
      ALGORITHM,
      TIME_DOMAIN_FORMULA,
      FREQUENCY_DOMAIN_FORMULA,
      MULTI_SEGMENT,

      NUM_MODES
    };


    /** Enumeration of the algorithms that can be used to create the waveforms. */
    enum algorithms
    {
      AMPLITUDE_MODULATION,
      OCTAVES_OF_PRIMES,
      PHASE_MODULATION,
      WINDOW_FUNCTION,
      //.... more to come

      NUM_ALGORITHMS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WaveformGenerator();          

    /** Destructor. */
    ~WaveformGenerator();         

    //---------------------------------------------------------------------------------------------
    // parameter-settings:



    /** Chooses one of the algorithms for wavefrom generation. @see: algorithms. */
    void setAlgorithm(int newAlgorithmIndex);

    /** Sets the shape of one of the waveforms that is to be used as basic waveforms in one of the
    algorithms. Depending on the chosen algorithm, this can have different meanings - for example, 
    in PHASE_MODULATION algorithm, index==0 would refer to the carrier-waveform and index==1 would 
    refer to the modulator waveform. @see basicWaveforms */
    void setBasicWaveshape(int index, int newShape);

    /** Sets one of the algorithm's parameters. Depending on the chosen algorithm, this can have 
    different meanings - for example, in PHASE_MODULATION algorithm, index==0 would refer to the 
    (relative) frequency of the carrier, index==1 refers to the (relative) frequency of the 
    modulator and index==2 refers to the modulation index. */
    void setParameter(int index, double value);

    /** Sets the shape of the carrier wave for the phase- and amplitude-modulation algorithms. */
    void setCarrierWaveform(int newWaveformIndex) { setBasicWaveshape(0, newWaveformIndex); }

    /** Sets the shape of the modulator wave for the phase- and amplitude-modulation algorithms. */
    void setModulatorWaveform(int newWaveformIndex) { setBasicWaveshape(1, newWaveformIndex); }

    /** Sets the relative frequency of the carrier wave for the phase- and amplitude-modulation
    algorithms. */
    void setRelativeCarrierFrequency(double newFrequency) { setParameter(0, newFrequency); }

    /** Sets the relative frequency of the modulator wave for the phase- and amplitude-modulation
    algorithms. */
    void setRelativeModulatorFrequency(double newFrequency) { setParameter(1, newFrequency); }

    /** Sets the modulation index for the phase- and amplitude-modulation algorithms. */
    void setModulationIndex(double newIndex) { setParameter(2, newIndex); }

    /** Sets the lowest prime-number for the OCTAVE_OF_PRIMES algorithm. */
    void setLowestPrime(int newLowestPrime) { setParameter(0, newLowestPrime); }

    /** Sets the highest prime-number for the OCTAVE_OF_PRIMES algorithm. */
    void setHighestPrime(int newHighestPrime) { setParameter(1, newHighestPrime); }

    //---------------------------------------------------------------------------------------------
    // waveform retrieval:

    /** Returns the length of the waveform. */
    int getTableLength() const { return tableLength; }

    /** Returns a pointer to the waveform. */
    double* getWaveForm() { return wave; } 

    //=============================================================================================

  protected:

    /** Fills the array with the waveform with all zeros. */
    void clearWaveform();

    /** Fills the array with the magnitude spectrum with all zeros. */
    void clearMagnitudes();

    /** Fills the array with the phase spectrum with 'phaseValue'. */
    void initPhases(double phaseValue);

    /** Renders the waveform according to the selected algorithm, its parameters and the chosen 
    basic waveshapes. */
    void renderWaveform();


    int    algorithmIndex;                // index of the selected algorithm
    static const int tableLength = 2048;  // length of the waveform
    double wave[tableLength];             // the actual waveform
    double magnitudes[tableLength/2];     // magnitude spectrum for the spectral algorithms
    double phases[tableLength/2];         // magnitude spectrum for the spectral algorithms
    static const int numParameters = 3;   // (maximum) number of parameters for the algorithms
    double p[numParameters];              // values of the parameters
    double (*wave1) (double);             // function pointers to the basic waveform-generation 
    double (*wave2) (double);             // functions

    FourierTransformerRadix2 fourierTransformer; 
      // used for the algorithms that operate in the frequency domain

    //---------------------------------------------------------------------------------------------
    // the actual waveform rendering algorithms implemented as member-functions:

    /** Renders a waveform via amplitude-modulation. p[0] is the relative frequency of the carrier, 
    p[1] the relative frequency of the modulator and p[3] is the modulation index. */
    void renderWaveformAmplitudeModulation();

    /** Renders a waveform that contains harmonics at octaves of prime numbers. If the chosen 
    waveform is not sinusoidal, the 'harmonics' also will actually be non-sinusoidal and thus, the 
    term does not really apply here. */
    void renderWaveformOctavesOfPrimes();

    /** Renders a waveform via phase-modulation. p[0] is the relative frequency of the carrier, 
    p[1] the relative frequency of the modulator and p[3] is the modulation index. */
    void renderWaveformPhaseModulation();

    /** Renders a waveform by applying a raised cosine window to a number of wavecycles of the 
    underlying carrier waveform. p[0] is the number of cycles in the carrier wave that would fit into 
    the cycle of the waveform, p[1] is the number of actually 'visible' cycles inside the windowing 
    cosine wave. */
    void renderWaveformWindowFunction();

  };

} // end namespace rosic

#endif // rosic_WaveformGenerator_h
