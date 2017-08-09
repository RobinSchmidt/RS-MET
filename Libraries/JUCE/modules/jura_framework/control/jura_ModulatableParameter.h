#ifndef jura_ModulatableParameter_h
#define jura_ModulatableParameter_h

/*

Idea to extend the same concept of MetaParameters to a modulation system for synthesizers:

ModulationSource:
-provide interface to pull out a current value
-subclasses can be fixed values (like gui parameter sliders), envelope generators, low frequency 
 oscillators, step sequencers, midi-note values, etc.
-at each sample, the first thing to do is to compute all current values of all existing modulation 
 sources

ModulationTarget:
-is typically some (lower level) dsp algorithm parameter such as a cutoff frequency of a filter
-keeps a list of assigned sources
-keeps a pointer to ModulationManager where it can sign up to receive inputs from sources
-when a sample is produced, pulls out values of all its connected sources and combines them and 
 then sets up the algorithm parameter accordingly
-perhaps combination can additively for some sources and multiplicatively for others: 1st add up 
 all additive sources, then multiply in all multiplicative sources - keep 2 lists
-can regulate the amount of the influence of each source
-maybe can apply different response curves to each source
-hmm...maybe 4 the last 2 points, we need an additional class: ModulationConnection

ModulationManager:
-allows ModulationsTargets to de/register themselves to ModulationSources

similarities to MetaParameters:
-needs similar sign-up system on the gui (right-click - assign, etc.)

differences to MetaParameters:
-each modulation target can be connected to any number of modulation sources (many-to-many instead 
 of one-to-many)
-uses pull rather than push mechanism to update the dependent parameter

feedback:
-it is desirable to be able for (parameters of) ModulationsSources themselves be modulated by 
 other ModulationSources or even by their own output
-maybe the parameters of a ModulationSource (like an LFO freq) should be subclasses of 
 ModulationTarget ...or just BE a ModulationTarget object - or maybe we'll need a class
 ModulatableParameter


maybe make the class hierarchy like this:
Parameter <- MetaControlledParameter <- ModulatableParameter <- PolyphonicParameter

*/


//=================================================================================================

/**   */

class JUCE_API ParameterModulator // : public Modulator - later, we may factor out a Modulator class that has nothing ot do with parameters
{

public:

  /** Constructor */
  ParameterModulator() {}


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterModulator)
};

//=================================================================================================

/**  */

class JUCE_API ModulatableParameter : public Parameter
{

public:

  /** Constructor */
  ModulatableParameter(const juce::String& name, double min = 0.0, double max = 1.0,
    double defaultValue = 0.5, int scaling = LINEAR, double interval = 0.0)
    : Parameter(name, min, max, defaultValue, scaling, interval) {}




protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ModulatableParameter)
};

//=================================================================================================

/**  */

class JUCE_API ParameterModulationManager
{

public:

  ParameterModulationManager() {};


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterModulationManager)
};

#endif
