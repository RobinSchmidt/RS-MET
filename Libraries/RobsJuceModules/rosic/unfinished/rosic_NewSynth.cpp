

/*
void rsNewSynth::getSampleFrameStereo(double* left, double* right)
{

}
*/


/*

General code architechture:


hmmm...maybe it would be a better idea to just somehow extend the existing modulation system to 
allow for polyphony and then implement the QuadSource and DualFilter as modules for ToolChain. 
Then, also implement a PolyModulator module that can be plugged in and that supports plugging in 
many modulations sources (similar to Elan's plugins).
-maybe all AudioModules should support dryGain, wetGain parameters where in the case of source
 modules, both are 1 by default and for effect plugins dry=0, wet=1 by default. The idea is that
 sources can be stacked and effects can be chained all within the general framework of ToolChain.
-To detect when a voice may be killed off, monitor its output after noteOff and if it stays below
 a certain threshold for some specified time, kill it off.



*/