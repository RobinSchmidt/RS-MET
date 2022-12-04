template<class TSig, class TPar>
rsFormantRemover<TSig, TPar>::rsFormantRemover(int newMaxOrder)
  : rsLinearPredictor<TSig, TPar>(newMaxOrder)
{
  reset();
}

template<class TSig, class TPar>
void rsFormantRemover<TSig, TPar>::reset()
{
  coeff   = 0.0;
  pastIn  = 0.0;
  pastOut = 0.0;
  rsLinearPredictor<TSig, TPar>::reset();
}


/*

Ideas:

Maybe try to implement a spectral dynamics algorithm based on adaptive filtering. See:

https://www.kvraudio.com/forum/viewtopic.php?t=590275
https://oeksound.com/plugins/soothe2/

Maybe use a filtered input signal to drive the analysis while still acting on the original input. 
The filter can be lowpass/peak/bandpass/highpass etc. to make the algorithm more sensivtive to 
certain frequency ranges. Maybe that's what soothe does with these EQ-style filter settings. We 
need to run the LMS algo on a pre-filtered signal and then run the resulting filter also on the 
original input. It's like side-chaining where the side-chain signal is a pre-filtered input. 
Attack/Release parameters can be related to learn- and forget-rate. Try using more advanced 
adaptive filter structures like recursive least squares (RLS), gradient adaptive lattice (GAL), 
least squares lattice (LSL), fast transversal filter (FTF), maybe even adaptive IIR filters.




*/