

//=================================================================================================
/*

ToDo:

- Drag in the code from FlatZapper and turn it into a class rsAllpassDisperser

- Maybe implement Thiran allpass interpolators.

- See also rosic::AllpassChain. IIRC, this was my first implementation that I used in my early
  phasers. It uses a DF1 biquad - which sucks! A state of art implementation should wap that out 
  for an SVF. The SVF should then also support 1st order allpass modes. And maybe also modes that
  combine 2 1st order allpasses into a 2nd order one. Maybe we should have an SVF-Chain that 
  implements a chain of a given number of equal SVFs. ...but actually, the limitation of the 
  filters to be all equal takes a lot of the potential fun away

- Implement creation of allpass filters from a given set of complex poles. The zeros are just the
  reciprocals of the poles (maybe with conjugation, not sure - but complex poles come pairwise 
  anyway). Use that to create Butterworth, Bessel, Papoulis, etc. allpases, i.e. allpases based on
  well known allpole lowpass designs. This is a general recipe for creating an allpass filter from
  any filter: just take its poles and use as zeros the reflected poles.

- Another general recipe (applicable to stable, minimum phase filters) could be: invert the filter 
  (swap numertaor and denominator), reverse the FIR part, take the product (i.e. series connection) 
  of the original and the invert-reversed one.

- Implement a phaser like allpass - maybe one in which all stages have the same coeffs.

- Maybe make versions of the classes that allow for fractional delaylines. Maybe that can be solved
  by templatizing on the delay-type (e.g. int vs double) and delayline type (basic, fractional, 
  etc.). Or maybe just allow that a fractional delayline can be used like an integer one, i.e. give
  it an "integer mode".

- Maybe add a module to ToolChain called AllpassZoo. Or maybe the different allpasses should be
  separate modules. But we may want to make a folder Allpass.

- Add getTransferFunctionAt functions. Do this also for the prototypes. Their implementations of it
  may look very different.


Interesting Resources:

  Frequency-Dependent Schroeder Allpass Filters (Sebastian Schlecht)
  https://www.mdpi.com/2076-3417/10/1/187

  ENERGY-PRESERVING TIME-VARYING SCHROEDER ALLPASS FILTERS (Kurt James Werner)
  https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_59.pdf



*/