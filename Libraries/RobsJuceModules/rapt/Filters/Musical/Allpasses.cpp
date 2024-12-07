

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



*/