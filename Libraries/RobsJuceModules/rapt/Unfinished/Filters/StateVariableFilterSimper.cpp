

//=================================================================================================
/*

Notes:

- On page 6, there's also the formula "k = 1/Q = 2 - 2*res". Does that mean there could be a 
  "resonance" parameter instead of Q? Figure out!


ToDo:

- Figure out and document what the coefficients and intermediate variables mean. Looking at 
  the scribble on the front page on the paper, it seem like k is the feedback factor after the
  1st integrator? And the a1, a2 are the gains of the integrators? And g is affecting them?
  v0 is the input voltage, v1, v2 the voltages after the 1st and 2nd integrator stage 
  representing bandpass and lowpass output? This can be inferred from the mixing coeffs 
  m0,m1,m2. They are 0,1,0 for bandpass and 0,0,1 for lowpass.

- Maybe have two template parameters TSig, TPar as in the other filters. I think,
  v0,v1,v2,v3,ic1eq,ic2eq must all be TSig, a1,a2,a3,m1,m2,m3 must be TPar

- Maybe move implementation into .cpp file ...but maybe not.

- Figure out the z-domain transfer function and implement a function 
  getTransferFunctionAt(rsComplex<TPar> z)

- Adjust the order of the modes in the enum to be the same as in the RBJ filter [DONE]...but RBJ
  has two bandpass variants. I think, this here is a const skirt gain bandpass. Maybe to obtain 
  const peak gain behavior, we just need to scale by k = 1/Q? [YES - DONE - seems OK]. The RBJ 
  filters are also missing a "peak" filter in the sense meant here. I think, it's just a 
  resonator? If so, try to introduce it in the RBJ filters as well. Maybe rename the mode to
  "Reson" or "Resonator". Maybe rename "Notch" to "Bandreject" for consistency. Or maybe call them
  Bandstop everywhere.  https://en.wikipedia.org/wiki/Band-stop_filter

*/