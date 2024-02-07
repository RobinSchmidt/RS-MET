


/*=================================================================================================

Ideas:

- Implement another kind of diffuser based on an allpass FDN, see:
  https://arxiv.org/pdf/2007.07337.pdf  Allpass Feedback Delay Networks (Sebastian J. Schlecht)
  it has accompanying MatLab code: https://github.com/SebastianJiroSchlecht/fdnToolbox
  - For stereo-in / stereo out implementation, consider two different SISO implementations with a
    common (given) feedback matrix A, two common (given) input vectors B, two different 
    (given) output vectors C, and two different (to be computed) feedforward matrices D.
  - As feedback matrix A, choose a Hadamard-matrix combined with per-delayline scale factors to
    ensure uniform mode decay.
  - As injection matrix B, choose a matrix that injects the input signal with factors +-1 into
    all delaylines.
  - Or...wait - I think we may not be able to prescribe A,B,C at will and then compute D. I think,
    we need to ensure the the 3 conditions of Eq. 25 hold:
    A * A^T  + B * B^T = I, A * C^T + B * D^T = 0, C * C^T + D * D^T = I
    I think, when A is a Hadamard matrix times a diagonal matrix, then A * A^T simplifies to a 
    diagonal matrix: Therefore, B * B^T = I - A * A^T also becomes diagonal
    Maybe we can prescribe D. I think, it should be non-singular (not sure, though)
  - Maybe work out an example with 2 inputs, 3 outputs and 4 delaylines, such that it is nicely
    general. JOS-Filters, pg. 301 says there should be at least as many outputs as inputs (why?).
    And the numbers aren't too artifical either: 2 inputs occurs with stereo signals, 3 outputs 
    could mean a stereo + center setup and 4 delaylines seems to be a nice choice, too. A 4x4
    Hadamard matrix is also pretty nice.

- Maybe try a nested allpass structure (a la Gardner) instead of a series (a la Schroeder). 
  see: http://gdsp.hf.ntnu.no/lessons/6/33/  ..has a nice block diagram for an implementation 
  structure of the nested allpasses. Ah - no - it's three allpass filters in series nested in 
  one outer allpass.

- A nested allpass design fetaures an increasing echo density over time whereas a series 
  conncetion featues a constant echo density. 
  see: https://valhalladsp.wordpress.com/tag/nested-allpass/ 


- Maybe to make the sound of a reverb algo independent from the sample rate, we should use the
  same physical delay-times (in seconds) at all sample rates. That requires interpolating the 
  delaylines. Looks like a perfect task for allpass interpolation. Maybe do this for the 
  pre-diffuser. For the FDN (with the longer delay times), rounding to integer should work well
  enough. The coeffs of the pre-diffuser should the also be computed from desired decay times.
  but we need to take care: we cannot call delayLines[i].readOutput() twice per sample when we want
  to do allpass interpolation on the delayline. We must read it once and store it in another temp
  variable.


See also:
https://ccrma.stanford.edu/~jos/fp/Allpass_Filters.html
https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters_I.html



*/