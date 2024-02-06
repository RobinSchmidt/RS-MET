


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

- Maybe try a nested allpass structure (a la Gardner) instead of a series (a la Schroeder). 
  see: http://gdsp.hf.ntnu.no/lessons/6/33/  ..has a nice block diagram for an implementation 
  structure of the nested allpasses. Ah - no - it's three allpass filters in series nested in 
  one outer allpass.

- A nested allpass design fetaures an increasing eacho density over time whereas a series 
  conncetion featues a constant echo density. 
  see: https://valhalladsp.wordpress.com/tag/nested-allpass/ 




See also:
https://ccrma.stanford.edu/~jos/fp/Allpass_Filters.html
https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters_I.html



*/