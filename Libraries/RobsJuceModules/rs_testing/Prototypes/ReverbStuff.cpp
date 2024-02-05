


/*=================================================================================================

Ideas:

- Implement another kind of diffuser based on an allpass FDN, see:
  https://arxiv.org/pdf/2007.07337.pdf  Allpass Feedback Delay Networks (Sebastian J. Schlecht)
  it has accompanying MatLab code: https://github.com/SebastianJiroSchlecht/fdnToolbox

- Maybe try a nested allpass structure (a la Gardner) instead of a series (a la Schroeder). 
  see: http://gdsp.hf.ntnu.no/lessons/6/33/  ..has a nice block diagram for an implementation 
  structure of the nested allpasses. Ah - no - it's three allpass filters in series nested in 
  one outer allpass.

- A nested allpass design fetaures an increasing eacho density over time whereas a series 
  conncetion featues a constant echo density. 
  see: https://valhalladsp.wordpress.com/tag/nested-allpass/ 



*/