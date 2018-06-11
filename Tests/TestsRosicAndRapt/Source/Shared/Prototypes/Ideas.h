#ifndef Ideas_h
#define Ideas_h


/*
-continuous/probabilistic logic functions that take in real numbers a,b (probabilities):
 -not(a)    = 1-a
 -and(a, b) = a*b;
 -or( a, b) = 1 - (1-a)*(1-b) = a + b - a*b
  here, (1-a)*(1-b) is the probability that a and b are both false
  nand(a, b) = 1-a*b
 -xor(a, b) = or(and(a,not(b)),and(b,not(a))) = a + b - a*b(3-a-b+a*b)
 -maybe these functions can be extended to take a 3rd argument c for the correlation between 
  a and b
  -and(a,b,c) should reduce to the form above when c=0: and(a,b,0) = a*b
   when c = +1: and(a=1,b=0,+1) = and(a=0,b=1,+1) = 0, and(a=1,b=1,+1) = a = b
   when c = -1: and(a=1,b=1,-1) = and(a=0,b=0,-1) = 0
   ...maybe we can find more conditions that must hold and derive a general formula for and(a,b,c)
   maybe and(a,b,c) = a*b-c*xor(a,b) or  and(a,b,c) = a*b-c*a*b*xor(a,b)?
   -check, if the result is always in 0..1
   -run numeric simulations using data of correlated 0s and 1s produced by some algorithm and 
    check, if the results predicted by the formula match the simulation result
   -hmm.actually in leupold 2, p355, there is a formula: 
    P(A and B) = P(A|B) * P(B) = P(B|A) * P(A)...well, kinda obvious..so it seems like it's enough 
    to know either P(A|B) or P(B|A) along with P(A),P(B) to compute P(A and B)
  -when a formula for "and" is established, a formula for "or" can be derived from 
   or(a,b) = not(and(not(a),not(b)))
  -maybe it should also output two new correlation values for (a,result), (b,result) as side 
   results
  -maybe, as a further generalization, not only (symmetric) correlations can be used but 
   conditional probabilities...maybe a correlation can be seen as a pair of symmetric conditational
   probabilities, where C = P(A|B) = P(B|A) but in general, these two conditional probabilities
   can be different -> look up - if that's the case, maybe bayes rule can be used for and(a,b)?
  -consider an example: produce random number x in 0..1, let: A: x < 0.5, B: 0.3 < x < 0.6,
   -> P(A)=0.5, P(B)=0.3, P(A|B)=1/3, P(B|A)=2/5=0.4 -> what's the correlation C(A,B)?
 -maybe a class for such variables can be written, using operators "*" for "and" and "+" for "or"
  -can we then also derive meaningful formulas for "/" and "-"?
 -...maybe using p,q instead of a,b is a nicer convention

-3 valued logic: true,false,unknown
 T and U = U, T or U = T, F and U = F, F or U = U, U and U = U, U or U = U
 the rest stays the same (F and T = F etc.)
 -can (?) be modelled by probabilistic logic by using 0.5 as "unknown" - maybe by clamping unknown
  results to 0.5

-use that probabilistic logic in a cellular automaton - the new activation of each cell is a 
 logical combination of its and neighbour activations
 -read out the celluar automaton one cell per sample, after each cycle compute new cell activations
 -have a global "fuzziness" parameter that crossfades all cell-activtions between hard 0/1 and 0.5
 -maybe lowpass filter cell contents during readout
 -maybe use a 2D automaton with various readout rules (scan lines, diagonal, hilbert-curve, etc.)
 -needs on-the-fly sample-rate conversion - read out the automaton always at its "natural" rate and
  then convert to target sample-rate
 -maybe neighbour correlations can be computed on the fly and be fed into the and/or operations 
  (maybe scaled by a factor)
 -the cells could be initialized with random values, maybe with an optional "balancing" of 
  0s and 1s: just count the 1s, and if there are n too many, change n randomly chosen 1-cells to 0
  (or the other way around)


Graph synth:
-there's a set of N vertices all of which have: 
 -x,y position
 -a list of connected other vertices, this list may contain duplications, for example, vertex 2 
  could have the list [0, 1, 4, 1]
 -a state k that indicates at which point in the list we currently are
-the vertices are traversed with uniform speed, adjusted so as to give the desired frequency
-when we arrive at vertex i, we look at list[k] for the next vertex to visit and increment the 
 state k (with wraparound) - when we arrive at the vertetex the next time, we will head for 
 list[k+1]
-maybe the vertices themselves could move around in space (this needs to be updated only as often 
 as we visit a vertex)

-parametric equations:
 4-leaf clover: 
 r\left(\theta\right)=3\left(1+0.3\sin^2\left(4\theta\right)\right)\sqrt[4]{\cos^2\left(2\theta\right)}
 \left(\cos\left(t\right)\cdot r\left(t\right),\ \sin\left(t\right)\cdot r\left(t\right)\right)
 https://www.desmos.com/calculator/hgg7kzkswn

Implicit equation solver:
-finds a set of solutions to the implicit equation f(x,y) = 0, i.e. fills x,y arrays with solutions
 to this equation
-each solution is constructed from an initial point x0,y0 in such a way that the distance between
 x0,y0 and the point x,y on the curve that satisfies f(x,y) = 0 is minimized
-the initial points x0,y0 for filling the array can be constructed by some formula, for example, we 
 could use points on the unit circle or unit square
-define the square of the distance between (x,y) and (x0,y0) as: 
 d2(x,y) := (x-x0)^2 + (y-y0)^2
-partial derivatives of d2 with respect to x and y are given by: dx = 2(x-x0), dy = 2(y-y0)
-or solution must satisfy the equation: F(x,y) := f(x,y) + dx^2 + dy^2 = 0, or:
 F(x,y) = f(x,y) + 4*((x-x0)^2 + (y-y0)^2)
-use 2D Newton iteration to solve F(x,y) = 0 using x0, y0 as initial guess
-actually, the factor 4 in this equation could be any nonzero factor, maybe it makes sense
 check, how the factor influences the convergence of the iteration, maybe it could be an optional
 parameter to the routine that defaults to 4...or maybe some other experimentally determined value
 that most often leads to fastest convergence
-..oh - but Newton iteration needs a formula for the gradient of F(x,y)...maybe we need a numerical
 approximation of that...or maybe an entirely different 2D root-finding algorithm can be used

3D to 2D projection:
-user specifies 8 parameters (degrees of freedom):
 -camera postion in cartesian (x,y,z) or spherical (r,phi,theta) coordinates,
  maybe allow also cylinder (r,phi,z) coordinates
 -camera direction or the point looked at (object/model position), also as coordinate 
  triple (x,y,z)
 -camera rotation alpha (around the line connecting the camera position and object position)
 -zoom or viewing angle (sort of solid angle fo what is finally visible)
-from the 8 parameters, another set of 8 values is computed that can be used in a 3D-to-2D 
 transform realized by the equation:
 |X| = |xx xy xz| * |x| + |dx|
 |Y|   |yx yy yz|   |y|   |dy|
                    |z|
 where (X,Y) is the 2D output vector, (x,y,z) the 3D input vector and the other 8 values 
 coefficients to be determined from our 8 user parameters
-setting the camera distance to infinity should give a parallel projection ("orthographic view")
-to derive the equations for the coeffs, maybe go through the usual 4D homogenous coordinate 
 approach used in computer graphics and boil it down at the end...or maybe use affine transofrms 
 in 3D
-i think, the zoom and cam-rotation together determine the distance and orientation of the 
 projection-plane with respect to the camera
-see:
 http://www.cs.trinity.edu/~jhowland/cs3353/intro/intro/  (but this uses more parameters)
-maybe the equation above should be interpreted as reduced version of:
 |X|   |xx xy xz|   |x|   |dx|
 |Y| = |yx yy yz| * |y| + |dy|
 |Z|   |zx zy zz|   |z|   |dz|
 and just discarding Z afterwards. the full equations is just an affine transform in 3D - maybe
 an rsAffineTransform3D can be used as baseclass
-maybe to derive the equations, we need to define a projection plane, i.e. a plane whose normal
 vector is given by the camera's aiming direction which itself is given by the difference between
 camera position and aiming point - maybe to make the plane unique, choose the one that is halfway
 between these two vectors - or maybe just choose a plane that is unit-distance away from the 
 camera - figure out, what's most convenient
 

Equation synth ("Equator", "Formula..", "Solv..."):
-use a func-shape like input to let the user define signal shapes via equations
-user can select between parametric, implicit and differential equations
-maybe the parametric and implicit solvers can run in realtime (the differential solver must run in
 realtime anyway, i guess), so the user can manipulate in realtime some constants in his equations
 -the implicit solver could use initial points on a (phase-shaped) unit circle
 

Neural Network Synthesizer ("Cortex")
-use a delayline of length L where L is the length of one cycle of the signal
-use M read pointers, placed strategically in the delayline, for example at m*(L/M), M = 1,..,M
-the output at sample instant n is obtained as a function f(...) of the M delayline outputs. this
 function wll be realized by a neural network
-for a strictly periodic waveform, we could use M = 1 and use the delayline output with delay L 
 as is (our f(...) must be capable of being an identity function..could be difficual with neural 
 networks)
-in order to simulate damping from one cycle to the next, each read pointer could take a weigthed
 sum of 2 or 3 samples around delay L
-to build models, we can analyze samples recorded at various pitches, the pitch could be one of the
 network inputs, so we continuously interpolate between pitches
-maybe the time should also be a network input
-this may require an ordering of the neuron and/or synapases...or maybe obtain the model for the 
 highest pitch first (this should be the simplest model), then use the network weights as initial
 guess for the next lower pitch and so on
-to analyze a sample, maybe go through the cycles from right to left (because later cycles are 
 simpler than earlier ones)

-make a general Analysis -> Visualize/Edit -> Resynthesis GUI/module framework in TooChain which
 can be used for neural, modal, spectral, etc. modeling


-rewrite EngineersFilter
 -PrototypeDesigner:
  -has only static functions to fill arrays of poles/zeros/gain like
   butterworth(int order, complex* poles, complex* zeros, T* gain, G = 1, G0 = 0),
   bessel, elliptic, etc. - all with the same signature
   -if G=1, G0=0 they may delegate to simpler implementations that create "pass" prototypes
   -if there are more poles that zeros, zeros will be set to infinity (real part = +inf)
  -it only computes the non-redundant poles/zeros in the upper-left quarter-plane
  -poles and zeros are ordered by ascending frequency (imaginary part) and if poles happen to have
   the same frequency, by ascending real part (increasing Q)
 -PoleZeroMapper
  -static functions like sProtoToLowpass, sProtoToBandpass - maybe they should also all have the 
   same signature such that a higher level class can keep a function pointer to the mapper function
 -make class rsFilterDesignerAnalog for s-plane design, this one should have state-variables for
  method, mode, ripple, rejection, gain
 -make a class rsFilterDesignerDigital that uses and rsFilterDesignerAnalog and does all the bilinear 
  (or whatever other) mapping
 -maybe make an alternative rsFilterDesignerDigital2 that uses a z-plane prototype and does only a 
  z-to-z frequency transformation when the cutoff frequency changes (saves all the pre-warping 
  computations)
 -the new filter should also be able to use a SVF insstead of a biquad for each stage
 -figure out, which order of poles has the best modulation response - maybe it's best to 
  have low-Q poles first in the biquad array?



Scripting Synthesizer ("Synthax")


wavefolding:
http://www.kvraudio.com/forum/viewtopic.php?f=33&t=501471
out = 4.0 * (std::abs(0.25 * in + 0.25 - std::round(0.25 * in + 0.25)) - 0.25);


for performance tests:
-make the same test M times
-order results and throw away first m and last m to get rid of outliers, sample-size is now
 N = M-2m
-find mean and variance
-maybe plot histograms
-(maybe): when comparing results from different implementations of the same algorithm and one 
 turns out to be faster than the other on the average - test if this advantage is statistically
 significant - a bit like testing medication against placebos


*/

#endif