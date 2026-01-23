Snowflake::Snowflake()
{
  commandsReady = false;

  // init to order 4 Koch snowflake:
  turnAngle = 60;
  turtle.setAngle(turnAngle);
  axiom = "F--F--F--";
  clearRules();
  addRule('F', "F+F--F+F");
  numIterations = 4;

  updateTurtleCommands();
  updateWaveTable();
}

void Snowflake::clearRules() 
{ 
  lindSys.clearRules(); 
  commandsReady = false;
  tableUpToDate = false; 
}

void Snowflake::addRule(char input, const std::string& output) 
{ 
  lindSys.addRule(input, output);
  commandsReady = false;
  tableUpToDate = false; 
}

void Snowflake::setAxiom(const std::string& newAxiom) 
{ 
  axiom = newAxiom; 
  commandsReady = false;
  tableUpToDate = false; 
}

void Snowflake::setNumIterations(int newNumIterations) 
{ 
  if(newNumIterations == numIterations)
    return;
  numIterations = newNumIterations; 
  commandsReady = false;
  tableUpToDate = false; 
}

void Snowflake::updateAllInternals()
{
  if(!commandsReady) updateTurtleCommands();
  if(!incUpToDate)   updateIncrement();
  if(!tableUpToDate) updateWaveTable();
}

int Snowflake::getNumTurtleLines()
{
  if(!commandsReady) updateTurtleCommands();
  return numLines;
}

void Snowflake::updateTurtleCommands()
{
  lindenmayerResult = lindSys.apply(axiom, numIterations);
  setTurtleCommands(turtle.extractCommands(lindenmayerResult));
  commandsReady = true;
}

//=================================================================================================
/*

Bugs:

- Empty seed: crash (just delete the seed in InitSquare)...seems also to happen when there's no 'F'
  in the seed, I was trying to enter this example from the Book "Lindenmayer Systems, Fractals and 
  Plants"
  ? : ?X ? ?X
  p : X ? XFX ? ?XFX
  ...LSFP, page 73...there are more designs - harvest them empty seed bug is fixed


Ideas: 

- Allows an initiator to use different angles than the generator(s) - for example, start with a 
  pentagon and replace edges with something triangular (as in the koch snowflake)

- To make that more efficiently parsable in realtime, replace the + by P and the - by M in the 
  realtime-parsed string - avoids to check, if next char is a number whenever a + or - is 
  encountered, P and M indicate that is is, + and - indicate that it isn't, maybe allow rational
  numbers to be entered, like 20/3 for 6.66666° but also allow to use decimal dot (for numerator
  and denominator), assume /1 when a / is encountered without any number thereafter

- Instead of having a single L-system, have a set of them, maybe named A,B,C,.., then instead of 
  setting a number of iterations, use a string AAABBCCA to mean: apply A 3 times, then B 2 times, 
  then  C 2 times then A once, allow also syntax like A3B2C2A1, where the final 1 is optional
  ...maybe allow also things like (AABA)^2 (BBA^3)...or ((AABA)2(BBA)3)2

- On wikipedia, there's a Sierpinski triangle version that requires G to interpreted like F - this
  could be realized by defining:
  Axiom: F-G-G, rules: (A: F=F-G+F+G-F; G=GG) (B: G=F) and then doing AAAB (or A3B) for 3rd order
  ...maybe the parser doesn't have to care about the parentheses - it may separate the system
  definition by the colon, but they make it more readable...but maybe newlines could also be used
  A: F=F-G+F+G-F; G=GG
  B: G=F

- Maybe this best realized by making a class LindenmayerSystemSet.

- Allow the left hand sides of rules to be strings instead of characters. This allows a context 
  sensitive replacement, for example a rule +F+=+F-F-F+ means: Replace F by F-F-F, if the F is 
  between + and + (the surrounding "context" plusses appear in the rhs, too - so the rule won't 
  have them removed them in the output string), or +F+ = +FF+: "make a line twice as long, if it is
  between two turns" - but it's more general than context sesintive replacement of F, we could also
  replace it by +F+ = -FF-, for example

- Maybe rename the "Axiom" back to "Seed". 

- Maybe allow for a non-integer number of iterations, 3.7 would crossfade between the result from
  3 and 4 iterations with weights 0.3 and 0.7. Maybe we should precompute and store the strings
  for all orders up to something so we don't have recompute strings on the fly when the user 
  modulates numIterations.

- The whole generator could have different modes: "L-System" (or "Lindenmayer"), Finite Automaton, 
  etc. (the latter could split into "Christoffel Words", "Paper Folding", etc.)

- Add a post processing stage that takes the output string S from the L-system, automaton, etc. and
  applies a transformation like S = S+S+S+. Maybe with this, the pseudo-koch Curves in PGCA on 
  pages 183ff can be turned into closed triangular shaped loops (maybe this would belong into the 
  TurtleSource?)


Extension to 3D:

- User defines 3 rotation angles ax, ay, az where az corresponds to our current 2D angle. These
  should be modulatable

- The rotate left/right commands now additionaly rotate around the other axes as well.

- This way of implementing 3D does not require any extension of the turtle syntax and we don't need
  specific 3D L-systems (a different way would be to have extra commands for the other rotations 
  which would require such things). But I'm not sure, if it will produce good results.

- It may be musically desirable to have harmonic relationships between the angles, in particular
  ay = az/2, ax = az/4 or ax = az/3

- Maybe try matrices that are not pure rotations but also include scaling and shearing.

- Maybe even allow 4D matrices - the 4th w-coordinate can be interpreted for perspective.

- Maybe we could have an array of different triplets of angles and have further turtle commands to
  skip to the next/previous angle triple (mayb N, P)

- Maybe extend the L-System syntax to include turning commands in all 3 dimensions. Maybe use L/R
  for left/right (alternatively to +/-), I/O or into/out-of the screen (or F/B for 
  forward/backward - but no: F is already taken) and maybe U/D for up/down.

- Maybe for the 3D mode, we could use 2 independent L-Systems, let's call them A and B, both with 
  their own rules, seed, frequency, etc. They would produce 4 output signals which we denotes as
  A.x, A.y, B.x, B.y. Now we have a 4D signal which we must somehow proejct down to 2D, maybe by
  first projecting from 4D to 3D and then from 3D to 2D (at least conceptually - 
  implementationally, the 4D to 2D projection can be done in one single step). The coneptual split
  between the projections allows us to work with (modulatable) Euler angles in the 3D domain. Maybe 
  the 4D to 3D projection can be parametrized in terms of rotations in pq-planes where p and q 
  independently traverse x,y,z,w (i.e. A.x, A.y, B.x, B.y). Maybe we can have separate filters for
  the 4 signals.

- Maybe try to convert x,y (cartesian) coordinates to r,p (polar) coordinates. Then either draw the
  waveform interpreting rp as x,y or use the intermediate representation for filtering the signals
  and then converting back. Or: just interpret the L-system output as r,p and convert these to x,y.
  Or the other way around.



- Call the whole synthesis method Fractal Geometric Synthesis (FG-synthesis), the extended 
  Lindenmayer/Turtle grammar Fractal Definition Language (FDL) or maybe fractal geometric synthesis
  language (FGSL)...or maybe Fractal Pattern Synthesis

- Write a tutorial:
  1: Turtle Graphics
  2: Lindenmayer Systems
  3: Fractal Definition Language
  3.1: Extensions to Turtle Syntax
  3.2: Extensions to Lindenmayer Syntax
  4: Fractal Geometry Synthesis
  4.1: Turning Angle Modulation ...maybe allow Step Size Modulation, too?
  4.2: Loop Modes
  4.3: Normalization Modes
  4.5: Subtractive Post Processing
  4.5.1 Regular Musical Filters
  4.5.2 Linear Phase Filters
  4.6: Rational Numbers and Harmonics
  4.5.2 Number of Segments (N) vs Segment Length (L)
  See:
  https://github.com/RobinSchmidt/RS-MET/blob/work/Documentation/Tutorials/FractalPatternSynthesis.md


- Sound Design Guide:

  - If you want the free-running behavior to match the reset-mode, you need to make sure that the
    turtle heads into the same direction as initially after it completed a cycle around the seed, 
    i.e, for a square seed (at 90°) don't use F+F+F+F but F+F+F+F+, the additional plus at the end
    makes the turtle look to the right again after completing the square.

  - A turn angle sligtly off from the ideal value lets the picture slowly rotate in free-running 
    mode

  - Rules with branches create nice overtone structures (branches introduce discontinuities in the 
    waveshape)

  - To analyze a patch:

    - Stop rotation (if any) by setting turning angle to precise value and/or use resetting
    - Turn the numIterations to zero to see the seed
    - Turn it to 1 to see the rule, maybe use a simple 'F' seed to see the pure rule



see here for inspiration for new curves
https://www.youtube.com/watch?v=RU0wScIj36o&t=53s

https://en.wikipedia.org/wiki/List_of_fractals_by_Hausdorff_dimension
https://en.wikipedia.org/wiki/N-flake#Hexaflake

maybe, these are created using using parameteric equations?
http://benice-equation.blogspot.de/
http://benice-equation.blogspot.de/2012/01/fractal-spirograph.html
https://medium.com/@Alikayaspor/essential-mathematical-gifs-that-will-make-mathematics-finally-make-sense-4873573f5883

more math fun:
https://medium.com/@Alikayaspor/essential-mathematical-gifs-that-will-make-mathematics-finally-make-sense-4873573f5883


Some open source fractal drawing software:
https://larryriddle.agnesscott.org/ifs/software/software.htm
https://github.com/mathriddle/IFS-Lsystem
The golden dragon curve looks nice: 
https://larryriddle.agnesscott.org/ifs/heighway/goldenDragon.htm
https://pythonturtle.academy/golden-dragon-curve-fractal-source-code/

*/