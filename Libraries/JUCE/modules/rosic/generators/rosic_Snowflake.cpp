Snowflake::Snowflake()
{
  // init to order 4 Koch snowflake:
  axiom = "F--F--F";
  clearRules();
  addRule('F', "F+F--F+F");
  numIterations = 4;
  updateTurtleCommands();
  updateWaveTable();

  // tableLength = 768, numPoints = 769 - is this right? hmm - the last point in the table doesn't 
  // count so maybe yes
}

void Snowflake::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  incUpToDate = false;
}

void Snowflake::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  incUpToDate = false;
}

void Snowflake::setRotation(double newRotation)
{
  rotator.setAngle((PI/180) * newRotation);
}

void Snowflake::setNumIterations(int newNumIterations) 
{ 
  if(newNumIterations == numIterations)
    return;
  numIterations = newNumIterations; 
  commandsReady = false;
  tableUpToDate = false; 
}

void Snowflake::setAngle(double newAngle) 
{ 
  turtle.setAngle(newAngle);

  renderer.setAngle(newAngle); // old
  tableUpToDate = false; 


  // todo: make the angle modulatable - don't re-render the table but only render the L-system 
  // output string, reduce it to the relevant characters and render one point at a time instead of
  // pre-rendering a wavetable
}

void Snowflake::clearRules() 
{ 
  lindSys.clearRules(); 
  commandsReady = false;

  renderer.clearRules(); 
  tableUpToDate = false; 
}

void Snowflake::Snowflake::addRule(char input, const std::string& output) 
{ 
  lindSys.addRule(input, output);
  commandsReady = false;

  renderer.addRule(input, output); 
  tableUpToDate = false; 
}

void Snowflake::setAxiom(const std::string& newAxiom) 
{ 
  axiom = newAxiom; 
  // turtleCommandsUpToDate = false;
  commandsReady = false;
  tableUpToDate = false; 
}

void Snowflake::updateWaveTable()
{
  renderer.render(axiom, numIterations, tableX, tableY);
  tableLength = (int)tableX.size()-1;
  tableUpToDate = true;
  updateIncrement();
  reset();
}

void Snowflake::updateTurtleCommands()
{
  lindenmayerResult = lindSys.apply(axiom, numIterations);
  turtleCommands = turtle.extractCommands(lindenmayerResult);
  numPoints = turtle.getNumberOfPoints(turtleCommands);
  commandsReady = true;
  updateIncrement();
  reset();
}

void Snowflake::reset()
{
  turtle.init(0, 0, 1, 0);

  pos = 0; // used in table-based synthesis

  // used in on-the-fly synthesis:
  commandIndex = 0;
  x[0] = 0;
  y[0] = 0;
  getNextPoint(&x[1], &y[1]);
  //turtle.evaluateCommand(turtleCommands, 0
}

void Snowflake::getNextPoint(double* x, double* y)
{
  if(turtleCommands.size() == 0)
    return;
  // maybe check also if numPoints >= 1 (i.e. there is at least one 'F' - because if not, we may
  // hang up in an infinite loop
  

  int i = commandIndex;
  while(! (turtleCommands[i] == 'F' || turtleCommands[i] == 'f') )
  {
    i++;
    if(i >= turtleCommands.size())
      i = 0;
  }

  // an 'F' or 'f' was found at i

  commandIndex = i;

  // maybe return true or false dependning on whether a line should or shouldn't be drawn (i.e. an
  // 'F' or 'f' is encountered)
}

void Snowflake::updateIncrement()
{
  inc = tableLength * frequency / sampleRate; // avoid division
  incUpToDate = true;
}

/*
Ideas: 

-make TurningAngle modulatable
 -optionally compensate for changes in turning-angle by a rotation of the whole image, i.e. rotate
  image by -TurningAngle ..or maybe -compensationAmount*TurningAngle
 -requires to not pre-render a wavetable but only render the L-system output string and interpret 
  it on the fly (the string should be stripped from meaningless characters, retaining only F,+,- 
  before)
 -may require to compute shift and scale (to remove DC and normalize) when a new angle is set up
  -maybe these values should be precomputed for various angles, for example 0,5,10,15,.. and 
   interpolation be used (to avoid to expensively compute the actually desired value whenever the
   angle changes)
-allow a turtle syntax like +30 or -45 to mean: turn 30° left or 45° right (instead of whatever 
 value the turtle drawer is set to)
 -allows an initiator to use different angles than the generator(s) - for example, start with a 
  pentagon and replace edges with something triangular (as in the koch snowflake)
 -to make that more efficiently parsable in realtime, replace the + by P and the - by M in the 
  realtime-parsed string - avoids to check, if next char is a number whenever a + or - is 
  encountered, P and M indicate that is is, + and - indicate that it isn't, maybe allow rational
  numbers to be entered, like 20/3 for 6.66666° but also allow to use decimal dot (for numerator
  and denominator), assume /1 when a / is encountered without any number thereafter
-make rotation angle quantizable to a set of numbers that the user can enter, for example
 120,90,72,60,45,36,30,24,22.5,20,15,12,10
-have different loop modes: forward, backward, alternating (avoids jumps for non-closed curves)
 -maybe in backward passes, (optionally) interpret angles as their negative (necessarry to actually
  "go back")
-have a "reset" parameter that determines the number of cycles (through the command list) after 
 which the turtle state is reset into its initial state, interpret "0" as "never" - relevant for
 on-the-fly rendering with "weird" turn angles
-maybe have an additional string of turtle commands that can be applied after each cycle has passed
 (empty by default)...or maybe even different command strings to be applied after different numbers 
 of cycles ..this defines the "turn/wrap around behavior"
-instead of having a single L-system, have a set of them, maybe named A,B,C,.., then instead of 
 setting a number of iterations, use a string AAABBCCA to mean: apply A 3 times, then B 2 times, 
 then  C 2 times then A once, allow also syntax like A3B2C2A1, where the final 1 is optional
 ...maybe allow also things like (AABA)^2 (BBA^3)...or ((AABA)2(BBA)3)2
 -on wikipedia, there's a sierpinski triangle version that requires G to interpreted like F - this
  could be realized by defining:
  axiom: F-G-G, rules: (A: F=F-G+F+G-F; G=GG) (B: G=F) and then doing AAAB (or A3B) for 3rd order
  ...maybe the parser doesn't have to care about the parentheses - it may separate the system
  definition by the colon, but they make it more readable...but maybe newlines could also be used
  A: F=F-G+F+G-F; G=GG
  B: G=F
 -maybe this best realized by making a class LindenmayerSystemSet
-allow the left hand sides of rules to be strings instead of characters. this allows a context 
 sensitive replacement, for example a rule +F+=+F-F-F+ means: replace F by F-F-F, if the F is 
 between + and + (the surrounding "context" plusses appear in the rhs, too - so the rule won't 
 have them removed them in the output string)
-maybe rename the "Axiom" back to "Seed"
-interpolation modes: left, right, nearest, linear, cubic
-allow a traversal speed (or better: duration/length) to be associated with each point (i.e. the 
 line segment that goes toward that point)...like  a normal 'F' would assume 1, an 'f' would 
 assume '0', but we could also have F2 to let it take 2 time units to traverse the segment
 the incremrent computation would then not use "numPoints" as multiplier but 
 sum-of-traversal-durations

-call the whole synthesis method Fractal Geometric Synthesis (FG-synthesis), the extended 
 Lindenmayer/Turtle grammar Fractal Definition Language (FDL) or maybe fractal geometric synthesis
 language (FGSL)
-write a tutorial:
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


see here for inspiration for new curves
https://www.youtube.com/watch?v=RU0wScIj36o&t=53s

*/