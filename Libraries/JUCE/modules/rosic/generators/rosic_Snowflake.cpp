Snowflake::Snowflake()
{
  // init to order 4 Koch snowflake:
  axiom = "F--F--F";
  clearRules();
  addRule('F', "F+F--F+F");
  numIterations = 4;
  updateWaveTable();
}

void Snowflake::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  incUpToDate = false;
  //updateIncrement();
}

void Snowflake::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  incUpToDate = false;
  //updateIncrement();
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
  tableUpToDate = false; 
  //updateWaveTable();
  // maybe don't update here, instead set a "dirty" flag, check it in getSample and update the
  // table there
}

void Snowflake::setAngle(double newAngle) 
{ 
  renderer.setAngle(newAngle); 
  tableUpToDate = false; 
  //updateWaveTable();
  // todo: make the angle modulatable - don't re-render the table but only render the L-system 
  // output string, reduce it to the relevant characters and render one point at a time instead of
  // pre-rendering a wavetable
}

void Snowflake::clearRules() 
{ 
  renderer.clearRules(); 
  tableUpToDate = false; 
}

void Snowflake::Snowflake::addRule(char input, const std::string& output) 
{ 
  renderer.addRule(input, output); 
  tableUpToDate = false; 
}

void Snowflake::setAxiom(const std::string& newAxiom) 
{ 
  axiom = newAxiom; 
  tableUpToDate = false; 
}

void Snowflake::updateWaveTable()
{
  renderer.render(axiom, numIterations, x, y);
  tableLength = (int)x.size()-1;
  tableUpToDate = true;
  updateIncrement();
  reset();
}

void Snowflake::reset()
{
  pos = 0;
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
  -maybe these values should be precomputes for various angles, for example 0,5,10,15,.. and 
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
  and denominator), assume /1 when a / is encountered with any number thereafter
-make rotation angle quantizable to a set of numbers that the user can enter, for example
 120,90,72,60,45,36,30,24,22.5,20,15,12,10
-have different loop modes: forward, backward, alternating (avoids jumps for non-closed curves)
-instead of having a single L-system, have a set of them, maybe named A,B,C,.., then instead of 
 setting a number of iterations, use a string AAABBCCA to mean: apply A 3 times, then B 2 times, 
 then  C 2 times then A once, allow also syntax like A3B2C2A1
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

-call the whole synthesis method Fractal Geometric Synthesis (FG-synthesis), the extended 
 Lindenmayer/Turtle grammar Fractal Definition Language (FDL)
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


see here for inspiration for new curves
https://www.youtube.com/watch?v=RU0wScIj36o&t=53s

*/