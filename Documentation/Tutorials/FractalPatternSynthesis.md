 Turtle Graphics
 ==============================================================================
 
 Turtle Graphics is a term for a conceptual plotting systems where we imagine a virtual pen that responds to commands, thereby drawing a picture. For some less than obvious reason, this pen is envisioned as a turtle. The commands that this turtle can respond to involve at least a command for moving one unit forward while drawing a line, one for turning left and one for turning right. Many implementations also offer a command to move a step forward without drawing a line. These commands are normally single characters, typical is:
 
 * F Move a step forward while drawing a line
 * f Move forward without drawing a line
 * + Turn one angular unit to the left (counterclockwise)
 * - Turn one angular unit to the right (clockwise)

The step size for forward steps and the turning angle are parameters that are typically set up once and for all in advance. At any moment, the turtle has a state consisting of its current position and the direction it is heading towards. Initially, this position is typically the origin `(x,y) = (0,0)` and the direction is given by the x-axis `(dx,dy) = (1,0)`. Starting from this initial state and assuming the angle increment to be set to 90°, the turtle can now interpret a string, like `FF+F-FF` by drawing to units to right, turning to the left (heading upwards now), draw one line segment upward, turn to the right and draw another two line segments rightward.



 
 
 Lindenmayer Systems
 ==============================================================================
 
 
 
 Fractal Synthesis Language
 ==============================================================================
 
 To build a musical sound synthesis algorithms from the concepts above
 
 Extensions to Turtle Syntax
 ---------------------------
 
 
 Extensions to Lindenmayer Syntax
 --------------------------------
  
  
 Fractal Pattern Synthesis
 ==========================
 
 
 Curve Traversal
 ---------------
 
 ###Cyclic Resets (or not)
 
 
 ###Loop Modes
 
 
   
 
 Turning Angle Modulation ...maybe allow Step Size Modulation, too?
 
 

  4.3: Normalization Modes
  4.5: Subtractive Post Processing
   4.5.1 Regular Musical Filters
   4.5.2 Linear Phase Filters
  4.6: Rational Numbers and Harmonics
   4.5.2 Number of Segments (N) vs Segment Length (L)

   
   
   
Sound Design Guide:
-if you want the free-running behavior to match the reset-mode, you need to make sure that the
 turtle heads into the same direction as initially after it completed a cycle around the seed, i.e,
 for a square seed (at 90°) don't use F+F+F+F but F+F+F+F+, the additionla plus at the end makes 
 the turtle look to the right again after completing the square
-an turn angle sligtly off from the ideal value lets the picture slowly rotate in free-running mode
-rules with branches create nice overtone structures (branches introduce discontinuities in the 
 waveshape)
-to analyze a patch:
 -stop rotation (if any) by setting turning angle to precise value and/or use resetting
 -turn the numIterations to zero to see the seed
 -turn it to 1 to see the rule, maybe use a simple 'F' seed to see the pure rule
 
 
 just for refernece while writing this
 http://commonmark.org/help/
 https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet