Fractal Pattern Synthesis 
=========================

 
Turtle Graphics
---------------
 
Turtle Graphics is a term for a conceptual plotting system where we imagine a virtual pen that responds to commands, thereby drawing a picture. For some less than obvious reason, this pen is envisioned as a turtle. The commands that this turtle can respond to involve at least a command for moving one unit forward while drawing a line, one for turning left and one for turning right. Many implementations also offer a command to move a step forward without drawing a line. These commands are normally single characters, typically `F`to move forward while drawing, `f` to move forward without drawing and `+` and `-` to turn one angular unit left or right respectively.

The step size for forward steps and the turning angle are parameters that are typically set up once and for all in advance. At any moment, the turtle has a state consisting of its current position and the direction it is heading towards. Initially, this position is typically the origin `(x,y) = (0,0)` and the direction is given by the x-axis `(dx,dy) = (1,0)`. Starting from this initial state and assuming the angle increment to be set to 90°, the turtle can now interpret a string of commands. For example, `FF+F-FF` would mean to  draw two units to the right, turn to the left (heading upwards now), draw one line segment upward, turn to the right and draw another two line segments rightward.

Another useful operation is to store the current state (position and direction), draw something, and then come back to the stored state and draw soemthing else. If this is supported, storing the state is typically invoked by the command `[` and restoring the turtle from most recently stored state that was not yet retrieved by `]`. In summary, a typical turtle graphics system will respond to the following commands:

 * `F` Move a step forward while drawing a line
 * `f` Move forward without drawing a line
 * `+` Turn one angular unit to the left (counterclockwise)
 * `-` Turn one angular unit to the right (clockwise)
 * `[` Save current state to stack
 * `]` Restore most recent, non restored state from stack
 
Depending on the concrete implementation, the drawing system may involve more commands, for example to increase or decrease the forward step size and turning angle, change line thickness, color, etc. It is also possible to have a turtle that moves around in 3 dimensions (it has probably left the beach and reached the water now).


Lindenmayer Systems
------------------- 

A Lindenmayer system, for short L-system, is a string rewriting system that can be used to generate strings that can be fed into a turtle graphics system in order to produce attractive drawings. Originally, these systems were conceived by theoretical botanist Aristid Lindenmayer to model the development of plants. An L-system consists of an initial string, known as "axiom" (I prefer to call it seed) and a set of replacement rules. A replacement rule, also know as production rule or just production, consists of a "predecessor" character and a "successor" string. For example, a rule: `F = F+F` would replace every occurrence of `F` inside a string to which it is applied by `F+F`. A Lindenmayer system applies the set of production rules to the initial string, then to the output of that application and so on - as often as desired. We can see that the number of `F`s doubles in each iteration which means, our initial string will grow exponentially with the number of iterations.


Generation of Fractals
----------------------

A turtle graphics interpreter and Lindenmayer systems are the two basic ingredients that we need to produce line drawings with self-similar ("fractal") features. It's intutively obvious that the results will show self-similarity, because in the iterative application of the L-system rules, smaller parts (typically `F`orward lines) are replaced by more complex structures, inside of which in the next iterations again the forward lines are replaced by the same structures and so on.



### Famous Examples



Koch curve, Hilbert curve, Peano Curve, Dragon curve


Sound Syntesis from Fractals
----------------------------

 
 
_______________________________________________________________________________
stuff below is just a skeleton and snippets
 
Fractal Synthesis Language
--------------------------
 
To build a musical sound synthesis algorithms from the concepts above
 
### Extensions to Turtle Syntax
 
 
### Extensions to Lindenmayer Syntax
  
  
Fractal Pattern Synthesis
-------------------------
 
 
### Curve Traversal
 
##### Cyclic Resets (or not)
 
 
##### Loop Modes
 
 

 
 
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