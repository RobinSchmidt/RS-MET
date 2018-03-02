Fractal Pattern Synthesis 
=========================

In this tutorial, I describe the application of a paradigm - which is widely known in the world of computer graphics and theoretical biology - to the synthesis of audio signals.

 
Turtle Graphics
---------------
 
Turtle Graphics is a term for a conceptual line drawing system where we imagine a virtual pen that responds to commands, thereby drawing a picture. For some less than obvious reason, this pen is envisioned as a turtle. The commands that this turtle can respond to involve at least a command for moving one unit forward while drawing a line, one for turning left and one for turning right. Many implementations also offer a command to move a step forward without drawing a line. These commands are normally single characters, typically `F` to move forward while drawing, `f` to move forward without drawing and `+`, `-` to turn one angular unit left or right respectively.

The step size for forward steps and the turning angle are parameters that are typically set up once and for all in advance. At any moment, the turtle has a state consisting of its current position and the direction it is heading towards. Initially, this position is typically the origin `(x,y) = (0,0)` and the direction is given by the x-axis `(dx,dy) = (1,0)`. Starting from this initial state, the turtle can now interpret a string of commands and leave a trail of lines. For example, assuming the angle increment being set to 90°, `FF+F-FF` would mean to  draw two units to the right, turn to the left (heading upwards now), draw one line segment upward, turn to the right and draw another two line segments rightward.

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

A Lindenmayer system, for short L-system, is a string rewriting system that was originally conceived by theoretical botanist Aristid Lindenmayer to model the development of plants. An L-system consists of an initial string, known as "axiom" (I prefer to call it seed) and a set of replacement rules. A replacement rule, also known as production rule or just production, consists of a "predecessor" character and a "successor" string. For example, a rule: `F = F+F`, when applied to a string, would replace every occurrence of `F` by `F+F`. A Lindenmayer system applies the set of production rules to the initial string, then to the output string of that application and so on - as often as desired. We can see that the number of `F`s doubles in each iteration. That implies that our string will grow exponentially with the number of iterations. All rules that can possibly be applied to a string are applied in parallel. This reflects the way in which plants grow by cell division and distinguishes an L-system from another widely used class of formal languages - namely Chomsky grammars - which apply one rule at a time in a single iteration.

For example, consider the simple L-system:

Axiom: A  
Rules: A = BA, B = AB

The string would develop as follows:

0. A
1. BA
2. ABBA (thank you for the music!)
2. BAABABBA
3. ABBABAABBAABABBA
4. BAABABBAABBABAABABBABAABBAABABBA

etc. With more complex sets of rules and/or more complex seeds, the resulting strings can get an interesting and complex internal structure.


Generation of Fractals
----------------------

A turtle graphics interpreter and Lindenmayer systems are the two basic ingredients that we need to produce line drawings with self-similar ("fractal") features. The strategy is - if you didn't already guess it - to use a Lindenmayer system to produce a string that contains turtle commands and throw that string at our turtle graphics interpreter. The turtle will respond only to the command characters and ignore all other characters in the string. It's intutively obvious that the results will show self-similarity, because in the iterative application of the L-system rules, smaller parts (typically `F`orward lines) are replaced by more complex structures, inside of which in the next iterations again the forward lines are replaced by the same structures and so on.


### Other Constructions of Fractals

You have probably seen pictures of the iconic Mandelbrot set and pehaps also of the related Julia sets. These belong to a different class of fractals which are generated by entirely different algorithms. In these other algorithms, you usually have a set of update rules that operate on numbers, like, for example:
```
a' = a^2 - b^2 + x
b' = 2 * a * b + y
```
where (a',b') denote the new, updated values and (a,b) the values before the update. You would start at an initial point (a0,b0) and iterate these rules again and again. When you do this, your point (a,b) may grow larger and larger without bounds. ...or it may not. If it will or will not explode depends on your choice of the initial point and also - in this case - the additive constants x,y. When you make a particular choice for (x,y), initialize (a,b) to (0,0) and iterate these rules and color all points black, for which the iteration doesn't explode and all points white for which it does explode, you will get a basic rendering of the Mandelbrot set. More sophisticated artistic renderings also take into account how fast the iteration diverges and choose a color according to that. As you can see, these kinds of algorithm do straightforwardly suggest a way to create images by directly coloring pixels according to the (non)divergence properties of the coordinates. But they do not produce any instructions for line drawings at all. Here, we are ultimately interested in sound synthesis and in this case we really want instructions for line drawing because those can straightforwardly be translated to audio signals. So, even though algorithms like one for the Mandelbrot set are great for visual art, we won't cosnider them here anymore.


### Examples

Below are few statdard examples:

Koch curve:

Angle: 60°
Axiom: F  
Rule:  F=F+F--F+  
The Koch curve is one of the earliest and most well known fractals. The constructione begin with a straight line segment, i.e. the seed is given by `F` and there's only a single rule that replaces each line segment by 4 line segments that go: right, up, down, right in a triangular shape with and angle of 60°

Sierpinski Triangle

Hilbert Curve:

The Hilbert curve starts by a U-shaped connenction of 3 lines with a turn angle of 90°. At each step

Dragon curve




Sound Syntesis from Fractals
----------------------------

Now that we understand how to make  fractal line drawings, the final step is to translate them to audio signals. To this end, we imagine that the turtle traverses the curve at uniform speed, thereby producing two functions of time x(t) and y(t). We simply use these two functions as our time domain signals for the left and right stereo channel. When feeding the so generated stereo signal into an oscilloscope in x/y mode, we will see the line drawing on its screen. We adjust the travesal speed such that the turtle walks through the curve once per cycle. When the turtle has reached the end of the curve, we loop back to the beginning and start anew. If we want to produce a note at 100Hz, our turtle would run through the curve 100 times per second. The jumping back to beginning at the end may introduce a discontinuity in the output signal unless the curve happens to have coinciding start and end points. With some care in setting up the L-system rules and seed, we can make sure that this is actually is the case. However, discontinuities in audio signals are not necessarily a bad thing either. In subtractive synthesis, the discontinuities in the sawtooth and square wave are precisely the features that give them their rich spectral content - so the discontinuities resulting from looping may actually be desirable. It depends on what kind of sound you want to design and it's good to have both options.


### Closed Loops and Free Running Turtles

First, we want to look at how we can make sure that the starting point of the curve coincides with the endpoint. One simple way to ensure that is to use a seed that describes itself a closed curve.






points to make:
-different constructions/curves with equal looking shapes (order of traversal matters)


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