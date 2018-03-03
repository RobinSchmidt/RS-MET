Fractal Pattern Synthesis 
=========================

In this tutorial, I describe the application of a generative algorithm which is widely known in the world of computer graphics and theoretical biology to the synthesis of audio signals. Everybody loves images of fractals for their visual beauty. It's about time to harness the underlying techniques to create beautiful sounds as well.

 
Turtle Graphics
---------------
 
Turtle Graphics is a term for a conceptual line drawing system where we imagine a virtual pen that responds to commands, thereby drawing a picture. For some less than obvious reason, this pen is envisioned as a turtle. The commands that this turtle can respond to involve at least a command for moving one unit forward while drawing a line, one for turning left and one for turning right. Many implementations also offer a command to move a step forward without drawing a line. These commands are normally single characters, typically `F` to move forward while drawing, `f` to move forward without drawing and `+`, `-` to turn one angular unit left or right respectively.

The step size for forward steps and the turning angle are parameters that are typically set up once and for all in advance. At any moment, the turtle has a state consisting of its current position and the direction it is heading towards. Initially, this position is typically the origin `(x,y) = (0,0)` and the direction is given by the x-axis `(dx,dy) = (1,0)`. Starting from this initial state, the turtle can now interpret a string of commands and leave a trail of lines. For example, assuming the angle increment being set to 90°, `FF+F-FF` would mean to  draw two units to the right, turn to the left (heading upwards now), draw one line segment upward, turn to the right and draw another two line segments rightward.

Another useful operation is to store the current state (position and direction), draw something, and then come back to the stored state and draw something else. Storing the state is typically invoked by the command `[` and restoring the turtle from most recently stored state that was not yet retrieved by `]`. In summary, a typical turtle graphics system will respond to the following commands:

 * `F` Move a step forward while drawing a line
 * `f` Move forward without drawing a line
 * `+` Turn one angular unit to the left (counterclockwise)
 * `-` Turn one angular unit to the right (clockwise)
 * `[` Save current state to stack
 * `]` Restore most recent, non restored state from stack
 
Depending on the concrete implementation, the drawing system may involve more commands, for example to increase or decrease the forward step size or turning angle, change drawing attributes like line thickness, color, etc. It is also possible to have a turtle that moves around in 3 dimensions instead of just 2 (it has probably passed the beach and reached the water now). 


Lindenmayer Systems
------------------- 

A Lindenmayer system, for short L-system, is a string rewriting system that was originally conceived by theoretical botanist Aristid Lindenmayer to model the development of plants. Besides their use as a scientific tool, they are also used in computer graphics to generate beautiful, natural looking models of plants and other organic structures. An L-system consists of an initial string, known as "axiom" (I prefer to call it seed) and a set of replacement rules. A replacement rule, also known as production rule or just production, consists of a "predecessor" character and a "successor" string. For example, a rule: `F = F+F`, when applied to a string, would replace every occurrence of `F` by `F+F`. A L-system applies the set of production rules to the initial string, then to the output string of that application and so on - as often as desired. We can see that with the simple example rule `F = F+F` the number of `F`s doubles in each iteration (assuming that this is the only rule that exists). That implies that our string will grow exponentially with the number of iterations. All rules that can possibly be applied to a string are applied in parallel. This reflects the way in which plants grow by cell division and distinguishes an L-system from another widely used class of formal languages - namely Chomsky grammars - which apply one rule at a time in a single iteration.

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

A turtle graphics interpreter and Lindenmayer systems are the two basic ingredients that we need to produce line drawings with self-similar ("fractal") features. The strategy is - if you didn't already guess it - to use a Lindenmayer system to produce a string that contains turtle commands and throw that string at our turtle graphics interpreter. The turtle will respond only to the command characters and ignore all other characters in the string. It's intutively obvious that the results will show self-similarity, because in the iterative application of the L-system rules, smaller parts (for example `F`orward lines) are replaced by more complex structures, inside of which in the next iterations again the forward lines are replaced by the same structures and so on. When a turtle interpreter is used to interpret an L-system output string, some implementations support the use of `G` and alternative to `F` as forward drawing command. The reason is that it may be useful to distinguish between edges that should be replaced by a production rule and others that shouldn't be replaced.


### Other Constructions of Fractals

You have probably seen pictures of the iconic Mandelbrot set and perhaps also of the related Julia sets. These belong to a different class of fractals which are generated by entirely different algorithms. In these other algorithms, you usually have a set of update rules that operate on numbers, like, for example:
```
a' = a^2 - b^2 + x
b' = 2 * a * b + y
```
where (a',b') denote the new, updated values and (a,b) the values before the update. You would start at an initial point (a0,b0) and iterate these rules again and again. When you do this, your point (a,b) may grow larger and larger without bounds. ...or it may not. If it will or will not explode depends on your choice of the initial point and also - in this case - the additive constants x,y. When you make a particular choice for (x,y), initialize (a,b) to (0,0) and iterate these rules and color all points (x,y) black, for which the iteration doesn't explode and all points white for which it does explode, you will get a basic black-and-white rendering of the Mandelbrot set. More sophisticated artistic renderings also take into account how fast the iteration diverges and choose a color according to that - or use even more sophisticated coloring algorithms taking into accout the whole trajectory or whatever. Anyway, as you can see, these kinds of algorithm do straightforwardly suggest a way to create images by directly coloring pixels according to the (non)divergence properties of the point at the pixel's coordinates. But they do not produce any instructions for line drawings at all. Here, we are ultimately interested in sound synthesis and in this case we really want instructions for line drawing because those can straightforwardly be translated to audio signals. So, even though algorithms like one for the Mandelbrot set are great for visual art, I not yet discovered any way to make them useful for sound synthesis - so we won't consider them here anymore. ...but for grahpics, they are *really* great - maybe better than this line drawing based stuff (if you want to learn more about visual fractal art, check out this amazing software: https://www.ultrafractal.com/).


### Examples

Below are few standard examples:

Koch curve:

Angle: 60°  
Axiom: F  
Rule:  F=F+F--F+  
The Koch curve is one of the earliest and most well known fractals. The construction begins with a straight line segment, i.e. the seed is given by `F` and there's only a single rule that replaces each line segment by 4 line segments that go: right, ramp-up, ramp-down, right in a triangular shape with and angle of 60°

Sierpinski Triangle (with `G`), without `G`

Hilbert Curve:

The Hilbert curve starts by a U-shaped connection of 3 lines with a turn angle of 90°. At each step

Dragon curve




Sound Synthesis from Fractals
-----------------------------

Now that we understand how to make  fractal line drawings, the final step is to translate them to audio signals. To this end, we imagine that the turtle traverses the curve at uniform speed, thereby producing two functions of time x(t) and y(t). We simply use these two functions as our time domain signals for the left and right stereo channel. When feeding the so generated stereo signal into an oscilloscope in x/y mode, we will see the line drawing on its screen. We adjust the traversal speed such that the turtle walks through the curve once per cycle. When the turtle has reached the end of the curve, we loop back to the beginning and start anew. If we want to produce a note at 100Hz, our turtle would run through the curve 100 times per second. The jumping back to beginning at the end may introduce a discontinuity in the output signal unless the curve happens to have coinciding start and end points. With some care in setting up the L-system rules and seed, we can make sure that this is actually is the case. However, discontinuities in audio signals are not necessarily a bad thing either. In subtractive synthesis, the discontinuities in the sawtooth and square wave are precisely the features that give them their rich spectral content - so the discontinuities resulting from looping may actually be desirable. It depends on what kind of sound you want to design and it's good to have both options.


### Closed Loops and Free Running Turtles

First, we want to look at how we can make sure that the starting point of the curve coincides with the endpoint. One simple way to ensure that is to use a seed that itself describes a closed curve. Assuming that the turn angle is set to 90°, a command string like `F+F+F+F` would draw a square. The endpoint of the turtle would coincide with the start point, so we have a closed curve. There's a little catch here, when we use this for sound synthesis. In the creation of graphical drawings, we usually trace out the curve just once - so in this setting, the command string above would be good enough to draw a simple square. But in sound synthesis, we'll trace out the curve again and again. Note that, after tracing out the curve, the turtle will be looking into the downward direction instead of rightward, as it initially did. So if we let the turtle run a second time to the command string without resetting its state before, we'll draw a second square below the first. In the third run, we'll draw a third square left to the second and in the fourth run we'll draw a fourth square on top of the third. Only after these four runs, the turtle will be again in its initial position *and* look into its initial direction. This can be useful, but you should be aware of it. If you want to trace out the same initial square again and again, one way to fix it would be to reset the turtle's state to the initial state after each run through the curve. Another way that doesn't require resetting would be to append another `+` to our command string, so it would become `F+F+F+F+`. This added `+` fixes the direction of the turtle after running to the command string. You need to watch out for this when using L-system definitions that you grab somehwere from the web or a book. There, a graphical setting is assumed in which the curve is traced out just once anyway, so there's no consideration of this effect and final "direction fixes" are usually absent.


### Resetting the Turtle

As mentioned in the previous section, it may be useful to reset the turtle's state after each run through the curve. But this is only the simplest strategy to reset the state and there's much more to say about it. In fact, I consider the resetting strategy as one of the major sound design parameters of this synthesis method.






