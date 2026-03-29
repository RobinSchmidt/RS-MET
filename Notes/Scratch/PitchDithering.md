
Background
-----------

A naively implemented digital oscillator produces a lot of aliasing. Various methods exist to mitigate the problem. Some of the methods are: mip-mapping, bleps and oversampling. This document describes yet another one of those methods. In my explanations of the method, I will take a sawtooth wave as example but the method can be applied to other waveforms as well.  
...TBC...


The Basic Idea
--------------

When the length of the cycles that we want to produce happens to be an integer number of samples, the aliasing frequencies happen to line up with the harmonics that are already there. In this case, the presence of aliasing frequencies does not introduce any undesired additional frequencies into the signal but instead just changes the amplitudes of the existing harmonics. This is a much less annoying kind of artifact which in this method, we will accept. Of course, the problem now is that we can only produce sawtooths with those fundamental frequencies whose pitch period happens to be an integer number of samples. If we just round the cycle length to the nearest integer, we would get considerable mistuning which would get worse towards higher pitches. When we have a sampling rate of $f_s$ and we want to produce a frequency $f$, then realtion between the cycle length $c$ in samples and frequency $f$ is given by:

$$\boxed{c = \frac{f_s}{f}, \quad f = \frac{f_s}{c}}$$

For example, if we assume a sampling rate of $f_s = 44100$ Hz and a cycle length of $c = 100$ samples, we would get a frequency of $f = 441$ Hz. Let's now assume that we want to produce a sawtooth with a cycle length of $100.3$ samples. We can't produce cycles with the non-integer length of $c = 100.3$ but we can produce cycles of length $c_l = 100$ and we can also produce cycles of length $c_h = 101$ where the subscripts $l,h$ stand for "low" and "high". What if we probabilistically alternate between these two integer cycle lengths $c_l, c_h$ in such a way that the _average_ cycle length comes out as our desired $c$? To achieve that, we would have to produce cycles of length $c_l = 100$ with a probability of $p_l = 0.7$ and cycles of length $c_h = 101$ with a probability of $c_h = 0.3$. As a general rule, we could always use two different cycle lengths $c_l = floor(c)$, $c_h = c_l + 1$, $c_f = c - c_l$, $p_l = 1 - c_f$, $p_h = c_f$ where $c_f$ is the fractional part of $c$.

...TBC...




The Refined Idea
----------------

With this rule as stated above, we would indeed always produce an average cycle length that is exactly as prescribed. But we have now introduced a new problem. Doing it like explained above does, of course, produce some sort of artifacts. Namely, we introduce a sort of frequency modulation by a random pulse wave signal. This random frequency modulation manifests itself as a sort of noise in the final output. The amount of this noise will depend on the particular setting of the desired cycle length $c$. If $c$ happens to be an exact integer, there will be no noise at all because the fractional part $c_f = 0$ is zero in this case and we will therefore produce cycles of length $c_l$ with probability $p_l = 1$. Apparently, we will get the greatest amount of noise when $c$ happens to be halfway between two integers, i.e. $c = xxx.5$ and no noise at all when c is an exact integer $c = xxx.0$. To develop a solution strategy, let's assume that our desired cycle length is $c = 100.0$. With the basic algorithm above, we would get a clean signal with no noise modulation at all. The idea is now to use cycles of the 3 lengths $c_l = 99, c_m = 100, c_h = 101$, where the subscripts $l,m,h$ stand for "low", "mid", "high" in such a way that the mean cycle length is also exactly $100$ and the variance of the probability distribution matches the variance that we would get in the worst case scenario, i.e. at the half-integers. ...TBC...


Derivation of the Probabilities
-------------------------------

It is apparent by now that the main task to make this work is to derive a formula or algorithm to compute the 3 desired cycle lengths $c_l, c_m, c_h$ along with their associated probabilities $p_l, p_m, p_h$ of producing cycles of these lengths from the given desired mean cycle length $c$.

...TBC...ToDo: Copy the solution formulas for the 3 probabilities $p_l, p_m, p_h$ from the code in the research repo into here. Maybe also copy the derivation.




Implementation
--------------

ToDo: Clean up the implementation in the research repo and move it into RAPT and link it here.


Experimental Results
--------------------

ToDo: Review the experimental results and maybe show some plots of spectra here. Maybe create audio examples and link them here. Maybe implement an interactive example implemenation using APE and produce a little demo video with it and link to it here.







<br><br><br><br><br><br><br><br><br><br>
----------------------------------------------------------------------------------------------------
Snippets

The half-integers are the worst case (why?). We volutarily introduce additional noise into the other cases to match the noise in the worst case. Instead of using two cycles (except at exact integers where we use one), we use three cycles (except at exact half-integres where we use two).

The derivation is in the file:
C:\Users\rob\GitRepos\RS-MET-Research\Notes\TempSketchPad.txt