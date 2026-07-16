
Pitch Dithering (Draft)
=======================

Background
----------

A naively implemented digital sawtooth oscillator produces a lot of aliasing. Various methods exist 
to mitigate the problem. Some of the methods are: mip-mapping, bleps and oversampling. This document
describes yet another one of those methods that I recently came up with. It's a way to replace the
annoying aliasing artifacts with a much more palatable kind of noise. In my explanations of the 
method, I will take a sawtooth wave as example but the method can be applied to other waveforms as 
well. In fact, if you scale and shift the output range of the sawtooth from -1...+1 to 0..1, you can
just use the sawtooth as "phasor" to produce any waveform from the saw and the result will be
likewise anti-aliased as well. As a corollary, you can even apply any waveshaper to the so produced
waveforms and the anti-aliasedness will "survive" the waveshaping which is a pretty unique feature
of this method.


The Initial Idea
----------------

When the length of the cycles that we want to produce happens to be an integer number of samples,
the aliasing frequencies happen to line up with the harmonics that are already there. In this case,
the resulting waveform will be perfectly periodic even in the discrete time sense such that 
`x[n] = x[n+P]` where P is the integer(!) period in samples. In this particular scenario, aliasing
is still present but the aliasing does not introduce any undesired additional frequencies into the
signal but instead just changes the amplitudes of the existing harmonics. This is a much less 
annoying kind of artifact which in this method, we will accept. Of course, the problem now is that
we can only produce sawtooths with those fundamental frequencies whose pitch period happens to be
an integer number of samples. If we just round the cycle length to the nearest integer, we would get
considerable mistuning which would get worse towards higher pitches. When we have a sampling rate of
$f_s$ and we want to produce a frequency $f$, the relation between the cycle length $c$ in samples
and frequency $f$ is given by:

$$\boxed{c = \frac{f_s}{f}, \quad f = \frac{f_s}{c}}$$

For example, if we assume a sampling rate of $f_s = 44100$ Hz and a cycle length of $c = 100$
samples, we would get a frequency of $f = 441$ Hz. Let's now assume that we want to produce a
sawtooth with a cycle length of $c = 100.3$ samples. We can't produce cycles with the non-integer
length of $100.3$ but we can produce cycles of length $c_1 = 100$ and we can also produce cycles of
length $c_2 = 101$. What if we probabilistically alternate between these two integer cycle lengths
$c_1, c_2$ in such a way that the _average_ cycle length comes out as our desired $c$? To achieve
that, we would have to produce cycles of length $c_1 = 100$ with a probability of $p_1 = 0.7$
and cycles of length $c_2 = 101$ with a probability of $c_2 = 0.3$. As a general rule, we could
always use two different cycle lengths $c_1 = floor(c)$, $c_2 = c_1 + 1$, $c_f = c - c_1$, 
$p_1 = 1 - c_f$, $p_2 = c_f$ where $c_f$ is the fractional part of $c$.

But there's a problem with this approach. With this rule as stated above, we would indeed always
produce an average cycle length that is exactly as prescribed. So we have solved the mistuning 
problem. But we have now introduced a new problem. It's arguably a less severe problem, so we 
actually did make some progress but it's still not good enough. Doing it like explained above does,
of course, produce some sort of artifacts. Namely, we introduce a sort of frequency modulation by a 
random pulse wave signal. This random frequency modulation manifests itself as a sort of noise in 
the final output. This noise in itself is something we are going to accept in this method. But what 
we don't want to accept is that the amount of this noise currently depends on the particular setting
of the desired cycle length $c$. If $c$ happens to be an exact integer, there will be no noise at 
all because the fractional part $c_f = 0$ is zero in this case and we will therefore produce cycles 
of length $c_1$ with probability $p_1 = 1$. Apparently, we will get the greatest amount of noise 
when $c$ happens to be halfway between two integers, i.e. $c = xxx.5$ and no noise at all when $c$ 
is an exact integer $c = xxx.0$. The amount of noise would vary as function of the fractional part 
of our desired cycle length. To have a consistent sound character of the oscillator, we don't want 
this. The amount of added noise should be the same regardless of how close to an integer our 
requested cycle length $c$ happens to be. We now want to equalize the noise, i.e. make it sound the 
same regardless of our value of $c_f$.


The Refined Idea
----------------

To develop a solution strategy, let's assume that our desired cycle length is $c = 100.0$. With the
basic algorithm above, we would get a clean signal with no noise modulation at all. We ask ourselves
how we would voluntarily introduce a noise into this signal that is statistically and sonically
similar to the noise in the $c = 100.5$ case. This is the worst case that would produce the greatest
amount of noise with the initial idea described above. It is "worst" because at the half integers, 
we are the farthest away possible from the "clean" case where $c$ is an exact integer. The new idea 
is now to also use random cycle lengths, even though we don't have to if the only goal would be to 
get the (average) cycle length right. Of course, we want to maintain an average cycle length of 
$100$. In order to achieve that, it is clear that we additionally need to use cycle lengths above 
_and_ below $100.0$. We need to use cycles of the 3 lengths $c_1 = 99, c_2 = 100, c_3 = 101$ in such
a way that the mean cycle length is still exactly $100$ and the variance of the probability 
distribution matches the variance that we would get in the worst case scenario, i.e. at the 
half-integers. It is apparent by now that the general task to make this work is to derive a formula 
or algorithm to compute the 3 desired cycle lengths $c_1, c_2, c_3$ along with their associated 
probabilities $p_1, p_2, p_3$ of producing cycles of these lengths. The input is the given desired 
mean cycle length $c$. As before, let $c_f = c - floor(c)$ denote the fractional part of our desired
(mean) cycle length $c$. If $c_f = 0.5$, we expect to be in an edge case where from the 3 lengths 
$c_1,c_2,c_3$ are only 2 actually used because one gets a probability of zero. This is our reference
case and we need to produce the values $c_1,c_2,c_3$ and $p_1,p_2,p_3$ for the other cases in such a
way, that the noise has always the same characteristics. We will use $c_2$ as our middle cycle 
length and we will always have $c_1 = c_2 - 1$ and $c_3 = c_2 + 1$. In the case where $c_f < 0.5$, 
we will need to use $c_2 = floor(c)$ and in the case where $c_f > 0.5$ we will need 
$c_2 = floor(c) + 1$. That this is right can most easily be understood from an example. If we have a
desired cycle length of $c = 100.5$ samples, we would use cycles of $100$ and $101$ samples with 
equal probability, namely with probability $0.5$. When the mean cycle length is lower, say 
$c = 100.3$, then we would expect to additionally also use cycles of length $99$ samples and when 
the mean cycle length is higher, say $100.7$, then we would additionally have to use cycles of 
length $102$. To summarize, for the 3 cycle lengths $c_1,c_2,c_3$ to be used, we use the following 
rule (in pseudocode):
```
ci = floor(c)         # Integer part of c
cf = c - ci           # Fractional part of c

if(cf < 0.5)
  c2 = ci             # Mid length when cf < 0.5
else
  c2 = ci + 1         # Mid length when cf >= 0.5

c1 = c2 - 1           # Short length
c3 = c2 + 1           # Long length
```
Now that we have determined the 3 cycle lengths $c_1,c_2,c_3$ to use, the next step is to determine
their associated probabilities $p_1,p_2,p_3$. To determine 3 values, we need 3 equations. The first
equation can be obtained from the requirement that our $p$ values have to add up to $1$ if we want
to interpret them as probabilities for 3 mutually exclusive events that together cover all the
possibilities, so we require: $p_1 + p_2 + p_3 = 1$. Next, we want to require that the mean cycle
length is our prescribed $c$, so we could use $p_1 c_1 + p_2 c_2 + p_3 c_3 = c$. However, for the
derivation, it turns out to be more convenient to express this equation in terms of the errors that
we make with our 3 cycle lengths. That is, we define the 3 errors $e_1 = c_1 - c, e_2 = c_2 - c,
e_3 = c_3 - c$ and require that the mean error is zero: $p_1 e_1 + p_2 e_2 + p_3 e_3 = 0$. The third
equation is obtained from our desire to always have the same variance $v$. The variance is the
expectation value of the squared errors, so we set $p_1 e_1^2 + p_2 e_2^2 + p_3 e_3^2 = v$. But what
value is that $v$? To figure that out, we turn again to our reference case where $c_f = 0.5$. In
that case, we know that we would only be dealing with two possible error values of $-0.5$ and $+0.5$
which both would occur with a probability of $0.5$. This gives the variance 
$v = 0.5 (-0.5)^2 + 0.5 (+0.5)^2 = 0.25$. So, our target value for $v$ is $1/4$. We can now give
these 3 equations to the computer algebra system SageMath using the following code:
```
var("e1 e2 e3 p1 p2 p3")
eq1 = 1   == p1       + p2       + p3
eq2 = 0   == p1*e1    + p2*e2    + p3*e3
eq3 = 1/4 == p1*e1*e1 + p2*e2*e2 + p3*e3*e3
solve([eq1,eq2,eq3],[p1,p2,p3])
```
which gives the result:
```
...something to do...
```




To summarize, the pseudocode to compute the 3
probabilities could look like:

...TBC...ToDo: Copy the solution formulas for the 3 probabilities $p_1, p_2, p_3$ from the code in
the research repo into here. Maybe also copy the derivation. Maybe try the pseudocode in Python. I
guess it could even work.




Implementation
--------------

A working implementation in C++ of an oscillator based on this idea can be found here:

https://github.com/RobinSchmidt/RS-MET/blob/work/Libraries/RobsJuceModules/rapt/Generators/PitchDitherOscs.h
https://github.com/RobinSchmidt/RS-MET/blob/work/Libraries/RobsJuceModules/rapt/Generators/PitchDitherOscs.cpp

ToDo: Give example code for how it can be used.


Experimental Results
--------------------

ToDo: Review the experimental results and maybe show some plots of spectra here. Maybe create audio
examples and link them here. Maybe implement an interactive example implementation using APE and
produce a little demo video with it and link to it here. Apply it to the supersaw. I think, it 
should be great for that because the introduced noise further thickens the spectrum and adds some
element of random pitch modulation to the indiviudal saws while also anti-aliasing them efficiently.
In fact, the idea for this method was born in a forum discussion about the JP-8000 supersaw here:  
https://www.kvraudio.com/forum/viewtopic.php?p=9189004#p9189004



Further Ideas
-------------

I think, the fact that we tune our probabilities in such a way to get the arithmetic(!) mean of the 
period right implies that in terms of frequencies, we hit the correct desired "mean frequency" only
when we interpret the "mean" as the harmonic mean. Maybe that is not exactly the right thing to do.
Maybe we should try to get the artithmetic mean frequency right which would imply that we would need
to get the harmonic mean of the periods right. So maybe the second equation (in its original form) 
$p_1 c_1 + p_2 c_2 + p_3 c_3 = c$ should be replaced by $p_1 / c_1 + p_2 / c_2 + p_3 / c_3 = 1 /c$?
Try that! Maybe then we should also define the 3 errors differently - namely as frequency errors
rather than period errors. This is basically a perceptual question: which frequency do we _perceive_
as the center frequency in a rapidly alternating jumble of frequencies? Maybe try to set up a
perceptual experiment to figure that out. Maybe we should also look into trying to fix the geometric
mean. See comments of this .md file (invisible in the rendered version) for more details.



<!---
<br><br><br><br><br><br><br><br><br><br>
----------------------------------------------------------------------------------------------------
Snippets

The half-integers are the worst case (why?). We volutarily introduce additional noise into the other
cases to match the noise in the worst case. Instead of using two cycles (except at exact integers
where we use one), we use three cycles (except at exact half-integres where we use two).

The derivation is in the file:
C:\Users\rob\GitRepos\RS-MET-Research\Notes\TempSketchPad.txt

ToDo:

- Explain what happens when we use a deterministic instead of a probabilistic algorithm to determine
  the next cycle length. In this case, the artifacts sound similar to aliasing.

- I think that using the geometric mean would mean to use c1^p1 * c2^p2 * c3^p3 = cbrt(c) as 2nd
  equation (i.e. in place of p1*c1 + p2*c2 + p3*c3 = c)

- Set up a perceptual experiment as follows: Alternate between cycles of length c1 = 100 and
  c2 = 200 and try to find a single length c that leads to the same pitch sensation. Will it be
  the arithmetic, geometric or harmonic mean of 100 and 200? That experiment will determine which
  of the means is the right one to match. Maybe the implementation should get other methods besides
  setMeanPeriod() namely: setMeanFrequency(), setMeanPitch(). We will need two other helper
  functions to calculate cycle distributions accoridng to the different rules. The part that
  computes the cycle lengths c1,c2,c3 will be the same everywhere so it should be factored out.

-->