
Background
-----------

A naively implemented digital oscillator produces a lot of aliasing. Various methods exists to mitigate the problem. Some of the methods are: mip-mapping, bleps and oversampling. This document describes yet another one of those methods. In my explanations of the method, I will take a sawtooth wvae as example but the method can be applied to other waveforms as well.  
...TBC...


Idea
----

When the length of the cycles that we want to produce happens to be an integer number of samples, the aliasing frequencies happen to line up with the harmonics that are already there. In this case, the presence of aliasing frequencies does not introduce any undesired additional frequencies into the signal but instead just changes the amplitudes of the existing harmonics. This is a much less annoying kind of artifact which in this method, we will accept. Of course, the problem now is that we can only produce sawtooths with those fundamental frequencies whose pitch period happens to be an integer number of samples. If we just round the cycle length to the nearest integer, we would get considerable mistuning which would get worse towards higher pitches. When we have a sampling rate of $f_s$ and we want to produce a frequency $f$, then realtion between the cycle length $c$ in samples and frequency $f$ is given by:

$$\boxed{c = \frac{f_s}{f}, \quad f = \frac{f_s}{c}}$$

For example, if we assume a sampling rate of $f_s = 44100$ Hz and a cycle length of $c = 100$ samples, we would get a frequency of $f = 441$ Hz. Let's now assume that we want to produce a sawtooth with a cycle length of $100.3$ samples. We can't produce cycles with the non-integer length of $c = 100.3$ but we can produce cycles of length $c_l = 100$ and we can also produce cycles of length $c_h = 101$ where the subscripts $l,h$ stand for "low" and "high". What if we probabilistically alternate between these two integer cycle lengths $c_l, c_h$ in such a way that the _average_ cycle length comes out as our desired $c$? 

...TBC...








