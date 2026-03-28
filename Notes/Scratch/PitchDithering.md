
Background
-----------

A naively implemented digital oscillator produces a lot of aliasing. Various methods exists to mitigate the problem. Some of the methods are: mip-mapping, bleps and oversampling. This document describes yet another one of those methods. In my explanations of the method, I will take a sawtooth wvae as example but the method can be applied to other waveforms as well.  
...TBC...


Idea
----

When the length of the cycles that we want to produce happens to be an integer number of samples, the aliasing frequencies happen to line up with the harmonics that are already there. In this case, the presence of aliasing frequencies does not introduce any undesired additional frequencies into the signal but instead just changes the amplitudes of the existing harmonics. This is a much less annoying kind of artifact which in this method, we will accept. Of course, the problem now is that we can only produce sawtooths with those fundamental frequencies whose pitch period happens to be an integer number of samples. If we just round the cycle length to the nearest integer, we would get considerable mistuning which would get worse towards higher pitches.  
...TBC...








