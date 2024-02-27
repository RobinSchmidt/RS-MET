A pack of samples that can be used as electronic bassdrums or, depending on the parameters, for 
"zap" sounds. The samples are all based on impulse responses of a specifically tuned chain of 
allpass filters. This means, the initial raw material has a *perfectly* white spectrum because the 
impulse itself is perfectly white and an allpass has a white magnitude spectrum by definition. The 
idea for this was inpired by a statement of Ady Scorb in this video at 9:57:

  https://www.youtube.com/watch?v=n3tmeChr7ec&t=9m57s

  "The philosophy I've come up with is: a kickdrum is like an impulse smeared in time. You want to 
  represent all the frequencies smoothly all the way down. This has the benefit of exciting every 
  frequency in a space and therefore you'll gonna get less nulls, less peaks and better translation 
  and better power in more spaces."

The technique of feeding an impulse into an allpass filter achieves exactly that - smearing an 
impulse in time - in a mathematically perfect way. No frequencies are boosted or attenuated. A bit 
of post-processing has been applied, though: the initially white output of the allpass chain has 
been filtered by a first order lowpass to give the spectrum a -6 dB/oct slope. Then, some subsonic 
rumble has been removed by a DC blocker highpass. Then a smooth fade-out was applied. 

The samples are named using the following scheme:

  ZappyKick_NS=50_FS=-2.wav

where the NS=50 means that the (N)umber of allpass (S)tages was set to 50 and some (F)requency 
(S)hape parameter was set to -2. The latter parameter controls, how the tuning frequencies of all 
the allpasses are distributed along the frequency axis. The allpass frequencies are distributed 
between 15 Hz and 8 kHz. These corner frequencies were fixed for this sample pack. For FS=0, the
allpass frequencies are distributed equidistantly in the pitch (i.e. log frequency) domain, for 
positive FS, there will be an upward (concave) curve, for negative FS a downward (convex) curve.

An interesting observation is that those samples with FS=0 feature a self-similar shape in the time 
domain. If you zoom in to the initial section, the samples will look basically the same as when you 
are zoomed out. This implies that it, to some extent, doesn't really matter if you play the sample 
back at its original speed or speed it up or slow it down - as samplers do, when using the same 
sample for different keys. When you do this with the FS=0 samples, the sound will not change. At 
least not very much. The residual changes are due to the bandwidth limitations of sampling and the 
two corner frequencies of the allpass chain. These limitations thwart perfect self similarity - but 
the self similarity is close enough to perfect such that you won't hear a big change for moderate 
playback speed changes.

Although the samples can be used as is, they are mainly meant to be used as raw material for 
further shaping. Out of the box, the samples as such are completely atonal. For example, if you like 
tuned bassdrums, you may want to process them further by perhaps boosting some frequency around your
fundamental with a peaking filter or maybe a resonant highpass. Or you may want to apply some 
distortion or layer them with some click samples. What is also interesting is to convolve them with 
a short noise burst.

The samples were created algorithmically by Robin Schmidt and are released under the Creative 
Commons CC BY-SA 4.0 License. See:

  https://creativecommons.org/licenses/by-sa/4.0/

The pack can be found as a 7-zip archive here:

  www.rs-met.com/sounds/samples/ZappyKicks.7z