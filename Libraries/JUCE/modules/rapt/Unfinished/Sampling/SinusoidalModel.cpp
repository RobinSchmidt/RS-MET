



/*

Ideas:

-frequency estimation:
-to find the actual frequency of a partial in an FFT frame, often parabolic interpolation of the
 log magnitude is used (the frequency is to be taken the maximum of the parabola that goes through
 3 adjacent bins)
-idea for refining this estimate: the transform a single sinusoid is given by an appropriately 
 shifted (and scaled) transform of the window function - try to find the optimum cross-correlation 
 lag, such that maximizes the cross-correlation between the frame's spectrum and the 
 window-transform...restricted to a certain correlation length (maybe 3...20)...maybe it could
 be useful to do it in the complex domain? maybe that would automatically account for phase
 information?
 -or: use a gaussian window - in this case the parabolic interpolation is actually exact:
  https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html

-peak picking/continuation:
-to match a peak in the k-th frame with exisitng active sinusoids, find the one that is closest in 
 frequency to a linear continuation of (k-1)-th and (k-2)-th frame 
-maybe special treatment has to be given to sine tracks that currently have only 1 frame - or maybe 
 rule them out by requiring that every stable partial must be at least 2 frames long to be 
 considered a partial at all...but that would require the peak picking algorithm not at one frame
 at a time but (at least) at 2 - but that may be a good idea anyway
-test the algorithm by using two sines that cross each other in the spectrogram (one ascending, one
 descending)

-create different variations of analysis/synthesis algorithms

*/