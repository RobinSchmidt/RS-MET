



/*

Ideas:
-to find the actual frequencyof a partial in an FFT frame, often parabolic interpolation of the
 log magnitude is used (the frequency is to be taken the maximum of the parabola that goes through
 3 adjacent bins)
-idea for refining this estimate: the transform a single sinusoid is given by an appropriately 
 shifted (and scaled) transform of the window function - try to find the optimum cross-correlation 
 lag, such that maximizes the cross-correlation between the frame's spectrum and the 
 window-transform...restricted to a certain correlation length (maybe 3...20)...maybe it could
 be useful to do it in the complex domain? maybe that would automatically account for phase
 information?


*/