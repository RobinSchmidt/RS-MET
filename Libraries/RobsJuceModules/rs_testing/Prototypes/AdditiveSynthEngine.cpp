





/*

Ideas:

-Maybe split the partials conceptually into groups and let the user access certain processing 
 options per group. For example, it may make sense to have one group for even and one for odd 
 partials and let the user set their gains individually. It may also make sense to split into 
 low/mid/high frequencies or groups containing harmonic and inharmonic content.
-Maybe the sizes of the groups should be independent from the SIMD vector size, and maybe even
 variable: e.g. group 1 has 10 harmonics, group 2 has 17, group 3 has 42, etc. Maybe we can take
 inspiration from how groups are used in sfz.


Notes:

Currently, the largest SIMD vector size available for single precision float in standard 
hardware is 16 (for example in AVX-512). We use twice that value just in case that something
like AVX-1024 will become available in the future, in which case the code may directly benefit
from it without any change.
...hmm...not sure, if that's a good idea...i think, from a perceptual point of view, smaller 
groups may make more sense...but then, when we provide an API for setting up settings for
groups, their size may actually be different from the simd vector size

*/