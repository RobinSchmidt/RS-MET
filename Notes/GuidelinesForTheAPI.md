API Guidelines (Draft)
======================

This document is still quite messy and unordered and the current API follows these guidelines that
are laid out here only partially. It is more like an aspirational document at the moment. I'm still
working on it.


Interface Considerations
-------------------------

- We should make the API consistent with regard to 

  - Naming conventions (including use of abbreviations and acronyms)

  - Parameter ordering (for functions with similar parameter lists)

  - Coefficient updating strategies (see below). Although, that's more an implementation rather than
    an interface consideration.

- Maybe we should use nested namespaces like RAPT::Data, RAPT::Math, RAPT::Filters, etc. 
  to get the functionality more ordered and also to render a better doxygen 
  documentation.

- We should anticipate using the C++20 module feature. Don't use it just yet (we want to 
  remain C++17 compatible for a while), but organize the structure in a way that 
  makes it easy to switch to using modules later  
  https://www.modernescpp.com/index.php/c-20-open-questions-to-modules  
  https://vector-of-bool.github.io/2019/01/27/modules-doa.html  
  Maybe the code can remain as is and the modules can be implemented on top of it rather similar
  like we now create the "JUCE module"? I'm not sure how modules interact with templates, though. 
  Maybe for a module definition, we need to commit to a specific type? In this case, we could
  perhaps create a reference module for RAPT using explicit instantiations with the types that make
  the most sense. That would also serve as documentation for the intended use.

- Anticipate using C++20 concepts at some point. It turns out that the following two concepts could
  make sense to define: Signal and Parameter. Observe the prevalence of the use of TSig and TPar
  throughout RAPT. These concepts are intertwined and should be defined in such a way that it is
  possible to perform arithmetic operations between signals and parameters, the result of which
  should be a signal. Rationale: Singals might be SIMD vectors and parameters the corresponding
  scalar type. Typically, both types represent continuous (i.e. real) numbers such as float, double,
  maybe also fixed point or multiprecision (I didn't try, though). Sometimes SIMD vectors make sense
  for signals _and_ parameters, sometimes only for signals.

- Maybe we should rename rosic to ramp (Rob's Audio and Music Processors)

- In rapt (but maybe not in rosic) consistently use omega = 2*pi*freq/sampleRate  in filters instead
  of frequency and sampleRate which is a redundant and inefficient parametrization (typically
  requires more computations to eventually get the coeffs). in rosic, a freq-and-sampleRate parameterization can be kept as convenience feature, but rapt is supposed to be more low-level
  and performance oriented.

- Maybe split the member variables of filter implementations into a "state" and "coeffs" part to
  allow for easier optimization of memory usage on higher levels. That may make a big difference
  especially when many filters with equal settings are used. Maybe even hae a 3rd class: Params for
  the user parameters. Example:
  ```
  class rsBiquad
  {
 
  public:
   
    struct Params
    {
      enum Type
      {
	      bypass = 0,
	      lowpass,
	      highpass,
	      // ...
	      numModes
	    };
   
      double sampleRate = 44100, frequency = 1000, quality = 1/sqrt(2);
	    Type type = bypass;
    };
   
    struct Coeffs
    {
      double b0, b1, b2, a1, a2;
    };
   
    struct StateDF1
    {
      double x1, x2, y1, y2;
    };
   
   
  protected:
   
    Coeffs   coeffs;
    StateDF1 state;
 
   };
   ```
 
   Maybe the classes for Params, Coeffs, State should be dragged out into classes in their own right 
   like rsBiquadParams, rsBiquadCoeffs, rsBiquadState. But that increases the surface area of the 
   library so it may be better to avoid it. Also: when templatizing, coeffs and params may need to 
   have a different template parameter than tse state (for example: Params/Coeffs use float and
   State uses rsFloat32x4), so we can't just propagate down temaplte params from rsBiquad to state -
   unless rsBiquad has two templated params TSig, TPar - which it probably should have anyway.
 
- Maybe preferably use `/**< ... */` instead of `/** ... */` for the doxygen stuff for member
  functions. It's just nicer when the function name appears as "headline" above the docstring in the
  code. But check, if intellisense still displays the docstrings when doing so. It also may not be
  so suitable for in-class function definitions (but maybe we should avoid them anyway unless they
  are one-liners we can still put them in the header file after the class declaration to allow (or
  even demand) inlining)
 

### Interface Consistency

- Template parameter names: 

  - Use `TSig` for the signal type. May later be replaced by a concept "Signal"

  - Use `TPar` for the parameter type. May later be replaced by a concept "Parameter" or "Param".

  - Some other names: `TPix` for pixels, `TCoef` for coefficients, `TVal` for values, `TArg` for
    function arguments, `Tx`, `Ty` for input and output types of functions  ...

- I use of camel case: My rationale for prefering `CamelCase` over `snake_case` is that the names
  tend to be shorter while the word separation via the capitalization is still obvious enough. The
  main argument for snake case is typically that the words are more clearly separated which is true
  but in my opionion, the separation in camel case is clear enough such that this marginal
  improvement does not justify the lengthening. Also, JUCE uses it as well although that had nothing
  to do with my decision. I made that long before even knowing JUCE. 
    
- I use `rs`-prefixing for class and function names. My rationale for prefixing everything with `rs`
  is to avoid name clashes. Yes, I know - that's what namespaces are there for but sometimes it is
  really inconvenient to always use `RAPT::someFunction()` instead of `rsSomeFunction()` and also I
  sometimes need my own versions of standard functions like `rsSin()` instead of `sin()` because
  the templatized nature of RAPT sometimes requires that I'm able to provide custom
  implementations of standard functions. Think, for example, a `sin()` function for SIMD vectors.
  If a class needs to compute a sine for some type `T`, I would just let it call `rsSin()` and
  provide a suitable explicit specialization of `rsSin()` for that type `T` and the templatized
  code of the class would the compile just fine when the class template is instantiated for type
  `T`. ...TBC...

- Class names 

  - I use `UpperCamelCase` with lower case prefix `rs` for classes. Example: `rsLadderFilter`. 
  
- Function names and signatures:

  - For free functions, also use camel case. Example `rsExp()`, `rsSin()`. It's actually supposed to
    be `lowerCamelCase` but the `rs`-prefix requires me to captitalize the function name anyway.

  - For member functions of classes I use bona fide `lowerCamelCase`. Example: `setCutoff()`.

  - ToDo: Have a `getSample()` and `processBlock()` member functions consistently in all DSP 
    classes.

  - If the module produces stereo samples, use processFrame because `getSample()` is only suitable
    for single channel (i.e. mono) processing. Maybe that makes it redundant but it is nice to be
    able to write things like:  
    `out = env.getSample() * filter.getSample(osc.getSample());`

  - ToDo: Consistently use pointers (not references) for output variables. It makes it visible at
    the call site, what is an output that may be (over)written.

  - The argument order should be consistent for functions that do similar things - especially in
    rsArrayTools (potentially dangerous (i.e. quietly breaking) change!)

  - The units (seconds, milliseconds) for parameters should be the same in all classes (dangerous
    change)

  - In rapt, we should probably not deal with physical units at all and instead use normalized units
    like samples or omega: $\omega = 2 \pi f / f_s$,  etc.

- Consistent use of enum class for choices

- (Maybe) avoid free functions. Wrap them into namespaces. Maybe use sub-namespaces
    (RAPT::Filters, RAPT::Generators, etc.). Or maybe (ab)use classes for collections of functions.
    The functions can then be static member functions. For some reason, I do not really like deeply
    nested namespaces. But maybe these are two separate things. Creating a sub-namespace just for a
    function collection may be "overkill" (I think of namespaces as big things) but for the bigger
    compartements of the library like "Filters", "Generators", etc., namespaces may actually be 
    appropriate.

  - I've been thinking about using std::vector for all i/o of the the non-realtime classes but some
    client code may use something else, so suing plain arrays is most flexible. Anyone can use it.
    So, the low-level number crunching code should be based on passing around C-style arrays.

  - Consistently make embedded objects accessible either via having them as public members or by
    providing getters. The former seems better - simpler client-side syntax and does the exact
    same thing. But maybe do this only in cases where the embedder does not need to maintain any 
    constraints on the embedded objects. Otherwise maybe let the embedder be a sort of facade (as in
    the facade design pattern).

  - Use consistently int - not size_t or something (in the sinusoidal model) ...or maybe size_t is
    better - but then we need a convention other than -1 to indicate things like "not found" in
    functions like findIndex - perhaps N would be the most obvious convention (N = length of array).
    That would also be consistent with STL conventions:  
    https://github.com/fish-shell/fish-shell/issues/3493


  - Use get/set consistently. Bad: Matrix::eigenvalues, Good: Matrix::getEigenvalues
    ...but only for non-static member functions - for static ones -> no get

  - Avoid heap allocations for temporary arrays - use workspace parameters instead but keep the
    functions that do heap allocations for convenience - but document all heap-allocations, i.e.
    write a warning. Annotate all functions that do heap allocation. Or maybe more generally:
    Annotate all functions that are unsafe to call in a relatime context. But maybe the default
    assumption should be that the function is unsafe and we should annotate those that are
    realtime-safe.

  - Use nouns like "getProduct" when the function returns an object (numbers qualify as well) but a
    verb like "multiply" when the function performs some action on passed inputs like multiplying
    array element-wise

  - Use abbreviations consistently in function names and their parameters. Use: Dist: Distortion, 
    Distro: Distribution, Freq: Frequency, Calc: Calculate, Coeff: Coefficient, Reso: Resonance, 
    Cyc: Cycle, Oct: Octave, Sec: Second,
 
  - See also: https://github.com/RobinSchmidt/RS-MET/wiki/Standards


### Use of Abbreviations and Acronyms

Generally we want names to be descriptive such code readers can correctly guess their meaning but we
also want them to be short in order to not blow up the verbosity of the code (which is bad for
 readibility). By their nature, descriptive names tend to be longer so we must strike a tradeoff
between descriptiveness and brevity. One way to do this is to introduce abbreviations and acronyms.
If these things are used in the API, they should be used consistently throughout the library. When
acronyms are used in the context of CamelCase

...TBC...


Implementation Considerations
-----------------------------

- Whether algo parameters are updated directly in setCutoff, setResonance, etc. or are updated in
  getSample based on an "upToDate" flag should be consistent. Or maybe that can be decided on a
  case-by-case basis? It currently is. But maybe it shouldn't. But enforced consistency may make it
  harder to optimize cases individually.

- What about thread safety? Should probably be completely abolished, at least in RAPT. 
  Synchronization should be dealt with on a higher level. Maybe in rosic, but probably jura.

- How polyphony can be handled: Each DSP class may have a simple, monophonic  implementation and
  some may optionally have a polyphonic version (maybe as subclass with suffix "Poly"). It would certainly be nice to have some sort of automatic way to turn monophonic DSP algorithms into 
  polyphonic ones. But that may lead to suboptimal implementation because typically, the settings
  and parts of the state can be shared among the voices. But which parts can and which can't be
  shared depends on the particluar algorithm
 

Performance Considerations
--------------------------

- Use const and constexpr whereever possible

- Always declare members in order of descencing size (reduces padding)

- Use workspaces for operations on arrays that need auxiliary memory to avoid heap allocations for
  temporary buffers. Maybe have convenience functions that allocate a workspace internally (e.g. 
  by declaring a local std::vector) and then calling the workspace based function with that.

- In low level DSP classes, try to avoid virtual functions as much as possible, especially in small
  objects (like 1st or 2nd order filters) where the vtable would significantly increase the size
  relatively to without vtable. So, as a general rule, use only compile time polymorphism, if 
  polymorphism is used at all. Exceptions must be justified on a case by case basis. See:  
  https://www.youtube.com/watch?v=wF6OfkqK6fE  Most C++ Virtual Calls Are NOT Dynamic Dispatch  
  https://www.youtube.com/watch?v=NWsdK_ly_fw  C++ Objects Are NOT What You Think (Inside vtables & ABI)


### Updating strategies

When a DSP algorithm has many parameters with corresponding setters, the question arises how to 
trigger the recalculation of the internal coefficients. Take, for example, a filter that has 
setters like setFrequency(double newFreq), setResonance(double newReso). The simplest strategy would 
be that both functions immediately trigger the coefficient recomputation. This would have the 
advantage that the coefficients would always be correct in the sense that they always reflect the 
desired user settings. The disadvantage is that in the event of calling setFrequency and 
setResonance in succession (a situation, that occurs often), they would be computed twice and only
the last computation is actually relevant - so it would be wasteful. With more than two setters, it
will get even worse. Another strategy is to have a single setup(double newFreq, double newReso)
function that takes both parameters at once and then triggers only one coefficient recomputation.
This solves the efficiency issue - but it may get unwieldy for classes with lots of setters. It will
also not play well with the callback system in jura::Parameter - such parameters expect to be wired
to a single argument setter. Yet another strategy is to go back to have single parameter setters but
do not let them trigger a recomputation right away. Instead, they just set a "dirty" flag. I assume 
this to be an atomic operation (verify!). In getSample(), this dirty-flag is checked and if it is 
true, the recomputation is triggered there. Advantages: It solves the efficiency issue and could 
additionally make the setters thread-safe because the recomputation is deferred to the audio-thread 
whenever the setters are called from a different thread (such as a GUI thread). The disadvantage is 
that after calling a setter, the coeffs do not reflect the user settings and we have one thing more 
to do in getSample. This overhead may not matter for complex algorithms - but for very simple ones, 
it may.


Safety Considerations
---------------------

- If a pointer cannot be a nullptr, say so by using not_null.

- Pay attention to edge cases. Sometimes they are handled automatically correct. Implement the
  general case first and test (and reason about) what that implementation will do in the edge case.


General Guidelines
------------------

A good api should be:

- Easy to learn an memorize.

- Hard to use incorrectly.

- Easy to extend.

- Complete with regard to features that client code may reasonably expect.

- Lead to readable client code. Library code readability could be sacrificed, if it leads to more
  readable client code. Client code readability is more important than library code readability.


Resources
---------

- Here are the guidelines recommended by the c++ commitee: 
  https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md

- "Modern c++" guidelines say, that we should avoid raw "new" and "delete" and use smart-pointers
   instead std::unique_ptr, std::shared_ptr

  - Try to do that, where reasonable (but always take potential performance degradation into
    account)

   - where not reasonable, convey information of ownership takeover in the respective function names
     instead of addChildFoo( new Foo ) use things like addOwnedChildFoo( new Foo )
 
- Some more advice:  
  https://www.linkedin.com/pulse/what-bad-code-how-write-clean-sultan-mahamud/  
  It recommends to avoid flag parameters. I use them in a couple of places. Maybe I should try to 
  change that.

- Seems like #pragma once is not really recommended as replacement for include guards - so we should
  stick to them, fo the time being  
  https://stackoverflow.com/questions/1143936/pragma-once-vs-include-guards

- This talk has interesting stuff about API design:  
  https://www.youtube.com/watch?v=5tg1ONG18H8  
  One thing: avoid adjacent parameters of the same type to minimize the chance of the caller getting
  the order wrong

- Some arguments about why free functions may be preferable over member functions:

  https://www.youtube.com/watch?v=WLDT1lDOsb4  CppCon 2017: Klaus Iglberger “Free Your Functions!”
  https://www.youtube.com/watch?v=nWJHhtmWYcY  Free your functions! - Klaus Iglberger - Meeting C++ 2017
  
  Examples for when free functions are better than member functions in the RAPT library could be 
  functions like rsAbs() for absolute value,  rsDot(), rsCross() for dot-and cross-product of 
  vectors, rsDet() for a determinant of a matrix etc. The rationale for prefering a free function 
  call like  d = rsDot(v, w)  over a member function call like  d = v.dot(w)  to compute the dot 
  product between two vectors v and w (for example, inside some geometry algorithm that should work 
  with generic(!) vectors) is that not all vector classes that we may want to use in our algo may 
  implement the dot() member function. For example, std::vector doesn't - and we also can't add one
  (unless we use a subclass of std::vector). However, a free function like rsDot() taking two 
  std::vectors could easily be added. Similar considerations apply to functions for computing a 
  norm, determinant, etc.

  He makes the argument, that convenience functions like clearAll() that internally just call lower 
  level clearing functions like clearThis(), clearThat() etc. should be pulled out of the class 
  because then we have one function less to worry about potentially messing up class invariances. In
  general, I think, this applies to all public member functions that just call other public member 
  functions. The disadvantage might be discoverability and generated documentation - it may not list
  the non-member function.



Misc
----





<br><br><br><br><br><br><br><br><br><br>
----------------------------------------------------------------------------------------------------
ToDo
----

- Re-organize the document. Maybe sections about: "Interface Consistency", "Implementation 
  Consistency", "Documentation Consistency", "Performance Considerations", "Safety Considerations", ....

 - Check use of std::list<> in rosic::PolymorphicInstrumentVoice

 - Define a consistent strategy for CameCasing words that can be seen as compound words or as two
   words like "sample rate" vs "samplerate" or "wave shape" vs "waveshape". I currently just 
   capitalize every part of a word that "could" be seen as a word in its own right. That is, I use
   things like WaveForm, WaveShape, SampleRate, OverSampling although some of them have become 
   single words already in common language use (like "oversampling" - nobody writes "over 
   sampling"). I'm not yet sure how to handle this best. Other examples: BandWidth, EigenValue,
   EigenVector, EigenSpace, EigenFunction ..what about LowPass, HighPass, BandPass, BandStop, etc.?
   It almost appears like compound words in English fall on a continuous spectrum between "should be
   considered one word" and "should be considered two words". Words like lowpass would be at the 
   "one word" end and words like "sample rate" more near the "two words" end. That makes it 
   difficult to opt for one single consistent rule. Up to now, I kinda was just winging it and 
   decided on a case by case basis but it seems like having a consistent rule in place may be a good
   idea. Maybe to captitalize consistently anything that could be regarded a word in its own right
   may be the way to go. It may be slightly "wrong" in some cases with regard to the English
   language but I think, it may be wrong in an acceptable way. It may remind the reader of the
   respective etymology of the word and thereby increase descriptiveness.


Ideas for abbreviations:

Amplitude:       Amp  ,
Argument:        Arg  ,
Analog:          Ana  ,
Attack:          Att  ,
Bandwidth:       Bw (in Hz: BwHz, in octaves: BwOct)  ,
Calculate:       Calc  ,
Coefficient:     Coeff  ,
Compare:         Comp  ,
Compartment:     Comp  ,
Complex:         Comp  ,
Component:       Comp  ,
Compression:     Comp  ,
Compute:         Comp  ,
Context:         Ctx  ,
Cycle:           Cyc  ,
Damping:         Damp  ,
Decay:           Dec  ,
Decibel:         Db  ,
Digital:         Digi  ,
Distance:        Dist  , 
Distortion:      Dist, Distort  ,
Distribution:    Dist, Distri, Distrib, Distro  ,
Eigenvalue:      EigVal  ,
Eigenvector:     EigVec  ,
Envelope:        Env  ,
Frequency:       Freq  ,
Filter:          Flt, Filt  ,
Function:        Fun, Func  ,
Hertz:           Hz  ,
Histogram:       Hist, Histo  ,
Index:           Idx  ,
Instance:        Inst  ,
Instantaneous:   Inst, Insta, Instant  ,
Instrument:      Inst, Instrum  ,
Matrix:          Mat  ,
Millisecond:     Ms  ,
Modifier:        Mod  ,
Modulation:      Mod  ,
Module:          Mod  ,
Modulus:         Mod  ,
Octave:          Oct  ,
Oscillator:      Osc  ,
Parameter:       Par, Param  ,
Phase:           Phs  ,
Release:         Rel  ,
Resolution:      Res, Reso, Resol  ,
Resonance:       Res, Reso, Reson  ,
Second:          Sec  ,
Singal:          Sig  ,
Spectrum:        Spec  ,
Specification:   Spec  ,
Sustain:         Sus  ,
Tensor:          Tens, Tns  ,
Value:           Val  ,
Vector:          Vec  ,

Aim for 1 or 2 syllables. Try to make the abbreviations unique - not like with "Comp". When writing
code and trying to find a suitable abbreviation, we should really first look up this dictionary to
see, if there already is an abbreviation used for that somewhere else in the library. We should
maintain a dictionary somewhere for this purpose.

Acronyms:

I think, we may capitalize only the first letter in an acronym, e.g. write Fft rather than FFT.
Examples:

Amplitude Modulation:               Am  ,
Attack Decay Sustain Relase:        Adsr  ,
Beats per Minute:                   Bpm  ,
Digital Signal Processing:          Dsp  ,
Fast Fourier Transform:             Fft  ,
Finite Impulse Response:            Fir  ,
Frequency Modulation:               Fm, FreqMod  ,
Graphical User Interface:           Gui  ,
Infinite Impulse Response:          Iir  ,
Low Frequency Oscillator:           Lfo, LowFreqOsc  ,
Lowpass filter:                     Lpf  ,
Linear Predictive Coding:           Lpc  ,
Multi Segment Envelope Generator:   Mseg  ,
Phase Modulation:                   Pm, PhaseMod  ,
Pulse Width Modulation:             Pwm  , 
Ring Modulation:                    Rm, RingMod  ,
Reverb Time 60:                     Rt60  ,
Semitone:                           St  ,


Symbols:

Normalized Radian Frequency:      omega: $\omega$  ,
Time Constant:                    tau: $\tau$  ,