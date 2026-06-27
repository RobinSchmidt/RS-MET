
Algorithms Overview
===================


Math
----

### Polynomials

A polynomial is a function of the form:

$$ f(x) = \sum_{n=0}^N a_n x^n $$

...TBC...

### Rational Functions


### Matrices and Linear Algebra

The class `rsMatrix` is the basis for everything that has to do with matrices and linear algebra. It
is templatized on the data type of the matrix entries, so you can instantiate it with simple built 
in types like `float`, `double`, `std::complex<double>` but also with more complicated types like
`rsRationalFunction<double>` to represent matrices of transfer functions. Such matrices occurr, for
example, in the context of MIMO filters and feedback delay networks (FDNs). The class has a
baseclass `rsMatrixView` which already provides a lot of the needed functionality without owning the
actual data such that you can use it to "wrap" an existing 1D array of numbers and treat it as a
matrix (for example, by passing it to linear algebra algorithms) while managing the memory yourself.
The classes `rsMatrix` and `rsMatrixView` provide the low level matrix operations like matrix
addition, multiplication, elementary row-manipulations, etc. 

To do actual linear algebra like solving linear systems of equations, inverting matrices, finding
eigenvalues and -vectors, etc. there is another class called `rsLinearAlgebra` which is a collection
of static functions that typically operate on `rsMatrixView` objects, ...TBC...




### Linear Transforms


Filters
-------

### One Pole Filter

### Biquad Chain

### State Variable Filter

### Ladder Filter

### Engineers Filter

### Quantile Filter


Generators
----------

### Table Lookup Oscillator

### Pitch Dithering

### Fractal Pattern Synthesis


Modulators
----------

### Simple Envelope Generators (AD, ADSR, ...)

<!---
### Low Frequency Oscillators (LFO)
-->

### Multi Segment Envelope Generator (MSEG)

<!---
### Random Modulation
-->





Audio Analysis
--------------



Modeling and Resynthesis
------------------------





Under Construction
==================

### ResoReplace Filter


Ideas
=====




<br><br><br><br>
----------------------------------------------------------------------------------------------------

#### ToDo

- Maybe change the order of the sections to put the more high level stuff at the top. This is the
  stuff that is most likely to be used by client code. The deeper depths of the library such as
  all the basic math stuff (matrices, polynomials, etc.) is perhaps less relevant to users in the 
  sense that they are less likely to _directly_ use it. These are more things that I use internally.