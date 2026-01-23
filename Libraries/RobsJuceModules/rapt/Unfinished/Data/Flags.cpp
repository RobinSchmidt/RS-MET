//using namespace RSLib;


//=================================================================================================
/*

ToDo:

- Implement a class rsMultiEnum than uses a single integer to represent the values of multiple
  enum variables in a similar way that the rsFlags8 class uses a single (8 bit) integer to 
  represent upt to 8 boolean variables. Use a template with template parameters being the 
  underlying integer type (e.g. rsUint32) and an variadic number of numeric template parameters to
  determine how many bits are being used for each enum. For example, an instantiaiton like
  rsMultiEnum<rsUint32, 5, 3, 2, 6, 1, 3, 8, 4> would represent the values from 8 enums where the
  first has up to 2^5 values, the second up to 2^3 values, the third up to 2^2 values and so on.
  To extract the values of the different enums, use bit-masking and -shifting. ...I'm not sure, if
  somehtin like that is possible, though. I think, it would be straightforward to do something like
  rsDoubleEnum, rsTripleEnum, etc. where we know in the implementation, how many enums the variable
  represents. But making it so flexible that one class can be used for all configurations - I'm not
  sure, how to do that. Maybe ask AI.


*/





