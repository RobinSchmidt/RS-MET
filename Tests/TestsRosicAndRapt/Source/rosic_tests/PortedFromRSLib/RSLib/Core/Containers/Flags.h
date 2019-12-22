#ifndef RS_FLAGS_H
#define RS_FLAGS_H

namespace RSLib
{

  /**

  This is a class to maintain a set of 8 boolean flags that can be individually accessed (to set 
  them and to retrieve them). The class is intended to save memory when a whole bunch of flags is 
  needed. Functionality wise, it's similar to the STL bitset.

  \todo: write classes rslFlags16, rslFlags32, rslFlags64 in an entirely analogous way - maybe 
  this can be done with a template class? ...with a mask that is T::maxValue/2
  or make a class rsFlagsN that can hold an arbitrary number of flags using an array of
  rsFlags8 - or do both

  \todo write a class rsFlagArray that contains an array of rsFlags8 and provides space 
  8*numArrayElements flags, access would work like that: 
  rsFlagArray::isFlagTrue(int index) 
  { 
    return flagArray[index/8].isFlagTrue[index%8]; 
  }
  where the flagArray member is of type rsArrayTools<rsFlags8>. to optimize, one could perhaps use
  flagArray[index>>3].isFlagTrue[index - index>>3]; -> verify this
  maybe call this class rsBitArray


  */

  class RSLib_API rsFlags8
  {

  public:


    /** \name Construction/Destruction */

    /** Standard constructor. Initializes all flags to zero. */
    rsFlags8() : flags(0) 
    { 
    
    }


    /** \name Operators */

    bool operator==(const rsFlags8& other) const
    {
      return flags == other.flags;
    }

    bool operator!=(const rsFlags8& other) const
    {
      return flags != other.flags;
    }


    /** \name Setup */

    /** Sets the flag with the given index to true. The index should be in the rnage 0...7, if 
    it's greater than that, the function will do nothing. */
    inline void setFlagTrue(int index) 
    { 
      flags |= ((unsigned char) 128 >> index); 
    }

    /** Sets the flag with the given index to false. */
    inline void setFlagFalse(int index) 
    { 
      flags &= ~( (unsigned char) 128 >> index ); 
    }

    /** Sets all flags to true. */
    inline void setAllTrue() 
    { 
      flags = 255;
    }

    /** Sets all flags to false. */
    inline void setAllFalse() 
    { 
      flags = 0;
    }

    /** Sets the flag with the given index to the passed new boolean value.  */
    void setFlag(int index, bool newValue);

    /** Toggles the flag into to the other state - a true flag becomes flase and vice versa. */
    void toggleFlag(int index);


    /** \name Inquiry */

    /** Returns whether the flag with the given index is true. The index should be in the range 
    0...7, if it's greater than that, the function will return false. */
    inline bool isFlagTrue(int index) const 
    { 
      return ((flags << index) & (unsigned char) 128) > 0; 
    }

    /** Returns whether the flag with the given index is false. The index should be in the range 
    0...7, if it's greater than that, the function will return false. */
    inline bool isFlagFalse(int index) const 
    { 
      return !isFlagTrue(index); 
    }

  protected:

        
    /** \name Data */

    rsUint8 flags;

  };




  // rename to rsBitSet, use setTrue/setFalse/switch, isTrue




  class RSLib_API rsFlagArray
  {

  public:


    /** \name Construction/Destruction */

    /** Standard constructor. Initializes all flags to zero. */
    rsFlagArray(rsUint64 numDesiredFlags);

    /** Destructor */
    ~rsFlagArray();


    /** \name Setup */

    /** Sets the flag with the given index to true. The index should be in the rnage 0...7, if 
    it's greater than that, the function will do nothing. */
    inline void setFlagTrue(rsUint64 index) 
    { 
      rsSetBitTrue(flagGroups64[index/64], index%64);
    }

    /** Sets the flag with the given index to false. */
    inline void setFlagFalse(rsUint64 index) 
    { 
      rsSetBitFalse(flagGroups64[index/64], index%64);
    }

    /** Sets all flags to true. */
    void setAllTrue() 
    { 
      memset(flagGroups64, -1, (numFullGroups+1)*sizeof(rsUint64));
    }

    /** Sets all flags to false. */
    void setAllFalse() 
    { 
      memset(flagGroups64, 0, (numFullGroups+1)*sizeof(rsUint64));
    }

    /** Sets the flag with the given index to the passed new boolean value.  */
    //void setFlag(rsUint64 index, bool newValue);

    /** Toggles the flag into to the other state - a true flag becomes flase and vice versa. */
    //void toggleFlag(rsUint64 index);


    /** \name Inquiry */




    /** Returns whether the flag with the given index is true. */
    inline bool isFlagTrue(rsUint64 index) const
    { 
      return rsIsBitTrue(flagGroups64[index/64], index%64);
    }

    /** Returns whether the flag with the given index is false. */
    inline bool isFlagFalse(rsUint64 index) const 
    { 
      return !isFlagTrue(index); 
    }


    /** Returns true, iff all the flags are true. */
    bool areAllFlagsTrue();

    /** Returns true, iff all the flags are false. */
    bool areAllFlagsFalse();


    rsUint64 getNumFlags() const
    {
      return numFlags;
    }

    /** Returns the number of flags which are true. */
    rsUint64 getNumTrueFlags() const;

    /** Returns the next index at which a true flag occurs, including the passed index in the 
    search. If no true flags are found beyond and including the passed index, the returned index
    will be equal the number of flags (which is one unit out of range). */
    rsUint64 getNextTrueFlag(rsUint64 index)  const
    {
      while( (index < numFlags) && (!isFlagTrue( index )) )   // get rid of (int)
        index++;
      return index;
    }



  protected:


    rsUint64 getNumTrueTailFlags() const;
        
    /** \name Data */

    rsUint64 numFlags;
    rsUint32 numFullGroups;
    rsUint32 tailLength;
    rsUint64 *flagGroups64;

  };

}

#endif
