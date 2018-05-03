#ifndef RS_INFINITEDATASTREAM_H
#define RS_INFINITEDATASTREAM_H

namespace RSLib
{

  /**

  This class represents an infinite array of data that conceptually extends from -infinity to 
  +infinity whereas the range of actually available data is of course confined to the interval
  0...length-1. Any attempt to read access of data outside this range will be automatically return 
  the zero value of the respective class and write access out of range will have no effect.

  \todo maybe rename to rsEndlessArray

  */

  template<class T>
  class rsInfiniteDataStream
  {

  public:

        
    /** \name Construction/Destruction */

    /** Constructor. The user should pass the pointer to the actual data values and the length of 
    the array. */
    rsInfiniteDataStream(T *dataPointer, int length)
    {
      setInputDataAddress(dataPointer, length);
    }


    /** \name Setup */

    /** Sets the pointer to the full input stream which is supposed to contain valid data of the 
    given length. */
    void setInputDataAddress(T *newAddress, int newLength)
    {
      data       = newAddress;
      dataLength = newLength;
      dummy      = T(0); 
    }


    /** \name Data Access */

    /** Returns the value at the given index or zero if the index is out of range. */
    T getValue(int index)
    {
      if( index < 0 || index > dataLength-1 )
        return T(0);
      return data[index];
    }

    /** Produces a buffer of output data of the given length starting at the given start index. If
    startIndex < 0 and/or startIndex+numValues is beyond the end of the underlying input data, the 
    returned buffer will be zero padded appropriately. */
    void getBuffer(T *buffer, int startIndex, int numValues)
    {
      int ir = startIndex;
      int iw = 0;
      while( ir < 0 && iw < numValues )
      {
        buffer[iw] = T(0);
        ir++;
        iw++;
      }
      while( ir < dataLength && iw < numValues )
      {
        buffer[iw] = data[ir];
        ir++;
        iw++;
      }
      while( iw < numValues )
      {
        buffer[iw] = T(0);
        iw++;
      }
      // \todo: optimize this - maybe precompute the required pre/post zero-padding lengths and use 
      // memset and memcpy - but be aware that via memcpy only "flat" copies will be passed into 
      // the oputput buffer - maybe make a separate function getBufferFlatCopy for that to make it 
      // clear to client code
    }

    /** Accesses the element at given index for reading and writing. For indices out of bounds, 
    zero will be returned for read access and write access will have no effect. */
    T& operator[](const int index)
    {
      if( index < 0 || index >= dataLength )
      {
        dummy = T(0); // required, because it may have been written to
        return dummy;
      }
      else
        return data[index];
    }

    // \todo setValue, setBufferValues

  protected:

    int dataLength;
    T *data;
    T dummy; // for handling out-of-bounds write access

  };

}

#endif
