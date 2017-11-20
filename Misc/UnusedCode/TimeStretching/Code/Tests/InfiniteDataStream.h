#ifndef RS_INFINITEDATASTREAM_H
#define RS_INFINITEDATASTREAM_H

namespace RSLib
{

  /**

  This class represents an array of data that conceptually extends from -infinity to +infinity 
  whereas the range of actually available data is of course confined to the interval 0...length-1. 
  Any accesses to data outside this range will be automatically return the zero value of the 
  respective class.

  */

  template<class T>
  class rsInfiniteDataStream
  {

  public:

        
    /** Constructor. Initializes the underlying data-pointer to NULL and the length to zero. */
    rsInfiniteDataStream()
    {
      data       = NULL;
      dataLength = 0;
    }


    /** \name Setup */

    /** Sets the pointer to the full input stream which is supposed to contain valid data of the 
    given length. */
    void setInputDataAddress(T *newAddress, int newLength)
    {
      data       = newAddress;
      dataLength = newLength;
    }


    /** \name Inquiry */

    /** Returns the length of the underlying data. */
    int getDataLength() const
    { 
      return dataLength; 
    }


    /** \name Data Access */

    /** Returns the value at the given index or zero if the index is out of range. */
    T getValue(int index) const
    {
      if( index < 0 || index > dataLength-1 )
        return T(0);
      return data[index];
    }

    /** Produces a buffer of output data of the given length starting at the given start index. If
    startIndex < 0 and/or startIndex+numValues is beyond the end of the underlying input data, the 
    returned buffer will be zero padded appropriately. */
    void getBuffer(T *buffer, int startIndex, int numValues) const
    {
      int ir = startIndex;
      int iw = 0;
      if( dataLength == 0 )
      {
        for(iw = 0; iw < numValues; iw++)
          buffer[iw] = T(0);
        return;
      }
      while(ir < 0 && iw < numValues)
      {
        buffer[iw] = T(0);
        ir++;
        iw++;
      }
      while(ir < dataLength && iw < numValues)
      {
        buffer[iw] = data[ir];
        ir++;
        iw++;
      }
      while(iw < numValues)
      {
        buffer[iw] = T(0);
        iw++;
      }
      // \todo: optimize this - maybe precompute the required pre/post zero-padding lengths and use 
      // memset and memcpy - but be aware that via memcpy only "flat" copies will be passed into 
      // the oputput buffer - maybe make a separate function getBufferFlatCopy for that to make it 
      // clear to client code
    }

    // \todo setValue, setBufferValues

  protected:

    int dataLength;
    T *data;

  };

}

#endif
