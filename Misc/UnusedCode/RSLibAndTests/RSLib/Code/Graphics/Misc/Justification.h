#ifndef RS_JUSTIFICATION_H
#define RS_JUSTIFICATION_H

namespace RSLib
{

  /**

  This class is used to represent different justifactions for graphics objects like top-left for
  standard text, etc. The default configuration when you create an rsJustification object is
  centered in both dimensions (horizontally and vertically). After creating the object, you can
  set the desired vertical and horizontal justifications with the "set" methods and inquire them
  with the "is" methods.

  \todo test this class

  */

  class RSLib_API rsJustification
  {

  public:

    /** \name Setup */

    void setJustifiedLeft();
    void setHorizontallyCentered();
    void setJustifiedRight();
    void setJustifiedTop();
    void setVerticallyCentered();
    void setJustifiedBottom();
    void setCentered();

    /** \name Inquiry */

    bool isJustifiedLeft() const;
    bool isHorizontallyCentered() const;
    bool isJustifiedRight() const;
    bool isJustifiedTop() const;
    bool isVerticallyCentered() const;
    bool isJustifiedBottom() const;
    bool isCentered() const;

    // todo provide factory-methods for commonly used justifictions like createTopLeft(), etc.


  private:

    /** \name Data */

    enum justificationFlagIndices
    {
      justifiedLeft = 0,
      justifiedRight,
      justifiedTop,
      justifiedBottom
    };

    rsFlags8 flags;

  };

}

#endif
