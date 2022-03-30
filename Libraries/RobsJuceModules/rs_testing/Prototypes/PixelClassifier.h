#ifndef RS_PIXELCLASSIFIER_H
#define RS_PIXELCLASSIFIER_H


/** A class for classifying pixels in images. It takes a reference to the image whose pixels we
want to classify in the constructor. The constructor also takes a second image of type char which 
serves to store the pixel classes, so it must have the same shape as the actual image being 
analyzed. Once constructed, the user can call various classification functions which assign the 
elements classes image/matrix/2D-array. We use char for the class labels because I think that 256 
different pixel classes should be more than enough in practice...but that can be upgraded to 
something like int, if needed. 

Classifying pixels according to various criteria can be important as a subroutine in certain
image processing tasks where pixels belonging to different classes need to be processed in 
different ways.

...tbc...  */

template<class TPix>
class rsPixelClassifier
{

public:

  rsPixelClassifier(const rsImage<TPix>& image, rsImage<char>& pixelClasses) 
    : img(image), classes(pixelClasses)
  {
    rsAssert(classes.hasSameShapeAs(img), 
      "The array for the pixel classes must have the same shape as the image");
  }

  /** Returns true, iff the pixel at coordinates x,y is an interior pixel of the given image, i.e.
  not a boundary pixel. */
  bool isInInterior(  int x, int y);

  bool isAtLeftEdge(  int x, int y);
  bool isAtRightEdge( int x, int y);
  bool isAtTopEdge(   int x, int y);
  bool isAtBottomEdge(int x, int y);
  bool isAtCorner(    int x, int y);




protected:

  const rsImage<TPix>& img;
  /**< Reference to the image, whose pixels we want to classify. Must be assigned at 
  construction. */

  rsImage<char>& classes;
  /**< Reference to an "image" which is actually just a 2D array of class labels for the pixels in
  our img member. */

};





#endif