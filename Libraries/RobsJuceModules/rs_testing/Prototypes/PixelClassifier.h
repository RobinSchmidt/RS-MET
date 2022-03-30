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


  /** Classifies a pixel as belonging to class c when it has at least one neighbor that satisfies 
  the given predicate. That means, if for the pixel at img(i,j) the predicate holds for the pixel
  and any of its neighbors, the corresponding C(i,j) element is set to c, otherwise C(i,j) is left
  as is where C refers to the pixelClasses matrix that was passed to the constructor. The predicate
  always takes two pixels as input: the pixel currently under investigation, i.e. img(i,j), and one
  of its neighbors, e.g. img(i,j+1) or img(i-1,j), etc. The 2nd argument to pred is supposed to be 
  any of img(i,j)'s 8 neighbors, the first is img(i,j) itself. Edge pixels have only 5 neighbors, 
  corner pixels have 3. The function takes this into account. */
  template<class P> 
  void classifyWhenHasNeighborWith(char c, P pred);
  // API: maybe swap parameters, use "If" instead of "When" (shorter), maybe use assign instead of
  // classify (also shorter)

  /** Classifies a pixel as belonging to class c when it has no neighbor that satisfies the given 
  predicate. Note that this is a different condition than having a neighbor that satisfies the
  negated predicate. This can be confusing. todo: verify...explain...give examples */
  template<class P> 
  void classifyWhenHasNoNeighborWith(char c, P pred);



protected:

  /** Checks, if the given predicate is true for any of the neighboring pixels of pixel x,y. The 
  predicate should take the center pixel's value as first argument and the neighbor pixel's value
  as second argument and return true, iff the predicate holds for this pair of pixels. This 
  function is meant to be used only for interior pixels, hence the suffix _I. */
  template<class P> bool hasNeighborWith_I(int x, int y, P pred); // P: predicate

  template<class P> bool hasNeighborWith_T(int x, int y, P pred); // variant for the top edge
  template<class P> bool hasNeighborWith_B(int x, int y, P pred); // variant for the bottom edge
  template<class P> bool hasNeighborWith_L(int x, int y, P pred); // variant for left edge
  template<class P> bool hasNeighborWith_R(int x, int y, P pred); // variant for right edge

  template<class P> bool hasNeighborWith_TL(P pred); // variant for top-left corner
  template<class P> bool hasNeighborWith_TR(P pred); // variant for top-right corner
  template<class P> bool hasNeighborWith_BL(P pred); // variant for bottom-left corner
  template<class P> bool hasNeighborWith_BR(P pred); // variant for bottom-right corner


  const rsImage<TPix>& img;
  /**< Reference to the image, whose pixels we want to classify. Must be assigned at 
  construction. */

  rsImage<char>& classes;
  /**< Reference to an "image" which is actually just a 2D array of class labels for the pixels in
  our img member. */

};


template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_I(int i, int j, P pred)
{
  rsAssert(isInInterior(i, j), "Function is made only for interior pixels."); 
  // Trying to use it for boundary pixels will lead to an access violation. For boundary pixels,
  // separate implementations exist and higher level code is supposed to use these.

  float p = img(i, j);  // pixel value

  // Check against direct neighbors (maybe factor out):
  if( pred(p, img(i-1,j)) ) return true;
  if( pred(p, img(i+1,j)) ) return true;
  if( pred(p, img(i,j-1)) ) return true;
  if( pred(p, img(i,j+1)) ) return true;

  // Check against diagonal neighbors (maybe factor out):
  if( pred(p, img(i-1,j-1)) ) return true;
  if( pred(p, img(i-1,j+1)) ) return true;
  if( pred(p, img(i+1,j-1)) ) return true;
  if( pred(p, img(i+1,j+1)) ) return true;

  // Predicate holds for all neighbor pixels in 3x3 neighborhood:
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_T(int i, int j, P pred)
{
  // Code is the same as in hasNeighborWith_I but all lines involving j-1 have been deleted 
  // because in the top row, j-1 is not a valid y-coordinate
  rsAssert(isAtTopEdge(i, j, img), "Made for top-edge pixels");
  rsAssert(!isAtCorner(i, j, img), "Not made for corner pixels");
  float p = img(i, j); 
  if( pred(p,img(i-1,j  )) ) return true;
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i-1,j+1)) ) return true;
  if( pred(p,img(i+1,j+1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_B(int i, int j, P pred)
{
  // Lines with j+1 have been removed:
  rsAssert(isAtBottomEdge(i, j, img), "Made for bottom-edge pixels");
  rsAssert(!isAtCorner(   i, j, img), "Not made for corner pixels");
  float p = img(i, j);
  if( pred(p,img(i-1,j  )) ) return true;
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j-1)) ) return true;
  if( pred(p,img(i-1,j-1)) ) return true;
  if( pred(p,img(i+1,j-1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_L(int i, int j, P pred)
{
  // Lines with i-1 have been removed:
  rsAssert(isAtLeftEdge(i, j, img), "Made for left-edge pixels");
  rsAssert(!isAtCorner( i, j, img), "Not made for corner pixels");
  float p = img(i, j);
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j-1)) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i+1,j-1)) ) return true;
  if( pred(p,img(i+1,j+1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_R(int i, int j, P pred)
{
  // Lines with i+1 have been removed:
  rsAssert(isAtRightEdge(i, j, img), "Made for right-edge pixels");
  rsAssert(!isAtCorner(  i, j, img), "Not made for corner pixels");
  float p = img(i, j);
  if( pred(p, img(i-1,j  )) ) return true;
  if( pred(p, img(i,  j-1)) ) return true;
  if( pred(p, img(i,  j+1)) ) return true;
  if( pred(p, img(i-1,j-1)) ) return true;
  if( pred(p, img(i-1,j+1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_TL(P pred)
{
  // Lines with i-1 and j-1 have been removed:
  int i = 0;
  int j = 0;
  float p = img(i, j);
  if( pred(p,img(i+1,j  )) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i+1,j+1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_TR(P pred)
{
  // Lines with i+1 and j-1 have been removed:
  int i = 0;
  int j = img.getWidth()-1;  // maybe use img.getRight()
  float p = img(i, j);
  if( pred(p,img(i-1,j  )) ) return true;
  if( pred(p,img(i,  j+1)) ) return true;
  if( pred(p,img(i-1,j+1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_BL(P pred)
{
  // Lines with i+1 and j-1 have been removed:
  int i = img.getHeight()-1;  // maybe use img.getBottom
  int j = 0;
  float p = img(i, j);
  if( pred(p, img(i-1,j  )) ) return true;
  if( pred(p, img(i,  j+1)) ) return true;
  if( pred(p, img(i-1,j+1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
bool rsPixelClassifier<TPix>::hasNeighborWith_BR(P pred)
{
  // Lines with i+1 and j+1 have been removed:
  int i = img.getHeight()-1; 
  int j = img.getWidth()-1; 
  float p = img(i, j);
  if( pred(p, img(i-1,j  )) ) return true;
  if( pred(p, img(i,  j-1)) ) return true;
  if( pred(p, img(i-1,j-1)) ) return true;
  return false;
}

template<class TPix>
template<class P> 
void rsPixelClassifier<TPix>::classifyWhenHasNeighborWith(char c, P pred)
{
  int w = img.getWidth();
  int h = img.getHeight();

  // Classify interior pixels:
  rsImage<char>& C = classes;  // shorthand
  for(int j = 1; j < h-1; j++) 
    for(int i = 1; i < w-1; i++) 
      if(hasNeighborWith_I(i, j, img, pred)) 
        C(i, j) = c;

  // Classify edge pixels (excluding corners):
  for(int i = 1; i < w-1; i++) { 
    if(hasNeighborWith_T(   i,   0,   img, pred)) C(i, 0  ) = c;    // top row
    if(hasNeighborWith_B(   i,   h-1, img, pred)) C(i, h-1) = c; }  // bottom row
  for(int j = 1; j < h-1; j++) {    
    if(hasNeighborWith_L(   0,   j,   img, pred)) C(0,   j) = c;    // left column
    if(hasNeighborWith_R(   w-1, j,   img, pred)) C(w-1, j) = c; }  // right column

  // Classify corner pixels:
  if(hasNeighborWith_TL(img, pred)) C(0,   0  ) = c;
  if(hasNeighborWith_TR(img, pred)) C(w-1, 0  ) = c;
  if(hasNeighborWith_BL(img, pred)) C(0,   h-1) = c;
  if(hasNeighborWith_BR(img, pred)) C(w-1, h-1) = c;
}

template<class P> 
void classifyWhenHasNoNeighborWith(char c, P pred)
{
  int w = img.getWidth();
  int h = img.getHeight();

  // Interior:
  rsImage<char>& C = classes;  // shorthand
  for(int j = 1; j < h-1; j++) 
    for(int i = 1; i < w-1; i++) 
      if(!hasNeighborWith_I(i, j, img, pred)) 
        C(i, j) = c;

  // Edges:
  for(int i = 1; i < w-1; i++) { 
    if(!hasNeighborWith_T(   i,   0,   img, pred)) C(i, 0  ) = c;    // top row
    if(!hasNeighborWith_B(   i,   h-1, img, pred)) C(i, h-1) = c; }  // bottom row
  for(int j = 1; j < h-1; j++) {    
    if(!hasNeighborWith_L(   0,   j,   img, pred)) C(0,   j) = c;    // left column
    if(!hasNeighborWith_R(   w-1, j,   img, pred)) C(w-1, j) = c; }  // right column

  // Corners:
  if(!hasNeighborWith_TL(img, pred)) C(0,   0  ) = c;
  if(!hasNeighborWith_TR(img, pred)) C(w-1, 0  ) = c;
  if(!hasNeighborWith_BL(img, pred)) C(0,   h-1) = c;
  if(!hasNeighborWith_BR(img, pred)) C(w-1, h-1) = c;
}







#endif