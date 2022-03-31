template<class TPix>
bool rsPixelClassifier<TPix>::isInInterior(int i, int j)
{
  return i > 0 && j > 0 && i < img.getWidth()-1 && j < img.getHeight()-1;
}

template<class TPix>
bool rsPixelClassifier<TPix>::isAtLeftEdge(int i, int j)
{
  return i == 0;
}

template<class TPix>
bool rsPixelClassifier<TPix>::isAtRightEdge(int i, int j)
{
  return i == img.getWidth()-1;
}

template<class TPix>
bool rsPixelClassifier<TPix>::isAtTopEdge(int i, int j)
{
  return j == 0;
}

template<class TPix>
bool rsPixelClassifier<TPix>::isAtBottomEdge(int i, int j)
{
  return j == img.getHeight()-1;
}

template<class TPix>
bool rsPixelClassifier<TPix>::isAtCorner(int i, int j)
{
  return i == 0                && j == 0                   // top-left
    ||   i == 0                && j == img.getHeight()-1   // bottom-left
    ||   i == img.getWidth()-1 && j == 0                   // top-right
    ||   i == img.getWidth()-1 && j == img.getHeight()-1;  // bottom-right
}




//=================================================================================================

/*

ToDo:
-Maybe cache w,h, for convenience so we don't have to write img.getWidth() etc. all the time 
 And/or use r,b fro right,bottom to also get rid of the -1 all the time
-Maybe move functions into.h file. The functions taking templatized predicates must go into the 
 .h file anyway so we may just as well drop everything into it. Or: use std::function for the 
 predicates and move all code into the cpp file.
-Use int for the classes. It's easily conceivable to need more than 256 classes.


*/
