template<class TPix, class TWgt, class TCor>
ImagePainter<TPix, TWgt, TCor>::ImagePainter(Image<TPix> *imageToPaintOn, Image<TWgt> *brushToUse)
{
  image = imageToPaintOn;
  brush = brushToUse;
}

