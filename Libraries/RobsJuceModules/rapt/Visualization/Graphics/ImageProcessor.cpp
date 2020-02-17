
template<class T>
rsImage<T> rsImageProcessor<T>::scaleUp(const rsImage<T>& img, int scl)
{
  int w = img.getWidth();
  int h = img.getHeight();
  rsImageF result(scl*w, scl*h);
  for(int x = 0; x < w; x++)  {
    for(int y = 0; y < h; y++) {
      for(int i = 0; i < scl; i++) {
        for(int j = 0; j < scl; j++) {
          result(scl*x+i, scl*y+j) = img(x, y); }}}}
  return result;
}
// todo:
// -allow different scaling factors for x and y
// -let the outer loop run over y and the inner over x