
template<class T>
void rsImageProcessor<T>::normalize(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  int N = img.getNumPixels();
  T min = rsArrayTools::minValue(p, N);
  for(int i = 0; i < N; i++)
    p[i] -= min;
  T max = rsArrayTools::maxValue(p, N);
  for(int i = 0; i < N; i++)
    p[i] /= max;
}
// this *may* be better numerically (less prone to roundoff errors) than the "fast" version - 
// needs test

template<class T>
void rsImageProcessor<T>::normalizeFast(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  int N = img.getNumPixels();
  T min = rsArrayTools::minValue(p, N);
  T max = rsArrayTools::maxValue(p, N);
  T scl = 1.f / (max-min);
  for(int i = 0; i < N; i++)
    p[i] = scl * (p[i] - min);
}

template<class T>
void rsImageProcessor<T>::normalizeJointly(rsImage<T>& img1, rsImage<T>& img2)
{
  //rsAssert(img2.hasSameShapeAs(img1));  // activate later
  using AT = rsArrayTools;
  int N = img1.getNumPixels();
  T* p1 = img1.getPixelPointer(0, 0);
  T* p2 = img2.getPixelPointer(0, 0);
  T min = rsMin(AT::minValue(p1, N), AT::minValue(p2, N));
  for(int i = 0; i < N; i++) {
    p1[i] -= min;
    p2[i] -= min; }
  T max = rsMax(AT::maxValue(p1, N), AT::maxValue(p2, N));
  for(int i = 0; i < N; i++) {
    p1[i] /= max;
    p2[i] /= max; }
}

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