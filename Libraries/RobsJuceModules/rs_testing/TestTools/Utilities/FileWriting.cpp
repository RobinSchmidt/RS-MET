#include "FileWriting.h" 

//#ifdef _MSC_VER
//#define _CRT_SECURE_NO_WARNINGS  // seems to have no effect - we still get the warning
//#endif

// see also: https://rosettacode.org/wiki/Bitmap/Write_a_PPM_file#C

bool writeImageToFilePPM(const char* path, unsigned char* buf, int w, int h)
{
  FILE* fd = fopen(path, "wb");  // "wb": write binary
  if(fd != NULL) {
    fprintf(fd, "P6\n%d %d\n255\n", w, h);
    fwrite(buf, 1, w*h*3, fd);
    fclose(fd); 
    return true; }
  else {
    rsError("Unable to open file");
    return false; }
}

bool writeImageToFilePPM(const rsImageF& img, const char* path)
{
  int w = img.getWidth();
  int h = img.getHeight();
  unsigned char* buf = new unsigned char[w*h*3];
  if(buf == nullptr) {
    rsError("Unable to allocate buffer"); 
    return false; }
  for(int y = 0; y < h; y++) {
    for(int x = 0; x < w; x++) {
      int i = y*w*3 + x*3;
      unsigned char gray = (unsigned char) (255 * img(x, y));
      buf[i+0] = gray;
      buf[i+1] = gray;
      buf[i+2] = gray; }}
  bool success = writeImageToFilePPM(path, buf, w, h);
  delete[] buf;
  return success;
}
// maybe add a check if(buf == NULL) { rsError("Unable to allocate buffer"); return false; }

bool writeImageToFilePPM(const rsImageF& R, const rsImageF& G, const rsImageF& B, 
  const char* path)
{
  int w = R.getWidth();
  int h = R.getHeight();
  rsAssert(G.getWidth() == w && G.getHeight() == h);
  rsAssert(B.getWidth() == w && B.getHeight() == h);
  unsigned char* buf = new unsigned char[w*h*3];
  if(buf == nullptr) {
    rsError("Unable to allocate buffer"); 
    return false; }
  for(int y = 0; y < h; y++) {
    for(int x = 0; x < w; x++) {
      int i = y*w*3 + x*3;
      buf[i+0] = (unsigned char) (255 * R(x, y));
      buf[i+1] = (unsigned char) (255 * G(x, y));
      buf[i+2] = (unsigned char) (255 * B(x, y));  }}
  bool success = writeImageToFilePPM(path, buf, w, h);
  delete[] buf;
  return success;
}

bool writeScaledImageToFilePPM(rsImageF& img, const char* path, int scl)
{
  // maybe factor out into a function magnify (or generally 
  // image.getResized(int newWidth, int newHeight, int interpolationMethod))
  int w = img.getWidth();
  int h = img.getHeight();
  rsImageF tmp(scl*w, scl*h);
  for(int x = 0; x < w; x++)  {
    for(int y = 0; y < h; y++) {
      for(int i = 0; i < scl; i++) {
        for(int j = 0; j < scl; j++) {
          tmp(scl*x+i, scl*y+j) = img(x, y); }}}}
  return writeImageToFilePPM(tmp, path);
}



void writeToMonoWaveFile(std::string path, float *signal, int numFrames, int sampleRate,
  int numBits)
{
  rosic::writeToMonoWaveFile(path.c_str(), signal, numFrames, sampleRate, numBits);
}

void writeToWaveFile(std::string path, const std::vector<std::complex<double>>& signal,
  int sampleRate)
{
  int N = (int)signal.size();
  std::vector<double> xL(N), xR(N);
  for(int n = 0; n < N; n++) {
    xL[n] = signal[n].real();
    xR[n] = signal[n].imag();
  }
  rosic::writeToStereoWaveFile(path.c_str(), &xL[0], &xR[0], N, sampleRate);
}
