#include "../Soundtouch/include/SoundTouch.h"
using namespace soundtouch;

#include "../Soundtouch/source/SoundStretch/WavFile.h"
#include "../ElastiquePro/zPlane/include/zPlane/elastiqueProAPI.h"

#define PI 3.1415926535897932384626433832795

int min(int x, int y);
int max(int x, int y);
double pitchOffsetToFreqFactor(double pitchOffset);
void clearBuffer(float *buffer, int length);
void copyBuffer(float *source, float *destination, int length);
double random(double min, double max, int seed = -1);
double cubicFadeIn(double x);
double cubicFadeOut(double x);