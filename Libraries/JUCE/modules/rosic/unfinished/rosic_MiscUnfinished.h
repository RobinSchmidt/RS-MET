#pragma once

namespace rosic
{

class rsEnvelopeExtractor
{

public:

  //void setInterpolationMode(int newMode);
  // linear, cubic

  void extractEnvelope(const double* input, int length, double* envelope);



  // convenience functions using settings suitable for various use cases:

  /** Function suitable for extracting the envelope of an extracted partial.. */
  static void extractPartialEnvelope(const double* input, int length, double* envelope);

protected:

  int mode = 0;



};




}