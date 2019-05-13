#pragma once

/** This file contains subclasses of classes to be tested. They provide some additonal 
functionality to facilitate certain tests that would be inconvenient of impossible without 
subclassing because they need access to protected members. */


class TestOscArrayPolyBlep : public rosic::rsOscArrayPolyBlep1
{

public:

  void initAmpArrays(double value = RS_NAN(double))
  {
    rsFill(ampsL, value);
    rsFill(ampsR, value);
  }

};