//-------------------------------------------------------------------------------------------------
// construction/destruction:

EllipticSubBandFilter::EllipticSubBandFilter()
{
  initBiquadCoeffs(); // initializes the filter-coefficents
  reset();            // resets the filters memory buffers (to 0)
  numStages = 6;      // this is a 12th order filter with 6 biquad-sections

  // design the analog prototype filter with unit cutoff frequency:
  rsPrototypeDesigner designer;
  designer.setApproximationMethod(rsPrototypeDesigner::ELLIPTIC);
  designer.setOrder(2*numStages);
  designer.setPassbandRipple(0.1);
  designer.setStopbandRejection(96.0);
  designer.getPolesAndZeros(prototypePoles, prototypeZeros);

  // now we have the poles and zeros of the prototype filter stored in our member variables, we
  // now make a digital halfband filter from them:
  setSubDivision(2.0);
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void EllipticSubBandFilter::setSubDivision(double newSubDivision)
{
  if( newSubDivision > 1.0 )
    subDivision = newSubDivision;
  else
  {
    // use a neutral filter when the subdivision is <= 1.0:
    initBiquadCoeffs();
    return;
    //DEBUG_BREAK;
  }

  // create temporary arrays for poles and zeros and copy the prototype poles and zeros into 
  // these:
  Complex poles[6];
  Complex zeros[6];
  for(int i=0; i<6; i++)
  {
    poles[i] = prototypePoles[i];
    zeros[i] = prototypeZeros[i];
  }

  // calculate the desired cutoff frequency of the denormalized analog filter
  double wd = (0.9*PI)/subDivision;        
   // normalized digital radian frequency - this factor of 0.9 is kinda arbitrary, ideally we 
   // would want to set a factor such that the stopband begins exactly at PI/subDivision (whereas 
   // here we are dealing with the end of the passband)

  //wd = PI/subDivision;    // for FreqShifter test

  double sampleRate = 44100.0; 
  // actually, this should be irrelevant but we need it as dummy. \todo make some tests if this is really irrelevant

  double wa = 2.0*sampleRate*tan(0.5*wd);  // pre-warped analog radian frequency 

  // obtain the denormalized s-domain poles and zeros by a lowpass->lowpass transform (this is done
  // in place on the temporary arrays):

  double gainDummy;
  rsPoleZeroMapper::prototypeToAnalogLowpass(poles, 6, zeros, 6, &gainDummy, wa);

  // transform poles and zeros from s-domain to z-domain via bilinear transform:
  rsPoleZeroMapper::bilinearAnalogToDigital(poles, 6, zeros, 6, sampleRate, &gainDummy);

  // convert the pole-zero representation to biquad cascade cofficients:
  FilterCoefficientConverter::polesAndZerosToBiquadCascade(poles, 6, zeros, 6, b0, b1, b2, a1, a2, false);

  // normalize DC-gain to unity for each stage:
  FilterCoefficientConverter::normalizeBiquadStages(b0, b1, b2, a1, a2, 0.0, 6);

  // use unity gain:
  //gain = 1.0;
}


