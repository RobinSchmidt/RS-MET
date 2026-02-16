// Maybe move this to .cpp file. But then we will need explicit instantiations.

template<class T>
void rsPitchDitherHelpers<T>::calcCycleDistribution(
  T period, T* midLength, T* probShort, T* probMid)
{
  // Compute lengths:
  T floorLength = rsFloor(period);
  T fracLength  = period - floorLength;
  T L1, L2, L3;
  if(fracLength < T(0.5))
    L1 = floorLength - T(1);
  else
    L1 = floorLength;
  L2 = L1 + T(1);
  L3 = L2 + T(1);

  // Compute intermediates:
  T e1 = L1 - period;
  T e2 = L2 - period;
  T e3 = L3 - period;
  T m1 = e1*e1;
  T m2 = e2*e2;
  T m3 = e3*e3;
  T M  = T(0.25);
  T M1 = M - m1;
  T M2 = M - m2; 
  T M3 = M - m3;
  T S  = T(1) / (e3*(m1-m2) - e2*(m1-m3) + e1*(m2-m3));

  // Compute outputs:
  *midLength = L2;
  *probShort = (M2*e3 - M3*e2) * S;
  *probMid   = (M3*e1 - M1*e3) * S;
  //*probLong  = (M1*e2 - M2*e1) * S;  // Would be redundant. See Notes
  
  // Notes:
  // 
  // - We don't have a probLong parameter because that would be redundant. It would always be given
  //   by 1 - (probShort + probMid).
  //
  // - The derivation of these formulas can be found in the textfile TempSketchPad.txt in the 
  //   research repo. ToDo: clean the derivation up and put it into its own dedicated textfile here
  //   in the main repo!
}