

void rsDualFilterPoly::computeCoeffs(double x, double y)
{
  // the filter equation is:
  // y1 = filter1(a*x + d*y2);    // output of 1st filter
  // y2 = filter2(b*x + c*y1);    // output of 2nd filter
  // y  = e*y1 + f*y2;            // total output is weighted sum

  // we require:
  // top-left: only output of filter 1 is audible:
  // a = 1, b = 0, c = 0, d = 0, e = 1, f = 0
  // top-right: only filter 2:
  // a = 0, b = 1, c = 0, d = 0, e = 0, f = 1
  // bottom-left: series: filter1 -> filter2 
  // a = 1, b = 0, c = 1, d = 0, e = 0, f = 1
  // bottom-right: filter2 -> filter1
  // a = 0, b = 1, c = 0, d = 1, e = 1, f = 0

  // The idea is to somehow have a fully general way of morphing between parallel and serial
  // configuration and also just having filter 1 or 2 audible...i don't know yet how to compute
  // the coeffs from user parameters but it would be nice if the configuration could also be done 
  // via vector a pad - maybe it should work like:
  // -top: parallel connection
  //  -left: only filter 1
  //  -right: only filter2
  //  -center: both filters equally loud
  // -bottom: serial connection
  //  -left: filter 1 first
  //  -right: filter 2 first
  //  -center: something in between - whatever that means
}