

void rsDualFilterPoly::computeCoeffs(double x, double y)
{
  // the filter equation is:
  // y1 = filter1(a*x + d*y2);
  // y2 = filter2(b*x + c*y1);
  // y  = e*y1 + f*y2;

  // we require:
  // top-left: only output of filter 1 is audible:
  // a = 1, b = 0, c = 0, d = 0, e = 1, f = 0
  // top-right: only filter 2:
  // a = 0, b = 1, c = 0, d = 0, e = 0, f = 1
  // bottom-left: filter1 -> filter2
  // a = 1, b = 0, c = 1, d = 0, e = 0, f = 1
  // bottom-right: filter2 -> filter1
  // a = 0, b = 1, c = 0, d = 1, e = 1, f = 0

}