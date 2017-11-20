
void RSLib::rsFadeOut(double *buffer, int start, int end)
{
  // rsAssert(end > start);
  int N = end - start;
  rsBreakpointModulator bm;
  bm.setSampleRate(1.0);
  bm.initialize();  
  bm.setBreakpointLevel(1, 0.0);
  bm.setBreakpointTime( 1, N);
  bm.setBreakpointShape(1, rsModBreakpoint::SMOOTH);  
  bm.noteOn();
  for(int n = start; n <= end; n++)
    buffer[n] *= bm.getSample();
}
