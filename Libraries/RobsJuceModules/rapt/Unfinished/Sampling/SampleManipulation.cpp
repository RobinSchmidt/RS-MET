template<class T>
void rsFadeOut(T* buffer, int start, int end)
{
  // rsAssert(end > start);
  int N = end - start;
  rsBreakpointModulator<T> bm;
  bm.setSampleRate(1.0);
  bm.initialize();  
  bm.setBreakpointLevel(1, 0.0);
  bm.setBreakpointTime( 1, N);
  bm.setBreakpointShape(1, rsModBreakpoint<T>::SMOOTH);  
  bm.noteOn();
  for(int n = start; n <= end; n++)
    buffer[n] *= bm.getSample();
}
