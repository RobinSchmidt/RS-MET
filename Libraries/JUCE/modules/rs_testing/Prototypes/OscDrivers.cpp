

template<class T, class TBlep>
void rsSyncOsc<T, TBlep>::reset()
{
  master.reset();
  slave.reset();
  blep.reset();
}



template<class T, class TBlep>
void rsDualBlepOsc<T, TBlep>::reset()
{
  osc1.reset();
  osc2.reset();
  blep1.reset();
  blep2.reset();
}