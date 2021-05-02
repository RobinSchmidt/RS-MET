
bool samplerEngineUnitTest()
{
  bool ok = true;

  using TSig = float;    // maybe use rsFloat64x2
  using TPar = float;    // maybe use double
  using TSmp = float;
  using SE   = rsSamplerEngine<TSig, TPar, TSmp>;

  SE se;

  // ToDo: 
  // -create a couple of simple samples (maybe sine-waves or something) and assign them to
  //  regions, trigger notes, record output and compare the produced output to what is expected
  // -set up performance tests, too


  return ok;
}