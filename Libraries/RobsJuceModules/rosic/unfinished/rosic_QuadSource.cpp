
rsQuadSourcePoly::rsQuadSourcePoly()
{
  for(int i = 0; i < numSources; i++)
    sources[i] = new rsPolyModule;
  mixer = new rsVectorMixerPoly;
}

rsQuadSourcePoly::~rsQuadSourcePoly()
{
  for(int i = 0; i < numSources; i++)
    delete sources[i];
  delete mixer;
}