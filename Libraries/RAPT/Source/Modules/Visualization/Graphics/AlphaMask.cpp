template<class TPix>
AlphaMask<TPix>::AlphaMask()
{
  size = 1;
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::setSize(double newSize)
{
  size = newSize;
  renderMask();
}

template<class TPix>
void AlphaMask<TPix>::renderMask()
{
  // something to do
}