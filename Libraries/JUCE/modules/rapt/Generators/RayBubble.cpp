template<class T>
rsRayBubble<T>::rsRayBubble()
{

}

// setup:


// processing:

template<class T>
void rsRayBubble<T>::getSampleFrame(T &x, T &y)
{

}

template<class T>
void rsRayBubble<T>::reset()
{
  xc = x0; 
  yc = y0;
  dx = speed * cos(angle);
  dy = speed * sin(angle);
}
