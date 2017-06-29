template<class T>
rsEllipse<T>::rsEllipse(T _scale, T _aspectRatio, T _angle, T _centerX, T _centerY)
  : scale(_scale), ratio(_aspectRatio), angle(_angle), centerX(_centerX), centerY(_centerY)
{
  updateCoeffs();
}

template<class T>
void rsEllipse<T>::setParameters(T newScale, T newAspectRatio, T newAngle, 
  T newCenterX, T newCenterY)
{
  scale   = newScale;
  ratio   = newAspectRatio;
  angle   = newAngle;
  centerX = newCenterX;
  centerY = newCenterY;
  updateCoeffs();
}

template<class T>
void rsEllipse<T>::getPointOnEllipse(T angle, T* x, T* y) const
{
  T s = sin(angle);
  T c = cos(angle);
  *x  = Axc*c + Axs*s + centerX;
  *y  = Ayc*c + Ays*s + centerY;
}

template<class T>
void rsEllipse<T>::getTangentCoeffsAt(T x, T y, T* A, T* B, T* C)
{
  // something to do
}

template<class T>
void rsEllipse<T>::updateCoeffs()
{
  // intermediate variables:
  T a  = sqrt(scale*ratio);
  T b  = sqrt(scale/ratio);
  T s  = sin(angle);
  T c  = cos(angle);
  T xc = centerX;    // shorthand for convenience
  T yc = centerY;
  T a2 = a*a;  
  T b2 = b*b;

  // parametric equation coeffs:
  Axc =  a*c;
  Axs = -b*s;
  Ayc =  a*s;
  Ays =  b*c;

  // implicit equation coeffs for conic:
  A = a2*s*s + b2*c*c;
  B = 2*(b2+a2)*s*c;
  C = a2*c*c + b2*s*s;
  D = -2*A*xc - B*yc;
  E = -2*C*yc - B*xc;
  F = A*xc*xc + B*xc*yc + C*yc*yc - a2*b2;
}
