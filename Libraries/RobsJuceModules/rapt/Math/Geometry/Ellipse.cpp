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
  this->A = a2*s*s + b2*c*c;
  this->B = 2*(b2-a2)*s*c;
  this->C = a2*c*c + b2*s*s;
  this->D = -2*this->A*xc - this->B*yc;
  this->E = -2*this->C*yc - this->B*xc;
  this->F = this->A*xc*xc + this->B*xc*yc + this->C*yc*yc - a2*b2;
  // formulas from here: https://en.wikipedia.org/wiki/Ellipse#General_ellipse
  // for some reason, gcc on linux and mac needs the this-pointers here - why?
  // here's some info:
  // https://stackoverflow.com/questions/4643074/why-do-i-have-to-access-template-base-class-members-through-the-this-pointer
  // this really sucks!
}





/*

Maybe get inspirations from existing 2D geometry libraries such as:

https://github.com/OneLoneCoder/olcUTIL_Geometry2D

*/