// Construction/Destruction:

template<class TSig, class TPar>
rsDoublePendulum<TSig, TPar>::rsDoublePendulum()
{
  // init inherited variables:
  this->x = 0.0;
  this->y.setDimensionality(4);
  setState(0.5*PI, PI, 0.0, 0.0); // 1st arm horizontal, 2nd arm upright

  // set up user parameters:
  setLength1(1.0);
  setLength2(1.0);
  setMass1(1.0);
  setMass2(1.0);
  setStepSize(0.01);
}

// Setup:

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::setLength1(TPar newLength)
{
  l1 = newLength;
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::setLength2(TPar newLength)
{
  l2 = newLength;
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::setMass1(TPar newMass)
{
  m1 = newMass;
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::setMass2(TPar newMass)
{
  m2 = newMass;
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::setStepSize(TPar newStepSize)
{
  h = newStepSize;
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::setState(TSig theta1, TSig theta2, TSig momentum1, TSig momentum2)
{
  this->y[0] = theta1;
  this->y[1] = theta2;
  this->y[2] = momentum1;
  this->y[3] = momentum2;
}

// Inquiry:

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::getState(TSig *t1, TSig *t2, TSig *p1, TSig *p2)
{
  *t1 = this->y[0];
  *t2 = this->y[1];
  *p1 = this->y[2];
  *p2 = this->y[3];
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::getAngles(TSig *a1, TSig *a2)
{
  *a1 = this->y[0];
  *a2 = this->y[1];
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::getMomenta(TSig *m1, TSig *m2)
{
  *m1 = this->y[2];
  *m2 = this->y[3];
}

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::getPendulumCoordinates(TSig *x1, TSig *y1, TSig *x2, TSig *y2)
{
  double t1, t2, p1, p2;
  getState(&t1, &t2, &p1, &p2);
  *x1 =  l1 * sin(t1);
  *y1 = -l1 * cos(t1);
  *x2 =  l2 * sin(t2) + *x1;
  *y2 = -l2 * cos(t2) + *y1;
}

// Processing:

template<class TSig, class TPar>
void rsDoublePendulum<TSig, TPar>::updateState()
{
  stepRungeKutta4(h);
  this->y[0] = rsWrapToInterval(this->y[0], -PI, PI);
  this->y[1] = rsWrapToInterval(this->y[1], -PI, PI);
}

template<class TSig, class TPar>
rsVector<TSig> rsDoublePendulum<TSig, TPar>::f(const TSig &x, const rsVector<TSig> &y)
{
  rsVector<TSig> v(4); // todo: avoid this memory allocation and copying...maybe let the
  // rsDifferentialEquationSystem baseclass use a std::function to compute the vector/array of
  // derivatives - it should get a reference (or pointer) to the currentt state and one (output)
  // reference for the array of derivatives

  // formulas taken from here (at the bottom, the Hamiltionian formulation):
  // http://scienceworld.wolfram.com/physics/DoublePendulum.html

  TSig g  = 9.81;  // gravitational constant - todo: make a file that collects physical constants

  // extract variables from state vector:
  TSig t1 = y[0];  // theta 1
  TSig t2 = y[1];  // theta 2
  TSig p1 = y[2];  // momentum 1
  TSig p2 = y[3];  // momentum 2

  // some intermediate variables that depend only on user parameters (their computation can be
  // done outside later, for optimization):
  TSig l12      = l1*l1;        // l1^2
  TSig l22      = l2*l2;        // l2^2
  TSig l1l2     = l1*l2;        // l1*l2
  TSig l12l22t2 = 2*l1l2*l1l2;  // 2 * l1^2 * l2^2
  l12l22t2        = 2*l12*l22;    // should stay the same - test in debugger

  // intermediate variables:
  TSig s  = sin(t1-t2);
  TSig c  = cos(t1-t2);
  TSig a  = m1 + m2*s*s;
  TSig C1 = (p1*p2*s) / (l1l2*a);
  TSig C2 = ((l22*m2*p1*p1 + l12*(m1+m2)*p2*p2 - l1l2*m2*p1*p2*c) * sin(2*(t1-t2)))
            / l12l22t2*a*a;

  //double C2 = (l12*m2*p1*p1 + l12*(m1+m2)*p2*p2 - l1l2*m2*p1*p2*c) * sin(2*(t1-t2)) / l12l22t2*a*a;
    // old - has mistakes - but sounded cool anyway

  // compute derivatives:
  v[0] = (l2*p1-l1*p2*c) / (l12*l2*a);                 // derivative with respect to theta1
  v[1] = (l1*(m1+m2)*p2 - l2*m2*p1*c) / (l1*l22*m2*a); // ...to theta2
  v[2] = -(m1+m2)*g*l1*sin(t1) - C1 + C2;              // ...to momentum 1
  v[3] = -m2*g*l2*sin(t2) + C1 - C2;                   // ...to momentum 2
  return v;
}

// Misc:


