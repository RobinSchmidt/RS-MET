using namespace RSLib;

// Construction/Destruction:

rsDoublePendulum::rsDoublePendulum()
{
  // init inherited variables:
  x = 0.0;
  y.setDimensionality(4);
  setState(0.5*PI, PI, 0.0, 0.0); // 1st arm horizontal, 2nd arm upright

  // set up user parameters:
  setLength1(1.0);
  setLength2(1.0);
  setMass1(1.0);
  setMass2(1.0);
  setStepSize(0.01);
}

// Setup:

void rsDoublePendulum::setLength1(double newLength)
{
  l1 = newLength;
}

void rsDoublePendulum::setLength2(double newLength)
{
  l2 = newLength;
}

void rsDoublePendulum::setMass1(double newMass)
{
  m1 = newMass;
}

void rsDoublePendulum::setMass2(double newMass)
{
  m2 = newMass;
}

void rsDoublePendulum::setStepSize(double newStepSize)
{
  h = newStepSize;
}

void rsDoublePendulum::setState(double theta1, double theta2, double momentum1, double momentum2)
{
  y[0] = theta1;
  y[1] = theta2;
  y[2] = momentum1;
  y[3] = momentum2;
}

// Inquiry:

void rsDoublePendulum::getState(double *t1, double *t2, double *p1, double *p2)
{
  *t1 = y[0];
  *t2 = y[1];
  *p1 = y[2];
  *p2 = y[3];
}

void rsDoublePendulum::getAngles(double *a1, double *a2)
{
  *a1 = y[0];
  *a2 = y[1];
}

void rsDoublePendulum::getMomenta(double *m1, double *m2)
{
  *m1 = y[2];
  *m2 = y[3];
}

void rsDoublePendulum::getPendulumCoordinates(double *x1, double *y1, double *x2, double *y2)
{
  double t1, t2, p1, p2;
  getState(&t1, &t2, &p1, &p2);
  *x1 =  l1 * sin(t1);
  *y1 = -l1 * cos(t1);
  *x2 =  l2 * sin(t2) + *x1;
  *y2 = -l2 * cos(t2) + *y1;
}

// Processing:

void rsDoublePendulum::updateState()
{
  stepRungeKutta4(h);
  y[0] = rsWrapToInterval(y[0], -PI, PI);
  y[1] = rsWrapToInterval(y[1], -PI, PI);
}

rsVector<double> rsDoublePendulum::f(const double &x, const rsVector<double> &y)
{
  rsVector<double> v(4);

  // formulas taken from here (at the bottom, the Hamiltionian formulation):
  // http://scienceworld.wolfram.com/physics/DoublePendulum.html

  double g  = 9.81;  // gravitational constant - todo: make a file that collects physical constants

  // extract variables from state vector:
  double t1 = y[0];  // theta 1
  double t2 = y[1];  // theta 2
  double p1 = y[2];  // momentum 1
  double p2 = y[3];  // momentum 2

  // some intermediate variables that depend only on user parameters (their computation can be
  // done outside later, for optimization):
  double l12      = l1*l1;        // l1^2
  double l22      = l2*l2;        // l2^2
  double l1l2     = l1*l2;        // l1*l2
  double l12l22t2 = 2*l1l2*l1l2;  // 2 * l1^2 * l2^2 
  l12l22t2        = 2*l12*l22;    // should stay the same - test in debugger

  // intermediate variables:
  double s  = sin(t1-t2);
  double c  = cos(t1-t2);
  double a  = m1 + m2*s*s;
  double C1 = (p1*p2*s) / (l1l2*a);
  double C2 = ((l22*m2*p1*p1 + l12*(m1+m2)*p2*p2 - l1l2*m2*p1*p2*c) * sin(2*(t1-t2))) 
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


