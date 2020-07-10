//#include <JuceHeader.h>
using namespace RAPT;

void plotSmoothEnvWithVariousSustains(double att, double dec)
{
  int N    = 1200;
  int nOff = 800;     // note-off sample instant
  int key  = 64;
  int vel  = 64;

  std::vector<double> y(N);

  GNUPlotter plt;

  rsAttackDecayEnvelope<double> env;
  env.setAttackSamples(att);
  env.setDecaySamples(dec);

  for(int i = 0; i <= 5; i++)
  {
    env.setSustain(i*0.2);
    env.reset();
    env.noteOn(key, vel);
    for(int n = 0; n < nOff; n++)
      y[n] = env.getSample();
    env.noteOff(key, vel);  // why does note-off need key and vel?
    for(int n = nOff; n < N; n++)
      y[n] = env.getSample();
    plt.addDataArrays(N, &y[0]);
  }

  plt.plot();
}

template<class T>
class rsAttackDecayFilterTest : public rsAttackDecayFilter<T>
{

public:

  void setCoeffs(T ca, T cd, T s)
  {
    this->ca = ca;
    this->cd = cd;
    this->s  = s;
    this->coeffsDirty = false;
  }

  void setState(T ya, T yd)
  {
    this->ya = ya;
    this->yd = yd;
  }

  T getX()
  {
    //T x = T(1);  // initial guess - maybe try to find a better guess
    //T tol = 1.e-13;


    T R  = cd/ca;
    auto logR = [=](T x)->T{ return rsLogB(x, R); }; // logarithm to basis R

    // objective function of which we want to find a zero:
    auto f = [=](T x)->T
    {
      T y =  (x+yd) * pow(cd, logR(((x+ya)*log(ca))/((x+yd)*log(cd))))
           - (x+ya) * pow(ca, logR(((x+ya)*log(ca))/((x+yd)*log(cd))))
           - 1/s;
      // todo: try to simplify the expression

      return y;
    };

    T x = rsRootFinder<T>::bisection(f, T(0), T(1), T(0));
    return x;

  }
  // todo: find a suitable name

};

template<class T>
void plotAttDecResponse(T ca, T cd, T ya, T yd, T s, T x, int N = 500)
{
  // just a plot to verify a formula...

  rsAttackDecayFilterTest<T> flt;
  flt.setCoeffs(ca, cd, s);
  flt.setState(ya, yd);

  std::vector<T> y(N);
  y[0] = flt.getSample(x);
  for(int n = 1; n < N; n++)
    y[n] = flt.getSample(0);

  // compute peak location np:
  T a0 = x + ya;
  T d0 = x + yd;
  T D  = d0*log(cd);
  T A  = a0*log(ca);
  T R  = cd/ca;
  T np = rsLogB(A/D, R);

  // compute peak height ep:
  T dp = d0*pow(cd, np);  // decay output at peak
  T ap = a0*pow(ca, np);  // attack output at peak
  T ep = s*(dp-ap);       // env output at peak

  // compute peak height again, using the monster formula without intermediate variables:
  auto logR = [=](T x)->T{ return rsLogB(x, R); }; // logarithm to basis R
  T ep2 =   (x+yd) * pow(cd, logR(((x+ya)*log(ca))/((x+yd)*log(cd)))) 
          - (x+ya) * pow(ca, logR(((x+ya)*log(ca))/((x+yd)*log(cd))));
  ep2 *= s;
  // OK ep2 == ep, so the formula is correct. Now, we must set that formula equal to 1 and solve 
  // for x, if possible - otherwise come up with an iterative algorithm to solve the implicit 
  // equation for x. We need to find x, such that:
  // 1/s =   (x+yd) * cd^(logR(((x+ya)*log(ca))/((x+yd)*log(cd))))
  //       - (x+ya) * ca^(logR(((x+ya)*log(ca))/((x+yd)*log(cd))))

  rsPlotVector(y);


  // now with compensation - the resulting function should have a peak with height 1:
  flt.setCoeffs(ca, cd, s);
  flt.setState(ya, yd);
  T x2 = flt.getX(); 
  // this function really needs a better name - maybe getCompensatedInput, 
  // getAccumulationCompensation, ...

  y[0] = flt.getSample(x2);
  for(int n = 1; n < N; n++)
    y[n] = flt.getSample(0);
  rsPlotVector(y);
  // ...hmm - it's not exactly 1 but 0.995 ...maybe the ya,yd states are taken one sample too late 
  // or early? also, the bisection method is far too slow to converge for production code - but we
  // are getting close...almost there


  int dummy = 0;
}

void attackDecayEnvelope()
{

  plotAttDecResponse(0.95, 0.99, 0.1, 0.4, 2.0, 1.0);



  int N    = 1200;
  int nOff = 800;     // note-off sample instant
  int key  = 64;
  int vel  = 64;


  std::vector<double> y(N);


  rsAttackDecayFilter<double> flt;
  rsAttackDecayEnvelope<double> env; 
  // maybe rename to rsSmoothADSR, when release has been made independent from decay


  
  // plot DC response of filter:
  flt.setAttackSamples(5);
  flt.setDecaySamples(15);
  for(int n = 0; n < N; n++)
    y[n] = flt.getSample(1.0);
  //rsPlotVector(y);
  double dcGain = flt.getGainAtDC();  // should be 22.335... - yep, works

  /*
  // plot a simple attack/decay envelope without sustain:
  env.setAttackSamples(50);
  env.setDecaySamples(100);
  env.setSustain(0.0);
  env.noteOn(key, vel);
  for(int n = 0; n < nOff; n++)
    y[n] = env.getSample();
  env.noteOff(key, vel);  // why does note-off need key and vel?
  for(int n = nOff; n < N; n++)
    y[n] = env.getSample();
  //rsPlotVector(y);
  */


  // plot the response that we get when we fire a succession of several note-ons at it:
  env.setAttackSamples(20);
  env.setDecaySamples(100);
  int dt = 1;  // delta-t between the input impulses
  env.reset();
  for(int n = 0; n < N; n++)
  {
    if(n % dt == 0)
    {
      env.noteOn(key, vel);
      env.noteOff(key, vel);  // should not matter, if we call noteOff or not
    }
    y[n] = env.getSample();
  }
  dcGain = env.getGainAtDC();
  rsPlotVector(y);
  // i think, we should attempt that the curve approaches 1 in a sort of saturation curve, when we
  // send a not at each sample


  // plot a family of envelopes with sustain settings 0.0,0.2,0.4,0.6,0.8,1.0:
  //plotSmoothEnvWithVariousSustains(50, 100);

  int dummy = 0;

  // Observations:
  // -the sustain level works but using nonzero sustain slightly changes the attack-time and 
  //  maximum peak level: the peak gets higher and occurs later with increasing sustain, for 
  //  example, with attack = 20, decay = 100, sustain = 0.5, the actual peak occurs at 24 samples
  //  and has a height of around 1.07 (as per the settings, it should occur at sample 20 witha 
  //  height of 1), with sustain = 1, we get apeak at sample 33 with height 1.16
  // -we actually have some sort of smooth ADSR now in which D==R and the smoothness is meant in 
  //  the sense that there are not abrupt changes in slope - maybe we should now somehow lift the
  //  D==R restriction to get a full smooth ADSR - maybe just switch the decay-coeff depending on
  //  whether the note is on or off - the smoothness comes from the attack filter
  

  // Notes:
  // -i think, with sustain==0, the resulting envelope is infinitely smooth everywhere, i.e. 
  //  infinitely often differentiable. With nonzero sustain, there will be a discontinuity in the
  //  2nd derivative at the transition from sustain to release, so it's only second order smooth
  //  at this point (-> verify this)
}
/*
trying to derive a formula to scale the input impulse to the filter when it has not yet decayed 
away completely so as to still reach 1.0 as peak height instead of overshooting it due to 
accumulation:

  an = a0 * ca^n      attack filter output
  dn = d0 * cd^n      decay filter output
  en = s*(dn-an)      envelope output

where:
 
  a0 = x + ca*ya
  d0 = x + cd*yd

i think, we need to: 

  -find the derivative of en with respect to n
  -set it to zero
  -solve the equation for n
  -compute en at that found value of n
  -set it to 1
  -solve the equation for x




*/
