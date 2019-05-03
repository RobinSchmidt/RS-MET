

void rsEllipseOscillator::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  updateOmega();
}

void rsEllipseOscillator::setSampleRate(double newSampleRate)
{
  omegaFactor = 2*PI/newSampleRate;
  updateOmega();
}

void updateOmega()
{

}