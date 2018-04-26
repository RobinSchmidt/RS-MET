using namespace RSLib;

// Construction/Destruction:

rsEnvelopeFollower::rsEnvelopeFollower()
{
  mode = MEAN_ABS;
}

rsEnvelopeFollower::~rsEnvelopeFollower()
{

}

// Setup:

void rsEnvelopeFollower::setMode(int Mode)
{
  if( Mode >= MEAN_ABS && Mode < NUM_MODES )
    mode = Mode;
}