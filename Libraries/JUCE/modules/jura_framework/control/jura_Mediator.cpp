
MediatedColleague::~MediatedColleague()
{
  if( mediator != nullptr )
    mediator->deRegisterColleague(this);
}
    
// setup/inquiry/mediation:

void MediatedColleague::setMediator(Mediator *newMediator) 
{ 
  if( mediator != nullptr )
    mediator->deRegisterColleague(this);
  mediator = newMediator; 
  if( mediator != nullptr )
    mediator->registerColleague(this);
}

void MediatedColleague::notifyMediator(int messageCode)
{
  if( mediator != nullptr )
    mediator->colleagueHasSentNotification(this, messageCode);
}

//=================================================================================================
// class Mediator:

Mediator::Mediator() 
{

}

Mediator::~Mediator()
{
  for(int i=0; i<colleagues.size(); i++)
    colleagues[i]->mediator = nullptr;
}

// setup/mediation:

void Mediator::registerColleague(MediatedColleague *colleagueToRegister)
{
  RAPT::rsAppendIfNotAlreadyThere(colleagues, colleagueToRegister);
}

void Mediator::deRegisterColleague(MediatedColleague *colleagueDeToRegister)
{
  RAPT::rsRemoveFirstOccurrence(colleagues, colleagueDeToRegister);
}

void Mediator::sendNotificationToColleagues(MediatedColleague *originatingColleague, 
  int messageCode)
{
  for(size_t i = 0; i < colleagues.size(); i++)
    colleagues[i]->mediatorHasSentNotification(originatingColleague, messageCode);
}

void Mediator::colleagueHasSentNotification(MediatedColleague *originatingColleague, 
  int messageCode)
{
  sendNotificationToColleagues(originatingColleague, messageCode);
}
