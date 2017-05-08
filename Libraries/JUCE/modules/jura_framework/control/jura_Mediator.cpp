
MediatedColleague::MediatedColleague() 
{
  mediator = NULL;
}

MediatedColleague::~MediatedColleague()
{
  if( mediator != NULL )
    mediator->deRegisterColleague(this);
}
    
// setup/inquiry/mediation:

void MediatedColleague::setMediator(Mediator *newMediator) 
{ 
  if( mediator != NULL )
    mediator->deRegisterColleague(this);
  mediator = newMediator; 
  if( mediator != NULL )
    mediator->registerColleague(this);
}

void MediatedColleague::notifyMediator(int messageCode)
{
  if( mediator != NULL )
    mediator->colleagueHasSentNotification(this, messageCode);
}

//=========================================================================================================================================
// class Mediator:

Mediator::Mediator() 
{

}

Mediator::~Mediator()
{
  for(int i=0; i<colleagues.size(); i++)
    colleagues.getUnchecked(i)->mediator = NULL;
}

// setup/mediation:

void Mediator::registerColleague(MediatedColleague *colleagueToRegister)
{
  colleagues.addIfNotAlreadyThere(colleagueToRegister);
}

void Mediator::deRegisterColleague(MediatedColleague *colleagueDeToRegister)
{
  colleagues.removeFirstMatchingValue(colleagueDeToRegister);
}

void Mediator::sendNotificationToColleagues(MediatedColleague *originatingColleague, 
  int messageCode)
{
  for(int i=0; i<colleagues.size(); i++)
    colleagues.getUnchecked(i)->mediatorHasSentNotification(originatingColleague, messageCode);
}

void Mediator::colleagueHasSentNotification(MediatedColleague *originatingColleague, 
  int messageCode)
{
  sendNotificationToColleagues(originatingColleague, messageCode);
}
