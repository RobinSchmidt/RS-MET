rsRepaintManager::rsRepaintManager()
{
  startTimerHz(30); // or maybe let this be controlled from ouside?
}

void rsRepaintManager::timerCallback()
{
  for(size_t i = 0; i < repaintees.size(); i++)
  {
    repaintees[i]->repaint(); // todo: repaint only conditionally
  }
}