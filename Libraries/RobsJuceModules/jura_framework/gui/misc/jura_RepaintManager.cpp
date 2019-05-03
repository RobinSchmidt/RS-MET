rsRepaintManager::rsRepaintManager()
{
  startTimerHz(30); // or maybe let this be controlled from ouside?
}

void rsRepaintManager::timerCallback()
{
  for(size_t i = 0; i < repaintees.size(); i++) {
    rsRepaintClient* rc = dynamic_cast<rsRepaintClient*>(repaintees[i]);
    if(rc != nullptr) {
      if(rc->needsRepaint())     // rsRepaintClients may suppress repainting conditionally
        repaintees[i]->repaint();
    }
    else
      repaintees[i]->repaint(); // all other Components always repaint
  }
}