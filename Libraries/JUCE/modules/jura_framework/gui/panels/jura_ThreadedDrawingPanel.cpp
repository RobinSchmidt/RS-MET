
//-------------------------------------------------------------------------------------------------
// construction/destruction:

ThreadedDrawingPanel::ThreadedDrawingPanel(TimeSliceThread* newThreadToUse) 
: ThreadedDrawingComponent(newThreadToUse)
{

}

void ThreadedDrawingPanel::resized()
{
  Panel::resized();
  ThreadedDrawingComponent::resized();
}

void ThreadedDrawingPanel::setDirty(bool shouldBeSetToDirty)
{
  Panel::setDirty(shouldBeSetToDirty);
  ThreadedDrawingComponent::setDirty(shouldBeSetToDirty);

}

void ThreadedDrawingPanel::paint(Graphics &g)
{
  ThreadedDrawingComponent::paint(g);
}