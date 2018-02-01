AudioModuleSelector::AudioModuleSelector(AudioModuleFactory* factoryToUse) 
  : RComboBox("ModuleSelector")
{
  setAudioModuleFactory(factoryToUse);
}

void AudioModuleSelector::setAudioModuleFactory(AudioModuleFactory* newFactory)
{
  moduleFactory = newFactory;
  RComboBox::clear();
  RTreeViewNode *parentNode = popUpMenu->getRootItem();
  RTreeViewNode *node;
  const std::vector<AudioModuleInfo>& infos = moduleFactory->getRegisteredModuleInfos();
  for(size_t i = 0; i < infos.size(); i++) {
    AudioModuleInfo info = infos[i];
    if(info.category == "")  { // uncategorized modules go into the top-level of the tree:
      node = new RTreeViewNode(info.type, int(i+1));
      popUpMenu->addTreeNodeItem(node); 
    }
    else {
      if(parentNode->getText() != info.category) {
        parentNode = popUpMenu->getItemByText(info.category);
        if(parentNode == nullptr) { // node for this category doesn't exist yet - add it:
          parentNode = new RTreeViewNode(info.category, -1, info.category);
          parentNode->setOpen(false);
          popUpMenu->addTreeNodeItem(parentNode);  
        }
      }
      parentNode->addChildNode(new RTreeViewNode(info.type, int(i+1)));
    }
  }
}

void AudioModuleSelector::drawHighlighted(bool shouldBeHighlighted)
{
  highlighted = shouldBeHighlighted;
  repaint();
}

void AudioModuleSelector::paint(Graphics& g)
{
  if(!highlighted)
    RComboBox::paint(g);
  else
  {
    g.fillAll(getHandleColour()); // only difference to baseclass version - maybe refactor
    g.setColour(getOutlineColour());
    g.drawRect(0, 0, getWidth(), getHeight(), 2);
    int x = 4;
    int y = getHeight()/2 - font->getFontAscent()/2;
    drawBitmapFontText(g, x, y, getSelectedItemText(), font, getTextColour());
  }
}
