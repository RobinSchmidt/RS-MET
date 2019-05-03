NodeShaperAudioModule::NodeShaperAudioModule(CriticalSection *newPlugInLock)
  : ModulatableAudioModule(newPlugInLock)
{
  setModuleTypeName("NodeShaper");
  createParameters();
}

NodeShaperAudioModule::~NodeShaperAudioModule()
{

}

AudioModuleEditor* NodeShaperAudioModule::createEditor(int type)
{
  return new jura::NodeShaperModuleEditor(this);
}

void NodeShaperAudioModule::createParameters()
{

}

//=================================================================================================

// construction/destruction:

NodeShaperModuleEditor::NodeShaperModuleEditor(NodeShaperAudioModule* newNodeShaperAudioModule)
  : AudioModuleEditor(newNodeShaperAudioModule)
{
  jassert(newNodeShaperAudioModule != NULL ); // you must pass a valid module here
  nodeShaperModule = newNodeShaperAudioModule;
  createWidgets();
  //updateWidgetsAccordingToState();
  setSize(300, 300);
}

NodeShaperModuleEditor::~NodeShaperModuleEditor()
{

}

void NodeShaperModuleEditor::createWidgets()
{
  addWidget(nodeEditor = new rsNodeBasedFunctionEditor(
    &nodeShaperModule->mapper, nodeShaperModule->getCriticalSection()) );
}

void NodeShaperModuleEditor::updateWidgetsAccordingToState()
{

}

void NodeShaperModuleEditor::resized()
{
  AudioModuleEditor::resized();
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  nodeEditor->setBounds(x, y, w, h);
}
