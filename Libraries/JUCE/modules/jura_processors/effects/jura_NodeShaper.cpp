NodeShaperAudioModule::NodeShaperAudioModule(CriticalSection *newPlugInLock)
  : ModulatableAudioModule(newPlugInLock)
{
  setModuleTypeName("NodeShaper");
  createParameters();
}

NodeShaperAudioModule::~NodeShaperAudioModule()
{

}

AudioModuleEditor* NodeShaperAudioModule::createEditor()
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
  setSize(480, 300);
}

NodeShaperModuleEditor::~NodeShaperModuleEditor()
{

}

void NodeShaperModuleEditor::createWidgets()
{

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
}
