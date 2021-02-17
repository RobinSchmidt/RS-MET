This is a customized version of the audio plugin host that comes with the JUCE distribution. It can
be used for plugin debugging and as starting point for other audio applications. The .jucer file
is customized in that it is configured to build for hosting vst2 plugins (and vst3 has been 
disabled). To make that work, i had to check the appropriate boxes in the menu for the
juce_audio_processors module and add the vst2 sdk search path to the "Header Search Paths" field in 
the main configuration. For more details about how to support vst2 and/or vst3, see the 
BuildNotes.txt in the ToolChain folder (the preliminaries for building vst plugins and hosts are 
the same).

The app stores it's settings in:

C:\Users\[UserName]\Documents\[SomePluginGraph].filtergraph

Sometimes it's necessary to rename or delete or move the most recently saved file, because when the 
host loads it auotmatically on start and the filtergraph immediately crashes the host, that's the 
way to recover