This is a customized version of the audio plugin host that comes with the JUCE distribution. It can
be used for plugin debugging and as starting point for other audio applications. The .jucer file
is customized in that it is configured to build for hosting vst2 plugins (and vst3 has been 
disabled). To make that work, i had to check the appropriate boxes in the menu for the
juce_audio_processors module and the following header search path:

../../../../../../RS-Private/Libaries/vstsdk2.4

has beed added to the "Header Search Paths" field in the main configuration. Obviously, the build 
relies on having that folder available. The RS-Private folder must sit right next to the main 
RS-MET folder and the long chain of "directory-ups" is due to the fact that the sdk files are 
included like this:

#include <pluginterfaces/vst2.x/aeffect.h>
#include <pluginterfaces/vst2.x/aeffectx.h>

in juce_VSTPluginFormat.cpp and the library search path is a relative path with respect to the 
directory where that file resides. However, as the name suggests, "RS-Private" is my private repo 
which contains stuff, i'm not allowed to redistribute and the vst2 sdk unfortunately belongs to 
these non-redistributable things. So to actually build it, you need to set your own vst2 sdk path 
there or create a dummy version of the RS-Private folder on your machine and put (a copy of) the 
vst2 sdk (which you have to obtain from somewhere) into the appropriate place and/or disable vst2 
support and enable vst3 instead (and add the appropriate vst3 search path pointing to your 
installation of the vst3 sdk).