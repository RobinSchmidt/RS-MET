


The codebase ist structured as follows: 

Libraries:

In the "Libraries" folder, there is a JUCE and a RAPT folder. The RAPT folder 
contains the "Rob's Audio Processing Toolkit" library - source code, project 
files, documentation and miscellaneous other stuff. The JUCE folder contains a
copy of the "modules" folder from the original JUCE distribution. That's where 
all the relevant JUCE library code sits. In addition to the modules that are 
part of JUCE the distribution, there are additional modules with the prefix 
jura_. These are my own juce based modules which are used in the products. 
These modules provide for the glue between JUCE's GUI and I/O framework and 
RAPT's signal processing algorithms - you may guess now, where the name 
"jura" comes from.

Products:




