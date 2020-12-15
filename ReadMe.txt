This is the codebase on which RS-MET products are based. My main and most 
important codebase - my lifework, so to speak. If you want to use the
code in an open source project, feel free to do so (but please notify me and 
give proper credits and if you use any of the JUCE based code, be sure to 
adhere to its licensing scheme, too). For closed source projects, you may 
purchase a commercial license. I negotiate the conditions individually, based 
on the size of the product/company, the role of my code within it, etc.


The codebase ist structured as follows: 

The "Libraries" folder contains a "RobsJuceModules" subfolder which contains a
couple of my own JUCE modules, conforming to the way, JUCE itself is organized 
into modules. The rapt module (Rob's Audio Processing Templates) is a template 
based library with rather low level code for math, number crunching and signal 
processing. It has no dependencies whatsoever (not even on juce_core). The 
rosic module (Rob's Signal Processing Classes), which depends only on rapt, is 
a bit more high-level and more convenient to use and even includes some 
framework'ish stuff (like thread-synchronization, polyphonic voice-management, 
etc.) to facilitate easy integration of the code into plugins. rapt and rosic, 
although conforming to juce's module organization, do not depend in any way on 
juce. They can be used in their own right and/or combined with other 
frameworks. jura_framework is my juce-based GUI and plugin framework and 
jura_processors is the glue that ties together the DSP code from rapt/rosic 
with the jura_framework based GUI code into actual plugins or sub-modules of 
plugins (such as oscillators, filters, effects, etc.).

The most important project that can actually be built by itself is in:
Products/AudioPlugins/ToolChain/
ToolChain is a plugin that is actually many plugins in one. You can create a 
chain of several sound processors (which i internally call AudioModules) that 
were previously distributed as plugins in their own right. Project management 
is just sooo much easier when everything is lumped into a single project. The 
code in this project is trivial because all the actual code is in the library. 
The other projects that can be built are mostly for development and testing.


Disclaimer:
I'm currently in the process of restructuring the codebase, merging code from
3 different codebases with *lots* of overlapping functionality but slightly 
different goals and interfaces. That is to say: it's rather messy at the 
moment and the API is still subject to change. Eventually, my goal is to 
provide a commercially viable DSP library for licensing to audio software 
companies - but as said: i need to clean up a lot of things, so don't take the 
messy and inconsistent API too seriously yet - i'm still working on it (we are 
pre version 1.0 at the moment). ...and if you have some special requirement 
that the library does not yet support - consider to hire me to add it. I'm 
generally available for freelance work on audio DSP algorithms with special 
interest in musical DSP.

for a quick overview for what i'm planning for rapt library, you may take a 
look at this horribly incomplete, outdated, skeletal document:
https://github.com/RobinSchmidt/RS-MET/blob/master/Documentation/RAPT/LaTeX/UserManual/UserManual.pdf





