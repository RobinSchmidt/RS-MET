This is the codebase on which RS-MET products are based. If you want to use the
code in an open source project, feel free to do so (but please notify me and 
give proper credits and if you use the JUCE code, be sure to adhere to its 
licensing scheme, too). For closed source projects, you may purchase a 
commercial license. I negotiate the conditions individually, based on the size 
of the product/company, the role of my code within it, etc.


The codebase ist structured as follows: 

In the "Libraries" folder, there is a JUCE folder which contains a "modules"
subfolder. In addition to the modules that are part of the JUCE distribution, 
it contains some additional modules. These are my own ones and on these, 
all the actual products are based. The rapt module (Rob's Audio Processing 
Templates) is a template based library with rather low level DSP code. The 
rosic module (Rob's Signal Processing Classes) is a bit more high-level and 
more convenient and often already facilitates easy integration of the code 
into plugins. jura_framework is my JUCE-based GUI and plugin framework and 
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







