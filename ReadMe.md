Welcome to the RS-MET Codebase
==============================

This is the codebase on which RS-MET (Robin Schmidt's Music Engineering Tools) products are based. 
My main and most important codebase. If you want to use the code in an open source project, feel
free to do so (but please notify me and give proper credits and if you use any of the JUCE based
code, be sure to adhere to its licensing scheme, too). For closed source projects, you may purchase
a commercial license. I negotiate the conditions individually, based on the size of the product 
and/or company, the role of my code within it, etc.


Repository Structure
--------------------

Shown below is a partial listing of the directory structure of the repo. There are many more folders
in the repo which are not shown, though. This is just an overview over the most important parts.

```
RS-MET/                        # Root folder of the repo
├── Libraries/                 # My own and 3rd party libraries
│   ├── JUCE/                  # Complete copy of the JUCE source tree (lib only, no add ons)
│   └── RobsJuceModules/       # My own libraries in the JUCE module format
│       ├── jura-framework/    # My JUCE-based GUI and plugin framework
│       ├── jura-processors/   # My audio plugin components with GUI and infrastructure
│       ├── rapt/              # Rob's Audio Processing Templates. Low level DSP/math/algorithms
│       └── rosic/             # Rob's Signal Processing Classes. Higher level plugin DSP components
├── Products/                  # JUCE-based buildable targets. Most are for internal use
│   ├── Applications/          # Standalone GUI apps
│   │   └── TestAppJURA/       # App containing unit tests for the jura-... stuff
│   └── AudioPlugins/          # VST/AU/... plugins
│       └── ToolChain/         # 👈 An all-in-one plugin. This is the MAIN GIG in this repo!
└── Tests/                     # Tests for the algorithms (unit tests, benchmarks, experiments, ...)
    └── TestsRosicAndRapt/     # Console app with tests for rapt and rosic
```

The repo is fully self contained so you don't need to worry about downloading any additional
dependencies. The "Libraries" folder contains a JUCE subfolder which contains a full copy of the 
source tree of the JUCE library. It's just the library itself without all the additional add-ons, 
example projects etc. The "RobsJuceModules" subfolder contains a couple of my own JUCE modules, 
conforming to the way, JUCE itself is organized into modules. The rapt module (Rob's Audio 
Processing Templates) is a template based library with rather low level code for math, number 
crunching and signal processing. It has no dependencies whatsoever (not even on juce_core). The
rosic module (Rob's Signal Processing Classes), which depends only on rapt, is a bit more high-level
and more convenient to use and even includes some framework'ish stuff that is out of the scope of 
rapt (such as thread-synchronization, polyphonic voice-management, etc.) to facilitate easy
integration of the code into plugins. The rapt and rosic libraries, although conforming to juce's
module organization, do not depend in any way on juce. They can be used in their own right and/or
combined with other frameworks. jura_framework is my juce-based GUI and plugin framework and
jura_processors is the glue that ties together the DSP code from rapt and rosic with the
jura_framework based GUI code into actual plugins or sub-modules of plugins (such as oscillators,
filters, effects, etc.). So that means that they depend on both rapt/rosic and juce. They are the
most high level modules and the ones from which actual plugins (with GUI and plugin-API plumbing) 
can be built.

The most important project that can actually be built by itself (i.e. is not just a library that is
supposed to be consumed by some other project) is ToolChain. It is a plugin that is actually many 
plugins in one. You can create a chain of several sound processors (which I internally call 
AudioModules) that were previously distributed as plugins in their own right. Project management is 
just sooo much easier when everything is lumped into a single project. The code in this project is 
trivial because all the actual code is in the library, more specifically, in the jura_processors 
module. This module is where all the high-level plugin code of ToolChain itself as well as all of 
its sub-plugins (aka "processors" aka "AudioModules") resides. The other projects that can be built 
are mostly for internal use, i.e. research, development, testing, debugging, etc. and should 
probably be ignored by people that are just interested in the ToolChain plugin which I assume to be 
the vast majority of visitors.


Disclaimer
----------

I'm currently in the process of restructuring the codebase, merging code from 3 different codebases
with *lots* of overlapping functionality but slightly different goals and interfaces. That is to
say: it's rather messy at the moment and the API is still subject to change. Eventually, my goal is
to provide a commercially viable DSP library for licensing to audio software companies - but as
said: i need to clean up a lot of things, so don't take the messy and inconsistent API too seriously
yet - i'm still working on it (we are pre version 1.0 at the moment). ...and if you have some
special requirement that the library does not yet support - consider to hire me to add it. I'm
generally available for freelance work on audio DSP algorithms with special interest in musical DSP.


<br><br><br><br>
----------------------------------------------------------------------------------------------------

#### ToDo

- Maybe make a section about ToolChain (with 2nd order headline) with some screenshots. They could
  be stored in a discussion thread in the GitHub repo in order to not bloat the repo itself. Explain
  how ToolChain itself can be viewed as a semi-modular synthesizer by letting the modules in
  different slots talk to one another via the modulation system. ToolChain has modules that can be
  used like full blown (effect- or instrument-) plugins in their own right but it has also simpler
  modules that make most sense in combination with other modules (like filters, oscillators, 
  envelope generators, etc.). Explain the ideas and concepts behind this.

- In the "Repository Structure" explain the dependencies. I already do to some extent but only 
  partially. Explain it more fully for ToolChain and the TestsRosicandRapt project. Maybe use a 
  top-down approach.

- Explain a bit the "TestsRosicAndRapt" project. It's the second most important one. It has all the
  unit tests for the (math, DSP, etc.) algorithms which makes it also kind of important. The 3rd
  most important project is probably the TestAppJURA. It has the unit tests for the higher level
  infrastructural stuff which is also kind of important. 
  
- Maybe also include a reference to the research repo. Some parts of the R&D I do also there and the
  criteria for what goes where are not very strict. They are roughly: What I suppose to eventually
  end up in the main repo, I may already initially develop within the main repo. The research repo
  also has a lot of "just for fun" stuff that will probably never make it into the main repo. 
  ...but who knows...

- Explain the branches. The most important ones are "master" and "work". Maybe rename "master" to
  "main" or "release" or "stable". It's supposed to be the cleanest branch. The "work" branch is 
  always the most up-to-date one, i.e. the one inside of which I usually work - with all the latest
  features but also the dirtiest one with the least thorough code review, testing and vetting. It's 
  currently set to be the default branch. Maybe include a warning about this or change it. ...but 
  I'm not sure how that will affect how GitHub will count my commit statistics.

- Maybe explain a bit why there is rapt and rosic with so much overlapping scope. The reasons are
  historical and it is expected that over time, much of the lower level stuff in rosic will be 
  promoted (or demoted?) into rapt. They may be "floated down" in the dependency chain.

- Maybe create a .md file with an overview over the most interesting DSP algorithms that are
  available inside the library. I started it here:
  https://github.com/RobinSchmidt/RS-MET/blob/work/Notes/Scratch/AlgorithmsOverview.md