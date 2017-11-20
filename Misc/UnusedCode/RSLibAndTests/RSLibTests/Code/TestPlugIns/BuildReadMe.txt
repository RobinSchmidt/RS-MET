To build a particular effect as test-plugin, you need to de-activate the "Excluded From Build"
setting for the corresponding .cpp file (in Source Files -> PlugIns) and activate the setting
for the .cpp for which it was formerly de-activated. This way, we can use the same project for
building all effects as plugin for testing purposes. The output-dll will always be named
TestPlugIn.dll (but you may rename it as you see fit). The unique IDs will be different for each
of the generated .dlls (they are determined in the code that you activate for compilation).
