#ifndef RS_RSCORE_H
#define RS_RSCORE_H


#include <RSBuildConfig.h> // import build settings from a (possibly customized) config file
                           // sitting in a directory listed in the compiler's include paths

#include "Definitions/Setup.h"
#include "Definitions/StandardHeaders.h"
#include "Definitions/ConstantDefinitions.h"
#include "Definitions/TypeDefinitions.h"
#include "Definitions/MacroDefinitions.h"

#include "Misc/ErrorHandling.h"

//#ifdef RS_DEBUG
//#include "Misc/GNUPlotCPP.h"  // so we can create plots when debugging
//#endif

#include "Utilities/NumberManipulations.h"
#include "Utilities/MathBasics.h"
#include "Utilities/MathBasics.inl"
#include "Utilities/ArrayFunctions.h"
#include "Utilities/ArrayFunctions.inl"
  // we need to include the .h before the .inl to avoid "was not declared in this scope [...]
  // declared here, later in translation unit" error in gcc
#include "Utilities/SortAndSearch.inl"
#include "Utilities/MatrixFunctions.h"
#include "Utilities/MatrixFunctions.inl"

#include "Containers/Array.h"
#include "Containers/Flags.h"
#include "Containers/Range.h"
#include "Containers/InfiniteDataStream.h"
#include "Containers/KeyValueMap.h"
#include "Containers/List.h"
#include "Containers/Tree.h"

#include "Text/String.h"
#include "Text/KeyValueStringPair.h"
#include "Text/KeyValueStringArray.h"
#include "Text/KeyValueStringTree.h"

#include "Files/File.h"
#include "Files/FileStream.h"
#include "Files/WaveFile.h"
#include "Files/FileInputOutput.h"

#include "Misc/Logger.h"
#include "Misc/SystemInformation.h"
#include "Misc/Callback.h"
#include "Misc/Parameter.h"

#endif
