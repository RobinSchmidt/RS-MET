#ifndef RS_STANDARDHEADERS_H
#define RS_STANDARDHEADERS_H


#if _MSC_VER
#pragma push_macro("_CRT_SECURE_NO_WARNINGS")  // by Vincent - what does this do?
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>      // STL algorithms
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#if _MSC_VER
#pragma pop_macro("_CRT_SECURE_NO_WARNINGS")
#endif


#endif
