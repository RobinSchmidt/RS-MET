/*
   AngelCode Scripting Library
   Copyright (c) 2003-2007 Andreas Jonsson

   This software is provided 'as-is', without any express or implied 
   warranty. In no event will the authors be held liable for any 
   damages arising from the use of this software.

   Permission is granted to anyone to use this software for any 
   purpose, including commercial applications, and to alter it and 
   redistribute it freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you 
      must not claim that you wrote the original software. If you use
      this software in a product, an acknowledgment in the product 
      documentation would be appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and 
      must not be misrepresented as being the original software.

   3. This notice may not be removed or altered from any source 
      distribution.

   The original version of this library can be located at:
   http://www.angelcode.com/angelscript/

   Andreas Jonsson
   andreas@angelcode.com
*/


//
// as_scriptfunction.h
//
// A container for a compiled script function
//



#ifndef AS_SCRIPTFUNCTION_H
#define AS_SCRIPTFUNCTION_H

#include "as_config.h"
#include "as_string.h"
#include "as_array.h"
#include "as_datatype.h"

BEGIN_AS_NAMESPACE

class asCScriptEngine;
class asCModule;

struct asSScriptVariable
{
	asCString name;
	asCDataType type;
	int stackOffset;
};

const int asFUNC_SYSTEM    = 0;
const int asFUNC_SCRIPT    = 1;
const int asFUNC_INTERFACE = 2;
const int asFUNC_IMPORTED  = 3;

struct asSSystemFunctionInterface;

class asCScriptFunction
{
public:
	asCScriptFunction(asCModule *mod);
	~asCScriptFunction();

	void AddVariable(asCString &name, asCDataType &type, int stackOffset);

	int GetSpaceNeededForArguments();
	int GetSpaceNeededForReturnValue();
	asCString GetDeclaration(asCScriptEngine *engine);
	int GetLineNumber(int programPosition);
	void ComputeSignatureId(asCScriptEngine *engine);

	int                          funcType;
	asCModule                   *module;
	asCString                    name;
	asCDataType                  returnType;
	asCArray<asCDataType>        parameterTypes;
	asCArray<int>                inOutFlags;
	int                          id;
	int                          scriptSectionIdx;
	asCArray<asDWORD>            byteCode;
	asCArray<asCObjectType*>     objVariableTypes;
	asCArray<int>	             objVariablePos;
	asCArray<int>                lineNumbers;
	int                          stackNeeded;
	bool                         isReadOnly;
	asCObjectType *              objectType;
	asCArray<asSScriptVariable*> variables;
	int                          signatureId;

	asSSystemFunctionInterface  *sysFuncIntf;
};

END_AS_NAMESPACE

#endif
