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
// as_scriptfunction.cpp
//
// A container for a compiled script function
//



#include "as_config.h"
#include "as_scriptfunction.h"
#include "as_tokendef.h"
#include "as_scriptengine.h"
#include "as_callfunc.h"

BEGIN_AS_NAMESPACE

asCScriptFunction::asCScriptFunction(asCModule *mod)
{
	funcType    = -1;
	module      = mod; 
	objectType  = 0; 
	name        = ""; 
	isReadOnly  = false;
	stackNeeded = 0;
	sysFuncIntf = 0;
	signatureId = 0;
}

asCScriptFunction::~asCScriptFunction()
{
	for( asUINT n = 0; n < variables.GetLength(); n++ )
	{
		DELETE(variables[n],asSScriptVariable);
	}

	if( sysFuncIntf )
	{
		DELETE(sysFuncIntf,asSSystemFunctionInterface);
	}
}

int asCScriptFunction::GetSpaceNeededForArguments()
{
	// We need to check the size for each type
	int s = 0;
	for( asUINT n = 0; n < parameterTypes.GetLength(); n++ )
		s += parameterTypes[n].GetSizeOnStackDWords();

	return s;
}

int asCScriptFunction::GetSpaceNeededForReturnValue()
{
	return returnType.GetSizeOnStackDWords();
}

asCString asCScriptFunction::GetDeclaration(asCScriptEngine *)
{
	asCString str;

	str = returnType.Format();
	str += " ";
	if( objectType )
	{
		if( objectType->name != "" )
			str += objectType->name + "::";
		else
			str += "?::";
	}
	if( name == "" )
		str += "?(";
	else
		str += name + "(";

	if( parameterTypes.GetLength() > 0 )
	{
		asUINT n;
		for( n = 0; n < parameterTypes.GetLength() - 1; n++ )
		{
			str += parameterTypes[n].Format();
			if( parameterTypes[n].IsReference() && inOutFlags.GetLength() > n )
			{
				if( inOutFlags[n] == 1 ) str += "in";
				else if( inOutFlags[n] == 2 ) str += "out";
				else if( inOutFlags[n] == 3 ) str += "inout";
			}
			str += ", ";
		}

		str += parameterTypes[n].Format();
		if( parameterTypes[n].IsReference() && inOutFlags.GetLength() > n )
		{
			if( inOutFlags[n] == 1 ) str += "in";
			else if( inOutFlags[n] == 2 ) str += "out";
			else if( inOutFlags[n] == 3 ) str += "inout";
		}
	}

	str += ")";

	return str;
}

int asCScriptFunction::GetLineNumber(int programPosition)
{
	if( lineNumbers.GetLength() == 0 ) return 0;

	// Do a binary search in the buffer
	int max = (int)lineNumbers.GetLength()/2 - 1;
	int min = 0;
	int i = max/2;

	for(;;)
	{
		if( lineNumbers[i*2] < programPosition )
		{
			// Have we found the largest number < programPosition?
			if( max == i ) return lineNumbers[i*2+1];
			if( lineNumbers[i*2+2] > programPosition ) return lineNumbers[i*2+1];

			min = i + 1;
			i = (max + min)/2; 
		}
		else if( lineNumbers[i*2] > programPosition )
		{
			// Have we found the smallest number > programPosition?
			if( min == i ) return lineNumbers[i*2+1];

			max = i - 1;
			i = (max + min)/2;
		}
		else
		{
			// We found the exact position
			return lineNumbers[i*2+1];
		}
	}
}

void asCScriptFunction::AddVariable(asCString &name, asCDataType &type, int stackOffset)
{
	asSScriptVariable *var = NEW(asSScriptVariable);
	var->name = name;
	var->type = type;
	var->stackOffset = stackOffset;
	variables.PushLast(var);
}

void asCScriptFunction::ComputeSignatureId(asCScriptEngine *engine)
{
	// This function will compute the signatureId based on the 
	// function name, return type, and parameter types. The object 
	// type for methods is not used, so that class methods and  
	// interface methods match each other.
	for( asUINT n = 0; n < engine->signatureIds.GetLength(); n++ )
	{
		if( name != engine->signatureIds[n]->name ) continue;
		if( returnType != engine->signatureIds[n]->returnType ) continue;
		if( isReadOnly != engine->signatureIds[n]->isReadOnly ) continue;
		if( inOutFlags != engine->signatureIds[n]->inOutFlags ) continue;
		if( parameterTypes != engine->signatureIds[n]->parameterTypes ) continue;

		signatureId = engine->signatureIds[n]->signatureId;
		return;
	}

	signatureId = id;
	engine->signatureIds.PushLast(this);
}

END_AS_NAMESPACE

