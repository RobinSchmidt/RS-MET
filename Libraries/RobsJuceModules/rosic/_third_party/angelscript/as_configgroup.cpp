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
// as_configgroup.cpp
//
// This class holds configuration groups for the engine
//



#include "as_config.h"
#include "as_configgroup.h"
#include "as_scriptengine.h"

BEGIN_AS_NAMESPACE

asCConfigGroup::asCConfigGroup()
{
	refCount = 0;
	defaultAccess = true;
}

asCConfigGroup::~asCConfigGroup()
{
}

int asCConfigGroup::AddRef()
{
	refCount++;
	return refCount;
}

int asCConfigGroup::Release()
{
	// Don't delete the object here, the engine will delete the object when ready
	refCount--;
	return refCount;
}

asCObjectType *asCConfigGroup::FindType(const char *obj)
{
	for( asUINT n = 0; n < objTypes.GetLength(); n++ )
		if( objTypes[n]->name == obj )
			return objTypes[n];

	return 0;
}

void asCConfigGroup::RefConfigGroup(asCConfigGroup *group)
{
	if( group == this || group == 0 ) return;

	// Verify if the group is already referenced
	for( asUINT n = 0; n < referencedConfigGroups.GetLength(); n++ )
		if( referencedConfigGroups[n] == group )
			return;

	referencedConfigGroups.PushLast(group);
	group->AddRef();
}

bool asCConfigGroup::HasLiveObjects(asCScriptEngine * /*engine*/)
{
	for( asUINT n = 0; n < objTypes.GetLength(); n++ )
		if( objTypes[n]->refCount != 0 )
			return true;

	return false;
}

void asCConfigGroup::RemoveConfiguration(asCScriptEngine *engine)
{
	assert( refCount == 0 );

	asUINT n;

	// Remove global variables
	for( n = 0; n < globalProps.GetLength(); n++ )
	{
		for( asUINT m = 0; m < engine->globalProps.GetLength(); m++ )
		{
			if( engine->globalProps[m] == globalProps[n] )
			{
				DELETE(engine->globalProps[m],asCProperty);
				engine->globalProps[m] = 0;
			}
		}
	}

	// Remove global behaviours
	for( n = 0; n < globalBehaviours.GetLength(); n++ )
	{
		int id = engine->globalBehaviours.operators[globalBehaviours[n]+1];
		engine->globalBehaviours.operators[globalBehaviours[n]] = 0;
		engine->globalBehaviours.operators[globalBehaviours[n]+1] = 0;
	
		// Remove the system function as well
		engine->DeleteScriptFunction(id);
	}

	// Remove global functions
	for( n = 0; n < scriptFunctions.GetLength(); n++ )
	{
		for( asUINT m = 0; m < engine->scriptFunctions.GetLength(); m++ )
		{
			if( engine->scriptFunctions[m] == scriptFunctions[n] )
			{
				engine->DeleteScriptFunction(m);
			}
		}
	}

	// Remove object types
	for( n = 0; n < objTypes.GetLength(); n++ )
	{
		for( asUINT m = 0; m < engine->objectTypes.GetLength(); m++ )
		{
			if( engine->objectTypes[m] == objTypes[n] )
			{
				DELETE(engine->objectTypes[m],asCObjectType);
				engine->objectTypes[m] = 0;
			}
		}
	}

	// Release other config groups
	for( n = 0; n < referencedConfigGroups.GetLength(); n++ )
		referencedConfigGroups[n]->refCount--;
	referencedConfigGroups.SetLength(0);
}

int asCConfigGroup::SetModuleAccess(const char *module, bool hasAccess)
{
	if( module == asALL_MODULES )
	{
		// Set default module access
		defaultAccess = hasAccess;
	}
	else
	{
		asCString mod(module ? module : "");
		asSMapNode<asCString,bool> *cursor = 0;
		if( moduleAccess.MoveTo(&cursor, mod) )
		{
			moduleAccess.GetValue(cursor) = hasAccess;
		}
		else
		{
			moduleAccess.Insert(mod, hasAccess);
		}
	}

	return 0;
}

bool asCConfigGroup::HasModuleAccess(const char *module)
{
	asCString mod(module ? module : "");
	asSMapNode<asCString,bool> *cursor = 0;
	if( moduleAccess.MoveTo(&cursor, mod) )
		return moduleAccess.GetValue(cursor);
	
	return defaultAccess;
}

END_AS_NAMESPACE
