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
// as_module.cpp
//
// A class that holds a script module
//

#include "as_config.h"
#include "as_module.h"
#include "as_builder.h"
#include "as_context.h"

BEGIN_AS_NAMESPACE

asCModule::asCModule(const char *name, int id, asCScriptEngine *engine)
{
	this->name     = name;
	this->engine   = engine;
	this->moduleID = id;

	builder = 0;
	isDiscarded = false;
	isBuildWithoutErrors = false;
	contextCount = 0;
	moduleCount = 0;
	isGlobalVarInitialized = false;
	initFunction = 0;
}

asCModule::~asCModule()
{
	Reset();

	if( builder ) 
	{
		DELETE(builder,asCBuilder);
		builder = 0;
	}

	// Remove the module from the engine
	if( engine )
	{
		if( engine->lastModule == this )
			engine->lastModule = 0;

		int index = (moduleID >> 16);
		engine->scriptModules[index] = 0;
	}
}

int asCModule::AddScriptSection(const char *name, const char *code, int codeLength, int lineOffset, bool makeCopy)
{
	if( !builder )
		builder = NEW(asCBuilder)(engine, this);

	builder->AddCode(name, code, codeLength, lineOffset, (int)builder->scripts.GetLength(), makeCopy);

	return asSUCCESS;
}

int asCModule::Build()
{
	assert( contextCount == 0 );

 	Reset();

	if( !builder )
		return asSUCCESS;

	// Store the section names
	for( size_t n = 0; n < builder->scripts.GetLength(); n++ )
	{
		asCString *sectionName = NEW(asCString)(builder->scripts[n]->name);
		scriptSections.PushLast(sectionName);
	}

	// Compile the script
	int r = builder->Build();
	DELETE(builder,asCBuilder);
	builder = 0;
	
	if( r < 0 )
	{
		// Reset module again
		Reset();

		isBuildWithoutErrors = false;
		return r;
	}

	isBuildWithoutErrors = true;

	engine->PrepareEngine();

	CallInit();

	return r;
}

int asCModule::ResetGlobalVars()
{
	if( !isGlobalVarInitialized ) return asERROR;

	CallExit();

	CallInit();

	return 0;
}

void asCModule::CallInit()
{
	if( isGlobalVarInitialized ) return;

	memset(globalMem.AddressOf(), 0, globalMem.GetLength()*sizeof(asDWORD));

	if( initFunction && initFunction->byteCode.GetLength() == 0 ) return;

	int id = asFUNC_INIT;
	asIScriptContext *ctx = 0;
	int r = engine->CreateContext(&ctx, true);
	if( r >= 0 && ctx )
	{
		// TODO: Add error handling
		((asCContext*)ctx)->PrepareSpecial(id, this);
		ctx->Execute();
		ctx->Release();
		ctx = 0;
	}

	isGlobalVarInitialized = true;
}

void asCModule::CallExit()
{
	if( !isGlobalVarInitialized ) return;

	for( size_t n = 0; n < scriptGlobals.GetLength(); n++ )
	{
		if( scriptGlobals[n]->type.IsObject() )
		{
			void *obj = *(void**)(globalMem.AddressOf() + scriptGlobals[n]->index);
			if( obj )
			{
				asCObjectType *ot = scriptGlobals[n]->type.GetObjectType();

				if( ot->beh.release )
					engine->CallObjectMethod(obj, ot->beh.release);
				else
				{
					if( ot->beh.destruct )
						engine->CallObjectMethod(obj, ot->beh.destruct);

					engine->CallFree(ot, obj);
				}
			}
		}
	}

	isGlobalVarInitialized = false;
}

void asCModule::Discard()
{
	isDiscarded = true;
}

void asCModule::Reset()
{
	assert( !IsUsed() );

	CallExit();


	// Free global variables
	globalMem.SetLength(0);
	globalVarPointers.SetLength(0);

	isBuildWithoutErrors = true;
	isDiscarded = false;

	if( initFunction )
	{
		engine->DeleteScriptFunction(initFunction->id);
		initFunction = 0;
	}

	size_t n;
	for( n = 0; n < scriptFunctions.GetLength(); n++ )
		engine->DeleteScriptFunction(scriptFunctions[n]->id);
	scriptFunctions.SetLength(0);

	for( n = 0; n < importedFunctions.GetLength(); n++ )
	{
		DELETE(importedFunctions[n],asCScriptFunction);
	}
	importedFunctions.SetLength(0);

	// Release bound functions
	for( n = 0; n < bindInformations.GetLength(); n++ )
	{
		int oldFuncID = bindInformations[n].importedFunction;
		if( oldFuncID != -1 )
		{
			asCModule *oldModule = engine->GetModuleFromFuncId(oldFuncID);
			if( oldModule != 0 ) 
			{
				// Release reference to the module
				oldModule->ReleaseModuleRef();
			}
		}
	}
	bindInformations.SetLength(0);

	for( n = 0; n < stringConstants.GetLength(); n++ )
	{
		DELETE(stringConstants[n],asCString);
	}
	stringConstants.SetLength(0);

	for( n = 0; n < scriptGlobals.GetLength(); n++ )
	{
		DELETE(scriptGlobals[n],asCProperty);
	}
	scriptGlobals.SetLength(0);

	for( n = 0; n < scriptSections.GetLength(); n++ )
	{
		DELETE(scriptSections[n],asCString);
	}
	scriptSections.SetLength(0);

	for( n = 0; n < classTypes.GetLength(); n++ )
		classTypes[n]->refCount--;
	classTypes.SetLength(0);

	// Release all used object types
	for( n = 0; n < usedTypes.GetLength(); n++ )
		usedTypes[n]->refCount--;
	usedTypes.SetLength(0);

	// Release all config groups
	for( n = 0; n < configGroups.GetLength(); n++ )
		configGroups[n]->Release();
	configGroups.SetLength(0);
}

int asCModule::GetFunctionIDByName(const char *name)
{
	if( isBuildWithoutErrors == false )
		return asERROR;
	
	// TODO: Improve linear search
	// Find the function id
	int id = -1;
	for( size_t n = 0; n < scriptFunctions.GetLength(); n++ )
	{
		if( scriptFunctions[n]->name == name )
		{
			if( id == -1 )
				id = scriptFunctions[n]->id;
			else
				return asMULTIPLE_FUNCTIONS;
		}
	}

	if( id == -1 ) return asNO_FUNCTION;

	return id;
}

int asCModule::GetMethodIDByDecl(asCObjectType *ot, const char *decl)
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	return engine->GetMethodIDByDecl(ot, decl, this);
}

int asCModule::GetImportedFunctionCount()
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	return (int)importedFunctions.GetLength();
}

int asCModule::GetImportedFunctionIndexByDecl(const char *decl)
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	asCBuilder bld(engine, this);

	asCScriptFunction func(this);
	bld.ParseFunctionDeclaration(decl, &func);

	// TODO: Improve linear search
	// Search script functions for matching interface
	int id = -1;
	for( asUINT n = 0; n < importedFunctions.GetLength(); ++n )
	{
		if( func.name == importedFunctions[n]->name && 
			func.returnType == importedFunctions[n]->returnType &&
			func.parameterTypes.GetLength() == importedFunctions[n]->parameterTypes.GetLength() )
		{
			bool match = true;
			for( asUINT p = 0; p < func.parameterTypes.GetLength(); ++p )
			{
				if( func.parameterTypes[p] != importedFunctions[n]->parameterTypes[p] )
				{
					match = false;
					break;
				}
			}

			if( match )
			{
				if( id == -1 )
					id = n;
				else
					return asMULTIPLE_FUNCTIONS;
			}
		}
	}

	if( id == -1 ) return asNO_FUNCTION;

	return id;
}

int asCModule::GetFunctionCount()
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	return (int)scriptFunctions.GetLength();
}

int asCModule::GetFunctionIDByDecl(const char *decl)
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	asCBuilder bld(engine, this);

	asCScriptFunction func(this);
	int r = bld.ParseFunctionDeclaration(decl, &func);
	if( r < 0 )
		return asINVALID_DECLARATION;

	// TODO: Improve linear search
	// Search script functions for matching interface
	int id = -1;
	for( size_t n = 0; n < scriptFunctions.GetLength(); ++n )
	{
		if( scriptFunctions[n]->objectType == 0 && 
			func.name == scriptFunctions[n]->name && 
			func.returnType == scriptFunctions[n]->returnType &&
			func.parameterTypes.GetLength() == scriptFunctions[n]->parameterTypes.GetLength() )
		{
			bool match = true;
			for( size_t p = 0; p < func.parameterTypes.GetLength(); ++p )
			{
				if( func.parameterTypes[p] != scriptFunctions[n]->parameterTypes[p] )
				{
					match = false;
					break;
				}
			}

			if( match )
			{
				if( id == -1 )
					id = scriptFunctions[n]->id;
				else
					return asMULTIPLE_FUNCTIONS;
			}
		}
	}

	if( id == -1 ) return asNO_FUNCTION;

	return id;
}

int asCModule::GetGlobalVarCount()
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	return (int)scriptGlobals.GetLength();
}

int asCModule::GetGlobalVarIDByName(const char *name)
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	// Find the global var id
	int id = -1;
	for( size_t n = 0; n < scriptGlobals.GetLength(); n++ )
	{
		if( scriptGlobals[n]->name == name )
		{
			id = (int)n;
			break;
		}
	}

	if( id == -1 ) return asNO_GLOBAL_VAR;

	return moduleID | id;
}

int asCModule::GetGlobalVarIDByDecl(const char *decl)
{
	if( isBuildWithoutErrors == false )
		return asERROR;

	asCBuilder bld(engine, this);

	asCProperty gvar;
	bld.ParseVariableDeclaration(decl, &gvar);

	// TODO: Improve linear search
	// Search script functions for matching interface
	int id = -1;
	for( size_t n = 0; n < scriptGlobals.GetLength(); ++n )
	{
		if( gvar.name == scriptGlobals[n]->name && 
			gvar.type == scriptGlobals[n]->type )
		{
			id = (int)n;
			break;
		}
	}

	if( id == -1 ) return asNO_GLOBAL_VAR;

	return moduleID | id;
}

int asCModule::AddConstantString(const char *str, size_t len)
{
	//  The str may contain null chars, so we cannot use strlen, or strcmp, or strcpy
	asCString *cstr = NEW(asCString)(str, len);

	// TODO: Improve linear search
	// Has the string been registered before?
	for( size_t n = 0; n < stringConstants.GetLength(); n++ )
	{
		if( *stringConstants[n] == *cstr )
		{
			DELETE(cstr,asCString);
			return (int)n;
		}
	}

	// No match was found, add the string
	stringConstants.PushLast(cstr);

	return (int)stringConstants.GetLength() - 1;
}

const asCString &asCModule::GetConstantString(int id)
{
	return *stringConstants[id];
}

int asCModule::GetNextFunctionId()
{
	return engine->GetNextScriptFunctionId();
}

int asCModule::GetNextImportedFunctionId()
{
	return FUNC_IMPORTED | (asUINT)importedFunctions.GetLength();
}

int asCModule::AddScriptFunction(int sectionIdx, int id, const char *name, const asCDataType &returnType, asCDataType *params, int *inOutFlags, int paramCount, bool isInterface, asCObjectType *objType)
{
	assert(id >= 0);

	// Store the function information
	asCScriptFunction *func = NEW(asCScriptFunction)(this);
	func->funcType   = isInterface ? asFUNC_INTERFACE : asFUNC_SCRIPT;
	func->name       = name;
	func->id         = id;
	func->returnType = returnType;
	func->scriptSectionIdx = sectionIdx;
	for( int n = 0; n < paramCount; n++ )
	{
		func->parameterTypes.PushLast(params[n]);
		func->inOutFlags.PushLast(inOutFlags[n]);
	}
	func->objectType = objType;

	scriptFunctions.PushLast(func);
	engine->SetScriptFunction(func);

	// Compute the signature id
	if( objType )
		func->ComputeSignatureId(engine);

	return 0;
}

int asCModule::AddImportedFunction(int id, const char *name, const asCDataType &returnType, asCDataType *params, int *inOutFlags, int paramCount, int moduleNameStringID)
{
	assert(id >= 0);

	// Store the function information
	asCScriptFunction *func = NEW(asCScriptFunction)(this);
	func->funcType   = asFUNC_IMPORTED;
	func->name       = name;
	func->id         = id;
	func->returnType = returnType;
	for( int n = 0; n < paramCount; n++ )
	{
		func->parameterTypes.PushLast(params[n]);
		func->inOutFlags.PushLast(inOutFlags[n]);
	}
	func->objectType = 0;

	importedFunctions.PushLast(func);

	sBindInfo info;
	info.importedFunction = -1;
	info.importFrom = moduleNameStringID;
	bindInformations.PushLast(info);

	return 0;
}

asCScriptFunction *asCModule::GetScriptFunction(int funcID)
{
	return engine->scriptFunctions[funcID & 0xFFFF];
}

asCScriptFunction *asCModule::GetImportedFunction(int funcID)
{
	return importedFunctions[funcID & 0xFFFF];
}

asCScriptFunction *asCModule::GetSpecialFunction(int funcID)
{
	if( funcID & FUNC_IMPORTED )
		return importedFunctions[funcID & 0xFFFF];
	else
	{
		if( (funcID & 0xFFFF) == asFUNC_INIT )
			return initFunction;
		else if( (funcID & 0xFFFF) == asFUNC_STRING )
		{
			assert(false);
		}

		return engine->scriptFunctions[funcID & 0xFFFF];
	}
}

int asCModule::AllocGlobalMemory(int size)
{
	int index = (int)globalMem.GetLength();

	size_t *start = globalMem.AddressOf();

	globalMem.SetLength(index + size);

	// Update the addresses in the globalVarPointers
	for( size_t n = 0; n < globalVarPointers.GetLength(); n++ )
	{
		if( globalVarPointers[n] >= start && globalVarPointers[n] < (start+index) )
			globalVarPointers[n] = &globalMem[0] + (size_t(globalVarPointers[n]) - size_t(start))/sizeof(void*);
	}

	return index;
}

int asCModule::AddContextRef()
{
	ENTERCRITICALSECTION(criticalSection);
	int r = ++contextCount;
	LEAVECRITICALSECTION(criticalSection);
	return r;
}

int asCModule::ReleaseContextRef()
{
	ENTERCRITICALSECTION(criticalSection);
	int r = --contextCount;
	LEAVECRITICALSECTION(criticalSection);

	return r;
}

int asCModule::AddModuleRef()
{
	ENTERCRITICALSECTION(criticalSection);
	int r = ++moduleCount;
	LEAVECRITICALSECTION(criticalSection);
	return r;
}

int asCModule::ReleaseModuleRef()
{
	ENTERCRITICALSECTION(criticalSection);
	int r = --moduleCount;
	LEAVECRITICALSECTION(criticalSection);

	return r;
}

bool asCModule::CanDelete()
{
	// Don't delete if not discarded
	if( !isDiscarded ) return false;

	// Are there any contexts still referencing the module?
	if( contextCount ) return false;

	// If there are modules referencing this one we need to check for circular referencing
	if( moduleCount )
	{
		// Check if all the modules are without external reference
		asCArray<asCModule*> modules;
		if( CanDeleteAllReferences(modules) )
		{
			// Unbind all functions. This will break any circular references
			for( size_t n = 0; n < bindInformations.GetLength(); n++ )
			{
				int oldFuncID = bindInformations[n].importedFunction;
				if( oldFuncID != -1 )
				{
					asCModule *oldModule = engine->GetModuleFromFuncId(oldFuncID);
					if( oldModule != 0 ) 
					{
						// Release reference to the module
						oldModule->ReleaseModuleRef();
					}
				}
			}
		}

		// Can't delete the module yet because the  
		// other modules haven't released this one
		return false;
	}

	return true;
}

bool asCModule::CanDeleteAllReferences(asCArray<asCModule*> &modules)
{
	if( !isDiscarded ) return false;

	if( contextCount ) return false;

	modules.PushLast(this);

	// Check all bound functions for referenced modules
	for( size_t n = 0; n < bindInformations.GetLength(); n++ )
	{
		int funcID = bindInformations[n].importedFunction;
		asCModule *module = engine->GetModuleFromFuncId(funcID);

		// If the module is already in the list don't check it again
		bool inList = false;
		for( size_t m = 0; m < modules.GetLength(); m++ )
		{
			if( modules[m] == module )
			{
				inList = true;
				break;
			}
		}

		if( !inList )
		{
			bool ret = module->CanDeleteAllReferences(modules);
			if( ret == false ) return false;
		}
	}

	// If no module has external references then all can be deleted
	return true;
}


int asCModule::BindImportedFunction(int index, int sourceID)
{
	// Remove reference to old module
	int oldFuncID = bindInformations[index].importedFunction;
	if( oldFuncID != -1 )
	{
		asCModule *oldModule = engine->GetModuleFromFuncId(oldFuncID);
		if( oldModule != 0 ) 
		{
			// Release reference to the module
			oldModule->ReleaseModuleRef();
		}
	}

	if( sourceID == -1 )
	{
		bindInformations[index].importedFunction = -1;
		return asSUCCESS;
	}

	// Must verify that the interfaces are equal
	asCModule *srcModule = engine->GetModuleFromFuncId(sourceID);
	if( srcModule == 0 ) return asNO_MODULE;

	asCScriptFunction *dst = GetImportedFunction(index);
	if( dst == 0 ) return asNO_FUNCTION;

	asCScriptFunction *src = srcModule->GetScriptFunction(sourceID);
	if( src == 0 ) return asNO_FUNCTION;

	// Verify return type
	if( dst->returnType != src->returnType )
		return asINVALID_INTERFACE;

	if( dst->parameterTypes.GetLength() != src->parameterTypes.GetLength() )
		return asINVALID_INTERFACE;

	for( size_t n = 0; n < dst->parameterTypes.GetLength(); ++n )
	{
		if( dst->parameterTypes[n] != src->parameterTypes[n] )
			return asINVALID_INTERFACE;
	}

	// Add reference to new module
	srcModule->AddModuleRef();

	bindInformations[index].importedFunction = sourceID;

	return asSUCCESS;
}

const char *asCModule::GetImportedFunctionSourceModule(int index)
{
	if( index >= (int)bindInformations.GetLength() )
		return 0;

	index = bindInformations[index].importFrom;

	return stringConstants[index]->AddressOf();
}

bool asCModule::IsUsed()
{
	if( contextCount ) return true;
	if( moduleCount ) return true;

	return false;
}

asCObjectType *asCModule::GetObjectType(const char *type)
{
	// TODO: Improve linear search
	for( size_t n = 0; n < classTypes.GetLength(); n++ )
		if( classTypes[n]->name == type )
			return classTypes[n];

	return 0;
}

asCObjectType *asCModule::RefObjectType(asCObjectType *type)
{
	if( !type ) return 0;

	// Determine the index local to the module
	for( size_t n = 0; n < usedTypes.GetLength(); n++ )
		if( usedTypes[n] == type ) return type;

	usedTypes.PushLast(type);
	type->refCount++;

	RefConfigGroupForObjectType(type);

	return type;
}

void asCModule::RefConfigGroupForFunction(int funcId)
{
	// Find the config group where the function was registered
	asCConfigGroup *group = engine->FindConfigGroupForFunction(funcId);
	if( group == 0 )
		return;

	// Verify if the module is already referencing the config group
	for( size_t n = 0; n < configGroups.GetLength(); n++ )
	{
		if( configGroups[n] == group ) 
			return;
	}

	// Add reference to the group
	configGroups.PushLast(group);
	group->AddRef();
}

void asCModule::RefConfigGroupForGlobalVar(int gvarId)
{
	// Find the config group where the function was registered
	asCConfigGroup *group = engine->FindConfigGroupForGlobalVar(gvarId);
	if( group == 0 )
		return;

	// Verify if the module is already referencing the config group
	for( size_t n = 0; n < configGroups.GetLength(); n++ )
	{
		if( configGroups[n] == group ) 
			return;
	}

	// Add reference to the group
	configGroups.PushLast(group);
	group->AddRef();
}

void asCModule::RefConfigGroupForObjectType(asCObjectType *type)
{
	if( type == 0 ) return;

	// Find the config group where the function was registered
	asCConfigGroup *group = engine->FindConfigGroupForObjectType(type);
	if( group == 0 )
		return;

	// Verify if the module is already referencing the config group
	for( size_t n = 0; n < configGroups.GetLength(); n++ )
	{
		if( configGroups[n] == group ) 
			return;
	}

	// Add reference to the group
	configGroups.PushLast(group);
	group->AddRef();
}

int asCModule::GetGlobalVarIndex(int propIdx)
{
	void *ptr = 0;
	if( propIdx < 0 )
		ptr = engine->globalPropAddresses[-int(propIdx) - 1];
	else
		ptr = globalMem.AddressOf() + (propIdx & 0xFFFF);

	for( int n = 0; n < (signed)globalVarPointers.GetLength(); n++ )
		if( globalVarPointers[n] == ptr )
			return n;

	globalVarPointers.PushLast(ptr);
	return (int)globalVarPointers.GetLength()-1;
}

void asCModule::UpdateGlobalVarPointer(void *pold, void *pnew)
{
	for( asUINT n = 0; n < globalVarPointers.GetLength(); n++ )
		if( globalVarPointers[n] == pold )
		{
			globalVarPointers[n] = pnew;
			return;
		}
}

END_AS_NAMESPACE

