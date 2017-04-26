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
// as_c.cpp
//
// A C interface to the library 
//

#include "as_config.h"

BEGIN_AS_NAMESPACE
#ifdef AS_C_INTERFACE

class asCBinaryStreamC : public asIBinaryStream
{
public:
	asCBinaryStreamC(asBINARYWRITEFUNC_t write, asBINARYREADFUNC_t read, void *param) {this->write = write; this->read = read; this->param = param;}

	void Write(const void *ptr, asUINT size) { write(ptr, size, param); }
	void Read(void *ptr, asUINT size) { read(ptr, size, param); }

	asBINARYREADFUNC_t read;
	asBINARYWRITEFUNC_t write;
	void *param;
};

int               asEngine_AddRef(asIScriptEngine *e)                                                                                                                                { return e->AddRef(); }
int               asEngine_Release(asIScriptEngine *e)                                                                                                                               { return e->Release(); }
int               asEngine_SetEngineProperty(asIScriptEngine *e, asDWORD property, asPWORD value)                                                                                    { return e->SetEngineProperty(property, value); }
asPWORD           asEngine_GetEngineProperty(asIScriptEngine *e, asDWORD property)                                                                                                   { return e->GetEngineProperty(property); }
int               asEngine_SetMessageCallback(asIScriptEngine *e, asFUNCTION_t callback, void *obj, asDWORD callConv)                                                                { return e->SetMessageCallback(asFUNCTION(callback), obj, callConv); }
int               asEngine_ClearMessageCallback(asIScriptEngine *e)                                                                                                                  { return e->ClearMessageCallback(); }
int               asEngine_RegisterObjectType(asIScriptEngine *e, const char *name, int byteSize, asDWORD flags)                                                                     { return e->RegisterObjectType(name, byteSize, flags); }
int               asEngine_RegisterObjectProperty(asIScriptEngine *e, const char *obj, const char *declaration, int byteOffset)                                                      { return e->RegisterObjectProperty(obj, declaration, byteOffset); }
int               asEngine_RegisterObjectMethod(asIScriptEngine *e, const char *obj, const char *declaration, asFUNCTION_t funcPointer, asDWORD callConv)                            { return e->RegisterObjectMethod(obj, declaration, asFUNCTION(funcPointer), callConv); }
int               asEngine_RegisterObjectBehaviour(asIScriptEngine *e, const char *datatype, asDWORD behaviour, const char *declaration, asFUNCTION_t funcPointer, asDWORD callConv) { return e->RegisterObjectBehaviour(datatype, behaviour, declaration, asFUNCTION(funcPointer), callConv); }
int               asEngine_RegisterGlobalProperty(asIScriptEngine *e, const char *declaration, void *pointer)                                                                        { return e->RegisterGlobalProperty(declaration, pointer); }
int               asEngine_RegisterGlobalFunction(asIScriptEngine *e, const char *declaration, asFUNCTION_t funcPointer, asDWORD callConv)                                           { return e->RegisterGlobalFunction(declaration, asFUNCTION(funcPointer), callConv); }
int               asEngine_RegisterGlobalBehaviour(asIScriptEngine *e, asDWORD behaviour, const char *declaration, asFUNCTION_t funcPointer, asDWORD callConv)                       { return e->RegisterGlobalBehaviour(behaviour, declaration, asFUNCTION(funcPointer), callConv); }
int               asEngine_RegisterInterface(asIScriptEngine *e, const char *name)                                                                                                   { return e->RegisterInterface(name); }
int               asEngine_RegisterInterfaceMethod(asIScriptEngine *e, const char *intf, const char *declaration)                                                                    { return e->RegisterInterfaceMethod(intf, declaration); }
int               asEngine_RegisterStringFactory(asIScriptEngine *e, const char *datatype, asFUNCTION_t factoryFunc, asDWORD callConv)                                               { return e->RegisterStringFactory(datatype, asFUNCTION(factoryFunc), callConv); }
int               asEngine_BeginConfigGroup(asIScriptEngine *e, const char *groupName)                                                                                               { return e->BeginConfigGroup(groupName); }
int               asEngine_EndConfigGroup(asIScriptEngine *e)                                                                                                                        { return e->EndConfigGroup(); }
int               asEngine_RemoveConfigGroup(asIScriptEngine *e, const char *groupName)                                                                                              { return e->RemoveConfigGroup(groupName); }
int               asEngine_SetConfigGroupModuleAccess(asIScriptEngine *e, const char *groupName, const char *module, bool haveAccess)                                                { return e->SetConfigGroupModuleAccess(groupName, module, haveAccess); }
int               asEngine_AddScriptSection(asIScriptEngine *e, const char *module, const char *name, const char *code, int codeLength, int lineOffset)                              { return e->AddScriptSection(module, name, code, codeLength, lineOffset); }
int               asEngine_Build(asIScriptEngine *e, const char *module)                                                                                                             { return e->Build(module); }
int               asEngine_Discard(asIScriptEngine *e, const char *module)                                                                                                           { return e->Discard(module); }
int               asEngine_GetFunctionCount(asIScriptEngine *e, const char *module)                                                                                                  { return e->GetFunctionCount(module); }
int               asEngine_GetFunctionIDByIndex(asIScriptEngine *e, const char *module, int index)                                                                                   { return e->GetFunctionIDByIndex(module, index); }
int               asEngine_GetFunctionIDByName(asIScriptEngine *e, const char *module, const char *name)                                                                             { return e->GetFunctionIDByName(module, name); }
int               asEngine_GetFunctionIDByDecl(asIScriptEngine *e, const char *module, const char *decl)                                                                             { return e->GetFunctionIDByDecl(module, decl); }
const char *      asEngine_GetFunctionDeclaration(asIScriptEngine *e, int funcID, int *length)                                                                                       { return e->GetFunctionDeclaration(funcID, length); }
const char *      asEngine_GetFunctionName(asIScriptEngine *e, int funcID, int *length)                                                                                              { return e->GetFunctionName(funcID, length); }
const char *      asEngine_GetFunctionModule(asIScriptEngine *e, int funcID, int *length)                                                                                            { return e->GetFunctionModule(funcID, length); }
const char *      asEngine_GetFunctionSection(asIScriptEngine *e, int funcID, int *length)                                                                                           { return e->GetFunctionSection(funcID, length); }
int               asEngine_GetMethodCount(asIScriptEngine *e, int typeId)																				                             { return e->GetMethodCount(typeId); }
int               asEngine_GetMethodIDByIndex(asIScriptEngine *e, int typeId, int index)															                                 { return e->GetMethodIDByIndex(typeId, index); }
int               asEngine_GetMethodIDByName(asIScriptEngine *e, int typeId, const char *name)														                               	 { return e->GetMethodIDByName(typeId, name); }
int               asEngine_GetMethodIDByDecl(asIScriptEngine *e, int typeId, const char *decl)															                             { return e->GetMethodIDByDecl(typeId, decl); }
int               asEngine_GetGlobalVarCount(asIScriptEngine *e, const char *module)                                                                                                 { return e->GetGlobalVarCount(module); }
int               asEngine_GetGlobalVarIDByIndex(asIScriptEngine *e, const char *module, int index)                                                                                  { return e->GetGlobalVarIDByIndex(module, index); }
int               asEngine_GetGlobalVarIDByName(asIScriptEngine *e, const char *module, const char *name)                                                                            { return e->GetGlobalVarIDByName(module, name); }
int               asEngine_GetGlobalVarIDByDecl(asIScriptEngine *e, const char *module, const char *decl)                                                                            { return e->GetGlobalVarIDByDecl(module, decl); }
const char *      asEngine_GetGlobalVarDeclaration(asIScriptEngine *e, int gvarID, int *length)                                                                                      { return e->GetGlobalVarDeclaration(gvarID, length); }
const char *      asEngine_GetGlobalVarName(asIScriptEngine *e, int gvarID, int *length)                                                                                             { return e->GetGlobalVarName(gvarID, length); }
void *            asEngine_GetGlobalVarPointer(asIScriptEngine *e, int gvarID)                                                                                                       { return e->GetGlobalVarPointer(gvarID); }
int               asEngine_GetImportedFunctionCount(asIScriptEngine *e, const char *module)                                                                                          { return e->GetImportedFunctionCount(module); }
int               asEngine_GetImportedFunctionIndexByDecl(asIScriptEngine *e, const char *module, const char *decl)                                                                  { return e->GetImportedFunctionIndexByDecl(module, decl); }
const char *      asEngine_GetImportedFunctionDeclaration(asIScriptEngine *e, const char *module, int importIndex, int *length)                                                      { return e->GetImportedFunctionDeclaration(module, importIndex, length); }
const char *      asEngine_GetImportedFunctionSourceModule(asIScriptEngine *e, const char *module, int importIndex, int *length)                                                     { return e->GetImportedFunctionSourceModule(module, importIndex, length); }
int               asEngine_BindImportedFunction(asIScriptEngine *e, const char *module, int importIndex, int funcID)                                                                 { return e->BindImportedFunction(module, importIndex, funcID); }
int               asEngine_UnbindImportedFunction(asIScriptEngine *e, const char *module, int importIndex)                                                                           { return e->UnbindImportedFunction(module, importIndex); }
int               asEngine_BindAllImportedFunctions(asIScriptEngine *e, const char *module)                                                                                          { return e->BindAllImportedFunctions(module); }
int               asEngine_UnbindAllImportedFunctions(asIScriptEngine *e, const char *module)                                                                                        { return e->UnbindAllImportedFunctions(module); }
int               asEngine_GetTypeIdByDecl(asIScriptEngine *e, const char *module, const char *decl)                                                                                 { return e->GetTypeIdByDecl(module, decl); }
const char *      asEngine_GetTypeDeclaration(asIScriptEngine *e, int typeId, int *length)                                                                                           { return e->GetTypeDeclaration(typeId, length); }
int               asEngine_GetSizeOfPrimitiveType(asIScriptEngine *e, int typeId)                                                                                                    { return e->GetSizeOfPrimitiveType(typeId); }
int               asEngine_SetDefaultContextStackSize(asIScriptEngine *e, asUINT initial, asUINT maximum)                                                                            { return e->SetDefaultContextStackSize(initial, maximum); }
asIScriptContext *asEngine_CreateContext(asIScriptEngine *e)                                                                                                                         { return e->CreateContext(); }
void *            asEngine_CreateScriptObject(asIScriptEngine *e, int typeId)                                                                                                        { return e->CreateScriptObject(typeId); }
void *            asEngine_CreateScriptObjectCopy(asIScriptEngine *e, void *obj, int typeId)                                                                                         { return e->CreateScriptObjectCopy(obj, typeId); }
void              asEngine_CopyScriptObject(asIScriptEngine *e, void *dstObj, void *srcObj, int typeId)                                                                              { e->CopyScriptObject(dstObj, srcObj, typeId); }
void              asEngine_ReleaseScriptObject(asIScriptEngine *e, void *obj, int typeId)                                                                                            { e->ReleaseScriptObject(obj, typeId); }
void              asEngine_AddRefScriptObject(asIScriptEngine *e, void *obj, int typeId)                                                                                             { e->AddRefScriptObject(obj, typeId); }
bool              asEngine_IsHandleCompatibleWithObject(asIScriptEngine *e, void *obj, int objTypeId, int handleTypeId)                                                              { return e->IsHandleCompatibleWithObject(obj, objTypeId, handleTypeId); }
int               asEngine_CompareScriptObjects(asIScriptEngine *e, bool &result, int behaviour, void *leftObj, void *rightObj, int typeId)                                          { return e->CompareScriptObjects(result, behaviour, leftObj, rightObj, typeId); }
int               asEngine_ExecuteString(asIScriptEngine *e, const char *module, const char *script, asIScriptContext **ctx, asDWORD flags)                                          { return e->ExecuteString(module, script, ctx, flags); }
int               asEngine_GarbageCollect(asIScriptEngine *e, bool doFullCycle)                                                                                                      { return e->GarbageCollect(doFullCycle); }
int               asEngine_GetObjectsInGarbageCollectorCount(asIScriptEngine *e)                                                                                                     { return e->GetObjectsInGarbageCollectorCount(); }
void              asEngine_NotifyGarbageCollectorOfNewObject(asIScriptEngine *e, void *obj, int typeId)                                                                              { e->NotifyGarbageCollectorOfNewObject(obj, typeId); }
void              asEngine_GCEnumCallback(asIScriptEngine *e, void *obj)                                                                                                             { e->GCEnumCallback(obj); }
int               asEngine_SaveByteCode(asIScriptEngine *e, const char *module, asBINARYWRITEFUNC_t outFunc, void *outParam)                                                         { asCBinaryStreamC out(outFunc, 0, outParam); return e->SaveByteCode(module, &out); }
int               asEngine_LoadByteCode(asIScriptEngine *e, const char *module, asBINARYREADFUNC_t inFunc, void *inParam)                                                            { asCBinaryStreamC in(0, inFunc, inParam); return e->LoadByteCode(module, &in); }

int              asContext_AddRef(asIScriptContext *c)                                                         { return c->AddRef(); }
int              asContext_Release(asIScriptContext *c)                                                        { return c->Release(); }
asIScriptEngine *asContext_GetEngine(asIScriptContext *c)                                                      { return c->GetEngine(); }
int              asContext_GetState(asIScriptContext *c)                                                       { return c->GetState(); }
int              asContext_Prepare(asIScriptContext *c, int funcID)                                            { return c->Prepare(funcID); }
int              asContext_SetArgByte(asIScriptContext *c, asUINT arg, asBYTE value)                           { return c->SetArgByte(arg, value); } 
int              asContext_SetArgWord(asIScriptContext *c, asUINT arg, asWORD value)                           { return c->SetArgWord(arg, value); } 
int              asContext_SetArgDWord(asIScriptContext *c, asUINT arg, asDWORD value)                         { return c->SetArgDWord(arg, value); } 
int              asContext_SetArgQWord(asIScriptContext *c, asUINT arg, asQWORD value)                         { return c->SetArgQWord(arg, value); }
int              asContext_SetArgFloat(asIScriptContext *c, asUINT arg, float value)                           { return c->SetArgFloat(arg, value); }
int              asContext_SetArgDouble(asIScriptContext *c, asUINT arg, double value)                         { return c->SetArgDouble(arg, value); }
int              asContext_SetArgAddress(asIScriptContext *c, asUINT arg, void *addr)						   { return c->SetArgAddress(arg, addr); }
int              asContext_SetArgObject(asIScriptContext *c, asUINT arg, void *obj)                            { return c->SetArgObject(arg, obj); }
void *           asContext_GetArgPointer(asIScriptContext *c, asUINT arg)                                      { return c->GetArgPointer(arg); }
int              asContext_SetObject(asIScriptContext *c, void *obj)										   { return c->SetObject(obj); }
asBYTE           asContext_GetReturnByte(asIScriptContext *c)                                                  { return c->GetReturnByte(); }
asWORD           asContext_GetReturnWord(asIScriptContext *c)                                                  { return c->GetReturnWord(); }
asDWORD          asContext_GetReturnDWord(asIScriptContext *c)                                                 { return c->GetReturnDWord(); }
asQWORD          asContext_GetReturnQWord(asIScriptContext *c)                                                 { return c->GetReturnQWord(); }
float            asContext_GetReturnFloat(asIScriptContext *c)                                                 { return c->GetReturnFloat(); }
double           asContext_GetReturnDouble(asIScriptContext *c)                                                { return c->GetReturnDouble(); }
void *           asContext_GetReturnAddress(asIScriptContext *c)											   { return c->GetReturnAddress(); }
void *           asContext_GetReturnObject(asIScriptContext *c)                                                { return c->GetReturnObject(); }
void *           asContext_GetReturnPointer(asIScriptContext *c)                                               { return c->GetReturnPointer(); }
int              asContext_Execute(asIScriptContext *c)                                                        { return c->Execute(); }
int              asContext_Abort(asIScriptContext *c)                                                          { return c->Abort(); }
int              asContext_Suspend(asIScriptContext *c)                                                        { return c->Suspend(); }
int              asContext_GetCurrentLineNumber(asIScriptContext *c, int *column)                              { return c->GetCurrentLineNumber(column); }
int              asContext_GetCurrentFunction(asIScriptContext *c)                                             { return c->GetCurrentFunction(); }
int              asContext_SetException(asIScriptContext *c, const char *string)                               { return c->SetException(string); }
int              asContext_GetExceptionLineNumber(asIScriptContext *c, int *column)                            { return c->GetExceptionLineNumber(column); }
int              asContext_GetExceptionFunction(asIScriptContext *c)                                           { return c->GetExceptionFunction(); }
const char *     asContext_GetExceptionString(asIScriptContext *c, int *length)                                { return c->GetExceptionString(length); }
int              asContext_SetLineCallback(asIScriptContext *c, asUPtr callback, void *obj, int callConv)      { return c->SetLineCallback(callback, obj, callConv); }
void             asContext_ClearLineCallback(asIScriptContext *c)                                              { c->ClearLineCallback(); }
int              asContext_SetExceptionCallback(asIScriptContext *c, asUPtr callback, void *obj, int callConv) { return c->SetExceptionCallback(callback, obj, callConv); }
void             asContext_ClearExceptionCallback(asIScriptContext *c)                                         { c->ClearExceptionCallback(); }
int              asContext_GetCallstackSize(asIScriptContext *c)                                               { return c->GetCallstackSize(); }
int              asContext_GetCallstackFunction(asIScriptContext *c, int index)                                { return c->GetCallstackFunction(index); }
int              asContext_GetCallstackLineNumber(asIScriptContext *c, int index, int *column)                 { return c->GetCallstackLineNumber(index, column); }
int              asContext_GetVarCount(asIScriptContext *c, int stackLevel)                                    { return c->GetVarCount(stackLevel); }
const char *     asContext_GetVarName(asIScriptContext *c, int varIndex, int *length, int stackLevel)          { return c->GetVarName(varIndex, length, stackLevel); }
const char *     asContext_GetVarDeclaration(asIScriptContext *c, int varIndex, int *length, int stackLevel)   { return c->GetVarDeclaration(varIndex, length, stackLevel); }
int              asContext_GetVarTypeId(asIScriptContext *c, int varIndex, int stackLevel)                     { return c->GetVarTypeId(varIndex, stackLevel); }
void *           asContext_GetVarPointer(asIScriptContext *c, int varIndex, int stackLevel)                    { return c->GetVarPointer(varIndex, stackLevel); }
int              asContext_GetThisTypeId(asIScriptContext *c, int stackLevel)                                  { return c->GetThisTypeId(stackLevel); }
void *           asContext_GetThisPointer(asIScriptContext *c, int stackLevel)                                 { return c->GetThisPointer(stackLevel); }
void *           asContext_SetUserData(asIScriptContext *c, void *data)                                        { return c->SetUserData(data); }
void *           asContext_GetUserData(asIScriptContext *c)                                                    { return c->GetUserData(); }

asIScriptEngine *asGeneric_GetEngine(asIScriptGeneric *g)                   { return g->GetEngine(); }
int              asGeneric_GetFunctionId(asIScriptGeneric *g)               { return g->GetFunctionId(); }
void *           asGeneric_GetObject(asIScriptGeneric *g)                   { return g->GetObject(); }
asBYTE           asGeneric_GetArgByte(asIScriptGeneric *g, asUINT arg)      { return g->GetArgByte(arg); }
asWORD           asGeneric_GetArgWord(asIScriptGeneric *g, asUINT arg)      { return g->GetArgWord(arg); }
asDWORD          asGeneric_GetArgDWord(asIScriptGeneric *g, asUINT arg)     { return g->GetArgDWord(arg); }
asQWORD          asGeneric_GetArgQWord(asIScriptGeneric *g, asUINT arg)     { return g->GetArgQWord(arg); }
float            asGeneric_GetArgFloat(asIScriptGeneric *g, asUINT arg)     { return g->GetArgFloat(arg); }
double           asGeneric_GetArgDouble(asIScriptGeneric *g, asUINT arg)    { return g->GetArgDouble(arg); }
void *           asGeneric_GetArgAddress(asIScriptGeneric *g, asUINT arg)   { return g->GetArgAddress(arg); }
void *           asGeneric_GetArgObject(asIScriptGeneric *g, asUINT arg)    { return g->GetArgObject(arg); }
void *           asGeneric_GetArgPointer(asIScriptGeneric *g, asUINT arg)   { return g->GetArgPointer(arg); }
int              asGeneric_GetArgTypeId(asIScriptGeneric *g, asUINT arg)    { return g->GetArgTypeId(arg); }
int              asGeneric_SetReturnByte(asIScriptGeneric *g, asBYTE val)   { return g->SetReturnByte(val); }
int              asGeneric_SetReturnWord(asIScriptGeneric *g, asWORD val)   { return g->SetReturnWord(val); }
int              asGeneric_SetReturnDWord(asIScriptGeneric *g, asDWORD val) { return g->SetReturnDWord(val); }
int              asGeneric_SetReturnQWord(asIScriptGeneric *g, asQWORD val) { return g->SetReturnQWord(val); }
int              asGeneric_SetReturnFloat(asIScriptGeneric *g, float val)   { return g->SetReturnFloat(val); }
int              asGeneric_SetReturnDouble(asIScriptGeneric *g, double val) { return g->SetReturnDouble(val); }
int              asGeneric_SetReturnAddress(asIScriptGeneric *g, void *addr) { return g->SetReturnAddress(addr); }
int              asGeneric_SetReturnObject(asIScriptGeneric *g, void *obj)  { return g->SetReturnObject(obj); }
void *           asGeneric_GetReturnPointer(asIScriptGeneric *g)            { return g->GetReturnPointer(); }

int         asStruct_AddRef(asIScriptStruct *s)                           { return s->AddRef(); }
int         asStruct_Release(asIScriptStruct *s)                          { return s->Release(); }
int         asStruct_GetStructTypeId(asIScriptStruct *s)                  { return s->GetStructTypeId(); }
int         asStruct_GetPropertyCount(asIScriptStruct *s)                 { return s->GetPropertyCount(); }
int         asStruct_GetPropertyTypeId(asIScriptStruct *s, asUINT prop)   { return s->GetPropertyTypeId(prop); }
const char *asStruct_GetPropertyName(asIScriptStruct *s, asUINT prop)     { return s->GetPropertyName(prop); }
void *      asStruct_GetPropertyPointer(asIScriptStruct *s, asUINT prop)  { return s->GetPropertyPointer(prop); }
int         asStruct_CopyFrom(asIScriptStruct *s, asIScriptStruct *other) { return s->CopyFrom(other); }

int    asArray_AddRef(asIScriptArray *a)                          { return a->AddRef(); }
int    asArray_Release(asIScriptArray *a)                         { return a->Release(); }
int    asArray_GetArrayTypeId(asIScriptArray *a)                  { return a->GetArrayTypeId(); }
int    asArray_GetElementTypeId(asIScriptArray *a)                { return a->GetElementTypeId(); }
asUINT asArray_GetElementCount(asIScriptArray *a)                 { return a->GetElementCount(); }
void * asArray_GetElementPointer(asIScriptArray *a, asUINT index) { return a->GetElementPointer(index); }
void   asArray_Resize(asIScriptArray *a, asUINT size)             { a->Resize(size); }
int    asArray_CopyFrom(asIScriptArray *a, asIScriptArray *other) { return a->CopyFrom(other); }

#endif
END_AS_NAMESPACE
