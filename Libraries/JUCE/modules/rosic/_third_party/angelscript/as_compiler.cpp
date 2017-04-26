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
// as_compiler.cpp
//
// The class that does the actual compilation of the functions
//

#include <math.h> // fmodf()

#include "as_config.h"
#include "as_compiler.h"
#include "as_tokendef.h"
#include "as_tokenizer.h"
#include "as_string_util.h"
#include "as_texts.h"
#include "as_parser.h"

BEGIN_AS_NAMESPACE

asCCompiler::asCCompiler(asCScriptEngine *engine) : byteCode(engine)
{
	builder = 0;
	script = 0;

	variables = 0;
	isProcessingDeferredParams = false;
	noCodeOutput = 0;
}

asCCompiler::~asCCompiler()
{
	while( variables )
	{
		asCVariableScope *var = variables;
		variables = variables->parent;

		DELETE(var,asCVariableScope);
	}
}

void asCCompiler::Reset(asCBuilder *builder, asCScriptCode *script, asCScriptFunction *outFunc)
{
	this->builder = builder;
	this->engine = builder->engine;
	this->script = script;
	this->outFunc = outFunc;

	hasCompileErrors = false;

	nextLabel = 0;
	breakLabels.SetLength(0);
	continueLabels.SetLength(0);

	byteCode.ClearAll();
	objVariableTypes.SetLength(0);
	objVariablePos.SetLength(0);

	globalExpression = false;
}

int asCCompiler::CompileDefaultConstructor(asCBuilder *builder, asCScriptCode *script, asCScriptFunction *outFunc)
{
	Reset(builder, script, outFunc);

	// Initialize the asCScriptStruct object then return
	byteCode.InstrPTR(BC_OBJTYPE, builder->module->RefObjectType(outFunc->objectType));
	byteCode.InstrSHORT(BC_PSF, 0);
	byteCode.Instr(BC_RDSPTR);
	byteCode.Call(BC_CALLSYS, engine->scriptTypeBehaviours.beh.construct, 2*PTR_SIZE);
	
	// Pop the object pointer from the stack
	byteCode.Ret(2*PTR_SIZE);

	byteCode.Finalize();

	// Copy byte code to the registered function
	outFunc->byteCode.SetLength(byteCode.GetSize());
	byteCode.Output(outFunc->byteCode.AddressOf());
	outFunc->stackNeeded = byteCode.largestStackUsed;
	outFunc->lineNumbers = byteCode.lineNumbers;
	outFunc->objVariablePos = objVariablePos;
	outFunc->objVariableTypes = objVariableTypes;

	return 0;
}

int asCCompiler::CompileFunction(asCBuilder *builder, asCScriptCode *script, asCScriptNode *func, asCScriptFunction *outFunc)
{
	Reset(builder, script, outFunc);

	int stackPos = 0;
	if( outFunc->objectType )
		stackPos = -PTR_SIZE; // The first parameter is the pointer to the object

	// Reserve a label for the cleanup code
	nextLabel++;

	// Add the first variable scope, which the parameters and
	// variables declared in the outermost statement block is
	// part of.
	AddVariableScope();

	//----------------------------------------------
	// Examine return type
	bool isConstructor = false;
	bool isDestructor = false;
	asCDataType returnType;
	if( func->firstChild->nodeType == snDataType )
	{
		returnType = builder->CreateDataTypeFromNode(func->firstChild, script);
		returnType = builder->ModifyDataTypeFromNode(returnType, func->firstChild->next, script, 0, 0);
	}
	else
	{
		returnType = asCDataType::CreatePrimitive(ttVoid, false);
		if( func->firstChild->tokenType == ttBitNot )
			isDestructor = true;
		else
			isConstructor = true;
	}

	//----------------------------------------------
	// Declare parameters
	// Find first parameter
	asCScriptNode *node = func->firstChild;
	while( node && node->nodeType != snParameterList )
		node = node->next;

	// Register parameters from last to first, otherwise they will be destroyed in the wrong order
	asCVariableScope vs(0);

	if( node ) node = node->firstChild;
	while( node )
	{
		// Get the parameter type
		asCDataType type = builder->CreateDataTypeFromNode(node, script);

		type = builder->ModifyDataTypeFromNode(type, node->next, script, 0, 0);

		// Is the data type allowed?
		if( type.GetSizeOnStackDWords() == 0 || type.IsReference() && type.GetSizeInMemoryBytes() == 0 )
		{
			asCString str;
			str.Format(TXT_PARAMETER_CANT_BE_s, type.Format().AddressOf());
			Error(str.AddressOf(), node);
		}

		// If the parameter has a name then declare it as variable
		node = node->next->next;
		if( node && node->nodeType == snIdentifier )
		{
			GETSTRING(name, &script->code[node->tokenPos], node->tokenLength);

			if( vs.DeclareVariable(name.AddressOf(), type, stackPos) < 0 )
				Error(TXT_PARAMETER_ALREADY_DECLARED, node);

			outFunc->AddVariable(name, type, stackPos);

			node = node->next;
		}
		else
			vs.DeclareVariable("", type, stackPos);

		// Move to next parameter
		stackPos -= type.GetSizeOnStackDWords();
	}

	int n;
	for( n = (int)vs.variables.GetLength() - 1; n >= 0; n-- )
	{
		variables->DeclareVariable(vs.variables[n]->name.AddressOf(), vs.variables[n]->type, vs.variables[n]->stackOffset);
	}

	// Is the return type allowed?
	if( (returnType.GetSizeOnStackDWords() == 0 && returnType != asCDataType::CreatePrimitive(ttVoid, false) )   || returnType.IsReference() && returnType.GetSizeInMemoryBytes() == 0 )
	{
		asCString str;
		str.Format(TXT_RETURN_CANT_BE_s, returnType.Format().AddressOf());
		Error(str.AddressOf(), node);
	}

	variables->DeclareVariable("return", returnType, stackPos);

	//--------------------------------------------
	// Compile the statement block

	// We need to parse the statement block now

	// TODO: PARSER: We can parse the statement block one statement at a time, thus save even more memory
	asCParser parser(builder);
	int r = parser.ParseStatementBlock(script, func->lastChild);
	if( r < 0 ) return -1;
	asCScriptNode *block = parser.GetScriptNode();

	bool hasReturn;
	asCByteCode bc(engine);
	LineInstr(&bc, func->lastChild->tokenPos);
	CompileStatementBlock(block, false, &hasReturn, &bc);
	LineInstr(&bc, func->lastChild->tokenPos + func->lastChild->tokenLength);

	// Make sure there is a return in all paths (if not return type is void)
	if( returnType != asCDataType::CreatePrimitive(ttVoid, false) )
	{
		if( hasReturn == false )
			Error(TXT_NOT_ALL_PATHS_RETURN, func->lastChild);
	}

	//------------------------------------------------
	// Concatenate the bytecode
	// Count total variable size
	int varSize = GetVariableOffset((int)variableAllocations.GetLength()) - 1;
	byteCode.Push(varSize);

	// Increase the reference for the object pointer, so that it is guaranteed to live during the entire call
	if( outFunc->objectType )
	{
		if( isConstructor )
		{
			// Initialize the asCScriptStruct object first
			byteCode.InstrPTR(BC_OBJTYPE, builder->module->RefObjectType(outFunc->objectType));
			byteCode.InstrSHORT(BC_PSF, 0);
			byteCode.Instr(BC_RDSPTR);
			byteCode.Call(BC_CALLSYS, engine->scriptTypeBehaviours.beh.construct, 2*PTR_SIZE);
		}

		byteCode.InstrSHORT(BC_PSF, 0);
		byteCode.Instr(BC_RDSPTR);
		byteCode.Call(BC_CALLSYS, outFunc->objectType->beh.addref, PTR_SIZE);
	}

	// Add the code for the statement block
	byteCode.AddCode(&bc);

	// Deallocate all local variables
	for( n = (int)variables->variables.GetLength() - 1; n >= 0; n-- )
	{
		sVariable *v = variables->variables[n];
		if( v->stackOffset > 0 )
		{
			// Call variables destructors
			if( v->name != "return" && v->name != "return address" )
				CompileDestructor(v->type, v->stackOffset, &byteCode);

			DeallocateVariable(v->stackOffset);
		}
	}

	// This is the label that return statements jump to
	// in order to exit the function
	byteCode.Label(0);

	// Release the object pointer again
	if( outFunc->objectType )
	{
		byteCode.InstrSHORT(BC_PSF, 0);
		byteCode.InstrPTR(BC_FREE, builder->module->RefObjectType(outFunc->objectType));
	}

	// Call destructors for function parameters
	for( n = (int)variables->variables.GetLength() - 1; n >= 0; n-- )
	{
		sVariable *v = variables->variables[n];
		if( v->stackOffset <= 0 )
		{
			// Call variable destructors here, for variables not yet destroyed
			if( v->name != "return" && v->name != "return address" )
				CompileDestructor(v->type, v->stackOffset, &byteCode);
		}

		// Do not deallocate parameters
	}

	// If there are compile errors, there is no reason to build the final code
	if( hasCompileErrors ) return -1;

	// At this point there should be no variables allocated
	assert(variableAllocations.GetLength() == freeVariables.GetLength());

	// Remove the variable scope
	RemoveVariableScope();

	byteCode.Pop(varSize);

	if( isConstructor )
		// Pop the extra object pointer as well
		byteCode.Ret(-stackPos + PTR_SIZE);
	else
		byteCode.Ret(-stackPos);

	// Tell the bytecode which variables are temporary
	for( n = 0; n < (signed)variableIsTemporary.GetLength(); n++ )
	{
		if( variableIsTemporary[n] )
			byteCode.DefineTemporaryVariable(GetVariableOffset(n));
	}

	// Finalize the bytecode
	byteCode.Finalize();

	// Compile the list of object variables for the exception handler
	for( n = 0; n < (int)variableAllocations.GetLength(); n++ )
	{
		if( variableAllocations[n].IsObject() && !variableAllocations[n].IsReference() )
		{
			objVariableTypes.PushLast(variableAllocations[n].GetObjectType());
			objVariablePos.PushLast(GetVariableOffset(n));
		}
	}

	if( hasCompileErrors ) return -1;

	return 0;
}




void asCCompiler::DefaultConstructor(asCByteCode *bc, asCDataType &type)
{
	assert(!type.IsObjectHandle());

	int func = 0;
	asSTypeBehaviour *beh = type.GetBehaviour();
	if( beh ) func = beh->construct;

	if( type.IsScriptArray() )
	{
		// The script struct constructor needs to know what type it is
		asCObjectType *objType = type.GetObjectType();
		bc->InstrPTR(BC_OBJTYPE, builder->module->RefObjectType(objType));

		bc->Alloc(BC_ALLOC, builder->module->RefObjectType(objType), func, 2*PTR_SIZE);
	}
	else
	{
		bc->Alloc(BC_ALLOC, builder->module->RefObjectType(type.GetObjectType()), func, PTR_SIZE);
	}
}

void asCCompiler::CompileConstructor(asCDataType &type, int offset, asCByteCode *bc)
{
	// Call constructor for the data type
	if( type.IsObject() && !type.IsObjectHandle() )
	{
		bc->InstrSHORT(BC_PSF, (short)offset);
		DefaultConstructor(bc, type);
	}
}

void asCCompiler::CompileDestructor(asCDataType &type, int offset, asCByteCode *bc)
{
	if( !type.IsReference() )
	{
		// Call destructor for the data type
		if( type.IsObject() )
		{
			// Free the memory
			bc->InstrSHORT(BC_PSF, (short)offset);
			bc->InstrPTR(BC_FREE, builder->module->RefObjectType(type.GetObjectType()));
		}
	}
}

void asCCompiler::LineInstr(asCByteCode *bc, size_t pos)
{
	int r, c;
	script->ConvertPosToRowCol(pos, &r, &c);
	bc->Line(r, c);
}

void asCCompiler::CompileStatementBlock(asCScriptNode *block, bool ownVariableScope, bool *hasReturn, asCByteCode *bc)
{
	*hasReturn = false;
	bool isFinished = false;
	bool hasWarned = false;

	if( ownVariableScope )
		AddVariableScope();

	asCScriptNode *node = block->firstChild;
	while( node )
	{
		if( !hasWarned && (*hasReturn || isFinished) )
		{
			hasWarned = true;
			Warning(TXT_UNREACHABLE_CODE, node);
		}

		if( node->nodeType == snBreak || node->nodeType == snContinue )
			isFinished = true;

		asCByteCode statement(engine);
		if( node->nodeType == snDeclaration )
			CompileDeclaration(node, &statement);
		else
			CompileStatement(node, hasReturn, &statement);

		LineInstr(bc, node->tokenPos);
		bc->AddCode(&statement);

		if( !hasCompileErrors )
			assert( tempVariables.GetLength() == 0 );

		node = node->next;
	}

	if( ownVariableScope )
	{

		// Deallocate variables in this block, in reverse order
		for( int n = (int)variables->variables.GetLength() - 1; n >= 0; n-- )
		{
			sVariable *v = variables->variables[n];

			// Call variable destructors here, for variables not yet destroyed
			// If the block is terminated with a break, continue, or
			// return the variables are already destroyed
			if( !isFinished && !*hasReturn )
				CompileDestructor(v->type, v->stackOffset, bc);

			// Don't deallocate function parameters
			if( v->stackOffset > 0 )
				DeallocateVariable(v->stackOffset);
		}

		RemoveVariableScope();
	}
}

int asCCompiler::CompileGlobalVariable(asCBuilder *builder, asCScriptCode *script, asCScriptNode *node, sGlobalVariableDescription *gvar)
{
	Reset(builder, script, 0);
	globalExpression = true;

	// Add a variable scope (even though variables can't be declared)
	AddVariableScope();

	asSExprContext ctx(engine);

	gvar->isPureConstant = false;

	// Parse the initialization nodes
	asCParser parser(builder);
	if( node )
	{
		int r = parser.ParseGlobalVarInit(script, node);
		if( r < 0 )
			return r;

		node = parser.GetScriptNode();
	}

	// Compile the expression
	if( node && node->nodeType == snArgList )
	{
		// Make sure that it is a registered type, and that it isn't a pointer
		if( gvar->datatype.GetObjectType() == 0 || gvar->datatype.IsObjectHandle() )
		{
			Error(TXT_MUST_BE_OBJECT, node);
		}
		else
		{
			// Compile the arguments
			asCArray<asSExprContext *> args;
			CompileArgumentList(node, args, &gvar->datatype);

			// Find all constructors
			asCArray<int> funcs;
			asSTypeBehaviour *beh = gvar->datatype.GetBehaviour();
			if( beh )
				funcs = beh->constructors;

			asCString str = gvar->datatype.Format();
			MatchFunctions(funcs, args, node, str.AddressOf());

			if( funcs.GetLength() == 1 )
			{
				// TODO: This reference is open while evaluating the arguments. We must fix this
				ctx.bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(gvar->index));

				PrepareFunctionCall(funcs[0], &ctx.bc, args);
				MoveArgsToStack(funcs[0], &ctx.bc, args, false);

				PerformFunctionCall(funcs[0], &ctx, true, &args, gvar->datatype.GetObjectType());
			}

			// Cleanup
			for( asUINT n = 0; n < args.GetLength(); n++ )
				if( args[n] ) 
				{
					DELETE(args[n],asSExprContext);
				}
		}
	}
	else if( node && node->nodeType == snInitList )
	{
		asCTypeInfo ti;
		ti.Set(gvar->datatype);
		ti.isVariable = false;
		ti.isTemporary = false;
		ti.stackOffset = (short)gvar->index;

		CompileInitList(&ti, node, &ctx.bc);

		node = node->next;
	}
	else
	{
		// Call constructor for all data types
		if( gvar->datatype.IsObject() && !gvar->datatype.IsObjectHandle() )
		{
			ctx.bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(gvar->index));
			DefaultConstructor(&ctx.bc, gvar->datatype);
		}

		if( node )
		{
			asSExprContext expr(engine);
			CompileAssignment(node, &expr);

			if( gvar->datatype.IsPrimitive() )
			{
				if( gvar->datatype.IsReadOnly() && expr.type.isConstant )
				{
					ImplicitConversion(&expr, gvar->datatype, node, false);

					gvar->isPureConstant = true;
					gvar->constantValue = expr.type.qwordValue;
				}

				asSExprContext lctx(engine);
				lctx.type.Set(gvar->datatype);
				lctx.type.dataType.MakeReference(true);
				lctx.type.dataType.MakeReadOnly(false);
				lctx.bc.InstrWORD(BC_LDG, (asWORD)builder->module->GetGlobalVarIndex(gvar->index));

				DoAssignment(&ctx, &lctx, &expr, node, node, ttAssignment, node);
			}
			else
			{
				asCTypeInfo ltype;
				ltype.Set(gvar->datatype);
				ltype.dataType.MakeReference(true);
				ltype.dataType.MakeReadOnly(false);
				ltype.stackOffset = -1;

				if( gvar->datatype.IsObjectHandle() )
					ltype.isExplicitHandle = true;

				// If left expression resolves into a registered type
				// check if the assignment operator is overloaded, and check
				// the type of the right hand expression. If none is found
				// the default action is a direct copy if it is the same type
				// and a simple assignment.
				asSTypeBehaviour *beh = 0;
				if( !ltype.isExplicitHandle )
					beh = ltype.dataType.GetBehaviour();
				bool assigned = false;
				if( beh )
				{
					// Find the matching overloaded operators
					int op = ttAssignment;
					asCArray<int> ops;
					asUINT n;
					for( n = 0; n < beh->operators.GetLength(); n += 2 )
					{
						if( op == beh->operators[n] )
							ops.PushLast(beh->operators[n+1]);
					}

					asCArray<int> match;
					MatchArgument(ops, match, &expr.type, 0);

					if( match.GetLength() > 0 )
						assigned = true;

					if( match.GetLength() == 1 )
					{
						// If it is an array, both sides must have the same subtype
						if( ltype.dataType.IsArrayType() )
							if( !ltype.dataType.IsEqualExceptRefAndConst(expr.type.dataType) )
								Error(TXT_BOTH_MUST_BE_SAME, node);

						asCScriptFunction *descr = engine->scriptFunctions[match[0]];

						// Add code for arguments
						MergeExprContexts(&ctx, &expr);

						PrepareArgument(&descr->parameterTypes[0], &expr, node, true, descr->inOutFlags[0]);
						MergeExprContexts(&ctx, &expr);


						asCArray<asSExprContext*> args;
						args.PushLast(&expr);
						MoveArgsToStack(match[0], &ctx.bc, args, false);

						// Add the code for the object
						ctx.bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(gvar->index));
						ctx.bc.Instr(BC_RDSPTR);

						PerformFunctionCall(match[0], &ctx, false, &args);

						ctx.bc.Pop(ctx.type.dataType.GetSizeOnStackDWords());

						ProcessDeferredParams(&ctx);
					}
					else if( match.GetLength() > 1 )
					{
						Error(TXT_MORE_THAN_ONE_MATCHING_OP, node);
					}
				}

				if( !assigned )
				{
					PrepareForAssignment(&ltype.dataType, &expr, node);

					// If the expression is constant and the variable also is constant
					// then mark the variable as pure constant. This will allow the compiler
					// to optimize expressions with this variable.
					if( gvar->datatype.IsReadOnly() && expr.type.isConstant )
					{
						gvar->isPureConstant = true;
						gvar->constantValue = expr.type.qwordValue;
					}

					// Add expression code to bytecode
					MergeExprContexts(&ctx, &expr);

					// Add byte code for storing value of expression in variable
					ctx.bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(gvar->index));

					PerformAssignment(&ltype, &expr.type, &ctx.bc, node->prev);

					// Release temporary variables used by expression
					ReleaseTemporaryVariable(expr.type, &ctx.bc);

					ctx.bc.Pop(expr.type.dataType.GetSizeOnStackDWords());
				}
			}
		}
	}

	// Concatenate the bytecode
	int varSize = GetVariableOffset((int)variableAllocations.GetLength()) - 1;
	byteCode.Push(varSize);

	byteCode.AddCode(&ctx.bc);

	// Deallocate variables in this block, in reverse order
	for( int n = (int)variables->variables.GetLength() - 1; n >= 0; --n )
	{
		sVariable *v = variables->variables[n];

		// Call variable destructors here, for variables not yet destroyed
		CompileDestructor(v->type, v->stackOffset, &byteCode);

		DeallocateVariable(v->stackOffset);
	}

	if( hasCompileErrors ) return -1;

	// At this point there should be no variables allocated
	assert(variableAllocations.GetLength() == freeVariables.GetLength());

	// Remove the variable scope again
	RemoveVariableScope();

	byteCode.Pop(varSize);

	return 0;
}

void asCCompiler::PrepareArgument(asCDataType *paramType, asSExprContext *ctx, asCScriptNode *node, bool isFunction, int refType, asCArray<int> *reservedVars)
{
	asCDataType param = *paramType;
	if( paramType->GetTokenType() == ttQuestion )
	{
		// Since the function is expecting a var type ?, then we don't want to convert the argument to anything else
		param = ctx->type.dataType;
		param.MakeHandle(ctx->type.isExplicitHandle);
		param.MakeReference(paramType->IsReference());
		param.MakeReadOnly(paramType->IsReadOnly());
	}
	else
		param = *paramType;

	asCDataType dt = param;

	// Need to protect arguments by reference
	if( isFunction && dt.IsReference() )
	{
		if( paramType->GetTokenType() == ttQuestion )
		{
			asCByteCode tmpBC(engine);

			// Place the type id on the stack as a hidden parameter
			tmpBC.InstrDWORD(BC_TYPEID, engine->GetTypeIdFromDataType(param));

			// Insert the code before the expression code
			tmpBC.AddCode(&ctx->bc);
			ctx->bc.AddCode(&tmpBC);
		}

		// Allocate a temporary variable of the same type as the argument
		dt.MakeReference(false);
		dt.MakeReadOnly(false);

		int offset;
		if( refType == 1 ) // &in
		{
			// If the reference is const, then it is not necessary to make a copy if the value already is a variable
			// Even if the same variable is passed in another argument as non-const then there is no problem
			if( dt.IsPrimitive() )
			{
				IsVariableInitialized(&ctx->type, node);

				if( ctx->type.dataType.IsReference() ) ConvertToVariable(ctx);
				ImplicitConversion(ctx, dt, node, false, true, reservedVars);

				if( !(param.IsReadOnly() && ctx->type.isVariable) )
					ConvertToTempVariable(ctx);

				PushVariableOnStack(ctx, true);
				ctx->type.dataType.MakeReadOnly(param.IsReadOnly());
			}
			else
			{
				IsVariableInitialized(&ctx->type, node);

				ImplicitConversion(ctx, param, node, false, true, reservedVars);

				// If the argument already is a temporary
				// variable we don't need to allocate another

				// If the parameter is read-only and the object already is a local
				// variable then it is not necessary to make a copy either
				if( !ctx->type.isTemporary && !(param.IsReadOnly() && ctx->type.isVariable))
				{
					// Make sure the variable is not used in the expression
					asCArray<int> vars;
					ctx->bc.GetVarsUsed(vars);
					if( reservedVars ) vars.Concatenate(*reservedVars);
					offset = AllocateVariableNotIn(dt, true, &vars);

					// Allocate and construct the temporary object
					asCByteCode tmpBC(engine);
					CompileConstructor(dt, offset, &tmpBC);

					// Insert the code before the expression code
					tmpBC.AddCode(&ctx->bc);
					ctx->bc.AddCode(&tmpBC);

					// Assign the evaluated expression to the temporary variable
					PrepareForAssignment(&dt, ctx, node);

					dt.MakeReference(true);
					asCTypeInfo type;
					type.Set(dt);
					type.isTemporary = true;
					type.stackOffset = (short)offset;

					if( dt.IsObjectHandle() )
						type.isExplicitHandle = true;

					ctx->bc.InstrSHORT(BC_PSF, (short)offset);

					PerformAssignment(&type, &ctx->type, &ctx->bc, node);

					ctx->bc.Pop(ctx->type.dataType.GetSizeOnStackDWords());

					ReleaseTemporaryVariable(ctx->type, &ctx->bc);

					ctx->type = type;

					ctx->bc.InstrSHORT(BC_PSF, (short)offset);
					if( dt.IsObject() && !dt.IsObjectHandle() )
						ctx->bc.Instr(BC_RDSPTR);

					if( paramType->IsReadOnly() )
						ctx->type.dataType.MakeReadOnly(true);
				}
			}
		}
		else if( refType == 2 ) // &out
		{
			// Make sure the variable is not used in the expression
			asCArray<int> vars;
			ctx->bc.GetVarsUsed(vars);
			if( reservedVars ) vars.Concatenate(*reservedVars);
			offset = AllocateVariableNotIn(dt, true, &vars);

			if( dt.IsPrimitive() )
			{
				ctx->type.SetVariable(dt, offset, true);
				PushVariableOnStack(ctx, true);
			}
			else
			{
				// Allocate and construct the temporary object
				asCByteCode tmpBC(engine);
				CompileConstructor(dt, offset, &tmpBC);

				// Insert the code before the expression code
				tmpBC.AddCode(&ctx->bc);
				ctx->bc.AddCode(&tmpBC);

				dt.MakeReference((!dt.IsObject() || dt.IsObjectHandle()));
				asCTypeInfo type;
				type.Set(dt);
				type.isTemporary = true;
				type.stackOffset = (short)offset;

				ctx->type = type;

				ctx->bc.InstrSHORT(BC_PSF, (short)offset);
				if( dt.IsObject() && !dt.IsObjectHandle() )
					ctx->bc.Instr(BC_RDSPTR);
			}

			// After the function returns the temporary variable will
			// be assigned to the expression, if it is a valid lvalue
		}
		else if( refType == 3 ) // &inout
		{
			// Only objects that support object handles 
			// can be guaranteed to be safe. Local variables are
			// already safe, so there is no need to add an extra
			// references
			if( !engine->allowUnsafeReferences &&
				!ctx->type.isVariable &&
				ctx->type.dataType.IsObject() &&
				!ctx->type.dataType.IsObjectHandle() &&
				ctx->type.dataType.GetBehaviour()->addref &&
				ctx->type.dataType.GetBehaviour()->release )
			{
				// Store a handle to the object as local variable
				asSExprContext tmp(engine);
				asCDataType dt = ctx->type.dataType;
				dt.MakeHandle(true);

				asCArray<int> vars;
				ctx->bc.GetVarsUsed(vars);
				if( reservedVars ) vars.Concatenate(*reservedVars);
				offset = AllocateVariableNotIn(dt, true, &vars);

				// Copy the handle
				ctx->bc.InstrWORD(BC_PSF, (asWORD)offset);
				ctx->bc.InstrPTR(BC_REFCPY, builder->module->RefObjectType(ctx->type.dataType.GetObjectType()));
				ctx->bc.Pop(PTR_SIZE);
				ctx->bc.InstrWORD(BC_PSF, (asWORD)offset);

				dt.MakeHandle(false);
				dt.MakeReference(true);
				ctx->type.SetVariable(dt, offset, true);
			}

			// Make sure the reference to the value is on the stack
			if( ctx->type.dataType.IsObject() && ctx->type.dataType.IsReference() )
				Dereference(ctx, true);
			else if( ctx->type.isVariable )
				ctx->bc.InstrSHORT(BC_PSF, ctx->type.stackOffset);
			else if( ctx->type.dataType.IsPrimitive() )
				ctx->bc.Instr(BC_PshRPtr);
		}
	}
	else
	{
		if( dt.IsPrimitive() )
		{
			IsVariableInitialized(&ctx->type, node);

			if( ctx->type.dataType.IsReference() ) ConvertToVariable(ctx);

			// Implicitly convert primitives to the parameter type
			ImplicitConversion(ctx, dt, node, false, true, reservedVars);

			if( ctx->type.isVariable )
			{
				PushVariableOnStack(ctx, dt.IsReference());
			}
			else if( ctx->type.isConstant )
			{
				ConvertToVariable(ctx);
				PushVariableOnStack(ctx, dt.IsReference());
			}
		}
		else
		{
			IsVariableInitialized(&ctx->type, node);

			// Implicitly convert primitives to the parameter type
			ImplicitConversion(ctx, dt, node, false, true, reservedVars);

			if( dt.IsObjectHandle() )
				ctx->type.isExplicitHandle = true;

			if( dt.IsObject() )
			{
				if( !dt.IsReference() )
				{
					// Objects passed by value must be placed in temporary variables
					// so that they are guaranteed to not be referenced anywhere else
					PrepareTemporaryObject(node, ctx);

					// The implicit conversion shouldn't convert the object to
					// non-reference yet. It will be dereferenced just before the call.
					// Otherwise the object might be missed by the exception handler.
					dt.MakeReference(true);
				}
				else
				{
					// An object passed by reference should place the pointer to
					// the object on the stack.
					dt.MakeReference(false);
				}
			}
		}
	}

	// Don't put any pointer on the stack yet
	if( param.IsReference() || param.IsObject() )
	{
		// &inout parameter may leave the reference on the stack already
		if( refType != 3 )
		{
			ctx->bc.Pop(PTR_SIZE);
			ctx->bc.InstrSHORT(BC_VAR, ctx->type.stackOffset);
		}

		ProcessDeferredParams(ctx);
	}
}

void asCCompiler::PrepareFunctionCall(int funcID, asCByteCode *bc, asCArray<asSExprContext *> &args)
{
	// When a match has been found, compile the final byte code using correct parameter types
	asCScriptFunction *descr = builder->GetFunctionDescription(funcID);

	// Add code for arguments
	asSExprContext e(engine);
	int n;
	for( n = (int)args.GetLength()-1; n >= 0; n-- )
	{
		// Make sure PrepareArgument doesn't use any variable that is already 
		// being used by any of the following argument expressions
		asCArray<int> reservedVars;
		for( int m = n-1; m >= 0; m-- )
			args[m]->bc.GetVarsUsed(reservedVars);

		PrepareArgument2(&e, args[n], &descr->parameterTypes[n], true, descr->inOutFlags[n], &reservedVars);
	}

	bc->AddCode(&e.bc);
}

void asCCompiler::MoveArgsToStack(int funcID, asCByteCode *bc, asCArray<asSExprContext *> &args, bool addOneToOffset)
{
	asCScriptFunction *descr = builder->GetFunctionDescription(funcID);

	int offset = 0;
	if( addOneToOffset )
		offset += PTR_SIZE;

	// Move the objects that are sent by value to the stack just before the call
	for( asUINT n = 0; n < descr->parameterTypes.GetLength(); n++ )
	{
		if( descr->parameterTypes[n].IsReference() )
		{
			if( descr->parameterTypes[n].IsObject() && !descr->parameterTypes[n].IsObjectHandle() )
			{
				if( descr->inOutFlags[n] != 3 )
					bc->InstrWORD(BC_GETOBJREF, (asWORD)offset);
				if( args[n]->type.dataType.IsObjectHandle() )
					bc->InstrWORD(BC_ChkNullS, (asWORD)offset);
			}
			else if( descr->inOutFlags[n] != 3 )
				bc->InstrWORD(BC_GETREF, (asWORD)offset);
		}
		else if( descr->parameterTypes[n].IsObject() )
		{
			bc->InstrWORD(BC_GETOBJ, (asWORD)offset);

			// The temporary variable must not be freed as it will no longer hold an object
			DeallocateVariable(args[n]->type.stackOffset);
			args[n]->type.isTemporary = false;
		}

		offset += descr->parameterTypes[n].GetSizeOnStackDWords();
	}
}

void asCCompiler::CompileArgumentList(asCScriptNode *node, asCArray<asSExprContext*> &args, asCDataType *type)
{
	assert(node->nodeType == snArgList);

	// Count arguments
	asCScriptNode *arg = node->firstChild;
	int argCount = 0;
	while( arg )
	{
		argCount++;
		arg = arg->next;
	}

	// The script array and script any needs to receive their object type
	if( type && type->IsScriptArray() )
	{
		argCount += 1;
	}

	// Prepare the arrays
	args.SetLength(argCount);
	int n;
	for( n = 0; n < argCount; n++ )
		args[n] = 0;

	n = argCount-1;
	if( type && type->IsScriptArray() )
	{
		args[n] = NEW(asSExprContext)(engine);
		args[n]->bc.InstrPTR(BC_OBJTYPE, builder->module->RefObjectType(type->GetObjectType()));
#ifndef AS_64BIT_PTR
		args[n]->type.Set(asCDataType::CreatePrimitive(ttInt, false));
#else
		args[n]->type.Set(asCDataType::CreatePrimitive(ttInt64, false));
#endif
		n--;
	}

	// Compile the arguments in reverse order (as they will be pushed on the stack)
	arg = node->lastChild;
	while( arg )
	{
		asSExprContext expr(engine);
		CompileAssignment(arg, &expr);

		args[n] = NEW(asSExprContext)(engine);
		MergeExprContexts(args[n], &expr);
		args[n]->type = expr.type;
		args[n]->exprNode = arg;

		n--;
		arg = arg->prev;
	}
}

void asCCompiler::MatchFunctions(asCArray<int> &funcs, asCArray<asSExprContext*> &args, asCScriptNode *node, const char *name, bool isConstMethod, bool silent, bool allowObjectConstruct)
{
	asUINT n;
	if( funcs.GetLength() > 0 )
	{
		// Check the number of parameters in the found functions
		for( n = 0; n < funcs.GetLength(); ++n )
		{
			asCScriptFunction *desc = builder->GetFunctionDescription(funcs[n]);

			if( desc->parameterTypes.GetLength() != args.GetLength() )
			{
				// remove it from the list
				if( n == funcs.GetLength()-1 )
					funcs.PopLast();
				else
					funcs[n] = funcs.PopLast();
				n--;
			}
		}

		// Match functions with the parameters, and discard those that do not match
		asCArray<int> matchingFuncs = funcs;

		for( n = 0; n < args.GetLength(); ++n )
		{
			asCArray<int> tempFuncs;
			MatchArgument(funcs, tempFuncs, &args[n]->type, n, allowObjectConstruct);

			// Intersect the found functions with the list of matching functions
			for( asUINT f = 0; f < matchingFuncs.GetLength(); f++ )
			{
				asUINT c;
				for( c = 0; c < tempFuncs.GetLength(); c++ )
				{
					if( matchingFuncs[f] == tempFuncs[c] )
						break;
				}

				// Was the function a match?
				if( c == tempFuncs.GetLength() )
				{
					// No, remove it from the list
					if( f == matchingFuncs.GetLength()-1 )
						matchingFuncs.PopLast();
					else
						matchingFuncs[f] = matchingFuncs.PopLast();
					f--;
				}
			}
		}

		funcs = matchingFuncs;
	}

	if( !isConstMethod )
		FilterConst(funcs);

	if( funcs.GetLength() != 1 && !silent )
	{
		// Build a readable string of the function with parameter types
		asCString str = name;
		str += "(";
		if( args.GetLength() )
			str += args[0]->type.dataType.Format();
		for( n = 1; n < args.GetLength(); n++ )
			str += ", " + args[n]->type.dataType.Format();
		str += ")";

		if( isConstMethod )
			str += " const";

		if( funcs.GetLength() == 0 )
			str.Format(TXT_NO_MATCHING_SIGNATURES_TO_s, str.AddressOf());
		else
			str.Format(TXT_MULTIPLE_MATCHING_SIGNATURES_TO_s, str.AddressOf());

		Error(str.AddressOf(), node);
	}
}

void asCCompiler::CompileDeclaration(asCScriptNode *decl, asCByteCode *bc)
{
	// Get the data type
	asCDataType type = builder->CreateDataTypeFromNode(decl->firstChild, script);

	// Declare all variables in this declaration
	asCScriptNode *node = decl->firstChild->next;
	while( node )
	{
		// Is the type allowed?
		if( type.GetSizeOnStackDWords() == 0 ||
			(type.IsObject() && !type.IsObjectHandle() && type.GetSizeInMemoryBytes() == 0) )
		{
			asCString str;
			// TODO: Change to "'type' cannot be declared as variable"
			str.Format(TXT_DATA_TYPE_CANT_BE_s, type.Format().AddressOf());
			Error(str.AddressOf(), node);
			
			// Use int instead to avoid further problems
			type = asCDataType::CreatePrimitive(ttInt, false);
		}

		// Get the name of the identifier
		GETSTRING(name, &script->code[node->tokenPos], node->tokenLength);

		// Verify that the name isn't used by a dynamic data type
		if( engine->GetObjectType(name.AddressOf()) != 0 )
		{
			asCString str;
			str.Format(TXT_ILLEGAL_VARIABLE_NAME_s, name.AddressOf());
			Error(str.AddressOf(), node);
		}

		int offset = AllocateVariable(type, false);
		if( variables->DeclareVariable(name.AddressOf(), type, offset) < 0 )
		{
			asCString str;
			str.Format(TXT_s_ALREADY_DECLARED, name.AddressOf());
			Error(str.AddressOf(), node);
		}

		outFunc->AddVariable(name, type, offset);

		node = node->next;
		if( node && node->nodeType == snArgList )
		{
			// Make sure that it is a registered type, and that is isn't a pointer
			if( type.GetObjectType() == 0 || type.IsObjectHandle() )
			{
				Error(TXT_MUST_BE_OBJECT, node);
			}
			else
			{
				// Compile the arguments
				asCArray<asSExprContext *> args;

				CompileArgumentList(node, args, &type);

				// Find all constructors
				asCArray<int> funcs;
				asSTypeBehaviour *beh = type.GetBehaviour();
				if( beh )
					funcs = beh->constructors;

				asCString str = type.Format();
				MatchFunctions(funcs, args, node, str.AddressOf());

				if( funcs.GetLength() == 1 )
				{
					asSExprContext ctx(engine);

					sVariable *v = variables->GetVariable(name.AddressOf());
					ctx.bc.InstrSHORT(BC_VAR, (short)v->stackOffset);

					PrepareFunctionCall(funcs[0], &ctx.bc, args);
					MoveArgsToStack(funcs[0], &ctx.bc, args, false);

					int offset = 0;
					asCScriptFunction *descr = builder->GetFunctionDescription(funcs[0]);
					for( asUINT n = 0; n < args.GetLength(); n++ )
						offset += descr->parameterTypes[n].GetSizeOnStackDWords();

					ctx.bc.InstrWORD(BC_GETREF, (asWORD)offset);

					PerformFunctionCall(funcs[0], &ctx, true, &args, type.GetObjectType());

					bc->AddCode(&ctx.bc);
				}

				// Cleanup
				for( asUINT n = 0; n < args.GetLength(); n++ )
					if( args[n] ) 
					{
						DELETE(args[n],asSExprContext);
					}
			}

			node = node->next;
		}
		else if( node && node->nodeType == snInitList )
		{
			sVariable *v = variables->GetVariable(name.AddressOf());

			asCTypeInfo ti;
			ti.Set(type);
			ti.isVariable = true;
			ti.isTemporary = false;
			ti.stackOffset = (short)v->stackOffset;

			CompileInitList(&ti, node, bc);

			node = node->next;
		}
		else
		{
			asSExprContext ctx(engine);

			// Call the default constructor here
			CompileConstructor(type, offset, &ctx.bc);

			// Is the variable initialized?
			if( node && node->nodeType == snAssignment )
			{
				// Compile the expression
				asSExprContext expr(engine);
				CompileAssignment(node, &expr);

				if( type.IsPrimitive() )
				{
					if( type.IsReadOnly() && expr.type.isConstant )
					{
						ImplicitConversion(&expr, type, node, false);

						sVariable *v = variables->GetVariable(name.AddressOf());
						v->isPureConstant = true;
						v->constantValue = expr.type.qwordValue;
					}

					asSExprContext lctx(engine);
					lctx.type.SetVariable(type, offset, false);
					lctx.type.dataType.MakeReadOnly(false);

					DoAssignment(&ctx, &lctx, &expr, node, node, ttAssignment, node);
				}
				else
				{
					// TODO: We can use a copy constructor here

					asCTypeInfo ltype;
					ltype.Set(type);
					ltype.dataType.MakeReference(true);
					// Allow initialization of constant variables
					ltype.dataType.MakeReadOnly(false);

					if( type.IsObjectHandle() )
						ltype.isExplicitHandle = true;


					// If left expression resolves into a registered type
					// check if the assignment operator is overloaded, and check
					// the type of the right hand expression. If none is found
					// the default action is a direct copy if it is the same type
					// and a simple assignment.
					asSTypeBehaviour *beh = 0;
					if( !ltype.isExplicitHandle )
						beh = ltype.dataType.GetBehaviour();
					bool assigned = false;
					if( beh )
					{
						// Find the matching overloaded operators
						int op = ttAssignment;
						asCArray<int> ops;
						asUINT n;
						for( n = 0; n < beh->operators.GetLength(); n += 2 )
						{
							if( op == beh->operators[n] )
								ops.PushLast(beh->operators[n+1]);
						}

						asCArray<int> match;
						MatchArgument(ops, match, &expr.type, 0);

						if( match.GetLength() > 0 )
							assigned = true;

						if( match.GetLength() == 1 )
						{
							// If it is an array, both sides must have the same subtype
							if( ltype.dataType.IsArrayType() )
								if( !ltype.dataType.IsEqualExceptRefAndConst(expr.type.dataType) )
									Error(TXT_BOTH_MUST_BE_SAME, node);

							asCScriptFunction *descr = engine->scriptFunctions[match[0]];

							// Add code for arguments
							MergeExprContexts(&ctx, &expr);

							PrepareArgument(&descr->parameterTypes[0], &expr, node, true, descr->inOutFlags[0]);
							MergeExprContexts(&ctx, &expr);

							asCArray<asSExprContext*> args;
							args.PushLast(&expr);
							MoveArgsToStack(match[0], &ctx.bc, args, false);

							// Add the code for the object
							sVariable *v = variables->GetVariable(name.AddressOf());
							ltype.stackOffset = (short)v->stackOffset;
							ctx.bc.InstrSHORT(BC_PSF, (short)v->stackOffset);
							ctx.bc.Instr(BC_RDSPTR);

							PerformFunctionCall(match[0], &ctx, false, &args);

							ctx.bc.Pop(ctx.type.dataType.GetSizeOnStackDWords());

							ProcessDeferredParams(&ctx);
						}
						else if( match.GetLength() > 1 )
						{
							Error(TXT_MORE_THAN_ONE_MATCHING_OP, node);
						}
					}

					if( !assigned )
					{
						PrepareForAssignment(&ltype.dataType, &expr, node);

						// If the expression is constant and the variable also is constant
						// then mark the variable as pure constant. This will allow the compiler
						// to optimize expressions with this variable.
						sVariable *v = variables->GetVariable(name.AddressOf());
						if( v->type.IsReadOnly() && expr.type.isConstant )
						{
							v->isPureConstant = true;
							v->constantValue = expr.type.qwordValue;
						}

						// Add expression code to bytecode
						MergeExprContexts(&ctx, &expr);

						// Add byte code for storing value of expression in variable
						ctx.bc.InstrSHORT(BC_PSF, (short)v->stackOffset);
						ltype.stackOffset = (short)v->stackOffset;

						PerformAssignment(&ltype, &expr.type, &ctx.bc, node->prev);

						// Release temporary variables used by expression
						ReleaseTemporaryVariable(expr.type, &ctx.bc);

						ctx.bc.Pop(expr.type.dataType.GetSizeOnStackDWords());

						ProcessDeferredParams(&ctx);
					}
				}

				node = node->next;
			}

			bc->AddCode(&ctx.bc);

			// TODO: Can't this leave deferred output params without being compiled?
		}
	}
}

void asCCompiler::CompileInitList(asCTypeInfo *var, asCScriptNode *node, asCByteCode *bc)
{
	if( var->dataType.IsArrayType() && !var->dataType.IsObjectHandle() )
	{
		// Count the number of elements and initialize the array with the correct size
		int countElements = 0;
		asCScriptNode *el = node->firstChild;
		while( el )
		{
			countElements++;
			el = el->next;
		}

		// Construct the array with the size elements

		// Find the constructor that takes an uint
		asCArray<int> funcs;
		funcs = var->dataType.GetBehaviour()->constructors;

		asCArray<asSExprContext *> args;
		asSExprContext arg1(engine), arg2(engine);
		arg1.bc.InstrDWORD(BC_PshC4, countElements);
		arg1.type.Set(asCDataType::CreatePrimitive(ttUInt, false));
		args.PushLast(&arg1);

		if( var->dataType.IsScriptArray() )
		{
			// Script arrays need the type id as well
			arg2.bc.InstrPTR(BC_OBJTYPE, builder->module->RefObjectType(var->dataType.GetObjectType()));
#ifndef AS_64BIT_PTR
			arg2.type.Set(asCDataType::CreatePrimitive(ttInt, false));
#else
			arg2.type.Set(asCDataType::CreatePrimitive(ttInt64, false));
#endif
			args.PushLast(&arg2);
		}


		asCString str = var->dataType.Format();
		MatchFunctions(funcs, args, node, str.AddressOf());

		if( funcs.GetLength() == 1 )
		{
			asSExprContext ctx(engine);

			if( var->isVariable )
				ctx.bc.InstrSHORT(BC_PSF, var->stackOffset);
			else
				ctx.bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(var->stackOffset));

			PrepareFunctionCall(funcs[0], &ctx.bc, args);
			MoveArgsToStack(funcs[0], &ctx.bc, args, false);

			PerformFunctionCall(funcs[0], &ctx, true, &args, var->dataType.GetObjectType());

			bc->AddCode(&ctx.bc);
		}
		else
			return;

		// Find the indexing operator that is not read-only that will be used for all elements
		asCDataType retType;
		if( var->dataType.IsScriptArray() )
		{
			retType = asCDataType::CreatePrimitive(ttInt, false);
			retType.MakeReference(true);
		}
		else
		{
			retType = var->dataType.GetSubType();
			retType.MakeReference(true);
			retType.MakeReadOnly(false);
		}
		int funcId = 0;
		asSTypeBehaviour *beh = var->dataType.GetBehaviour();
		for( asUINT n = 0; n < beh->operators.GetLength(); n += 2 )
		{
			if( ttOpenBracket == beh->operators[n] )
			{
				asCScriptFunction *desc = builder->GetFunctionDescription(beh->operators[n+1]);
				if( !desc->isReadOnly &&
					 desc->parameterTypes.GetLength() == 1 &&
					 (desc->parameterTypes[0] == asCDataType::CreatePrimitive(ttUInt, false) ||
					  desc->parameterTypes[0] == asCDataType::CreatePrimitive(ttInt,  false)) &&
					 desc->returnType == retType )
				{
					funcId = beh->operators[n+1];
					break;
				}
			}
		}

		if( funcId == 0 )
		{
			Error(TXT_NO_APPROPRIATE_INDEX_OPERATOR, node);
			return;
		}

		// Add a reference for the configuration group
		builder->module->RefConfigGroupForFunction(funcId);

		asUINT index = 0;
		el = node->firstChild;
		while( el )
		{
			if( el->nodeType == snAssignment || el->nodeType == snInitList )
			{
				asSExprContext lctx(engine);
				asSExprContext rctx(engine);

				if( el->nodeType == snAssignment )
				{
					// Compile the assignment expression
					CompileAssignment(el, &rctx);
				}
				else if( el->nodeType == snInitList )
				{
					int offset = AllocateVariable(var->dataType.GetSubType(), true);

					rctx.type.Set(var->dataType.GetSubType());
					rctx.type.isVariable = true;
					rctx.type.isTemporary = true;
					rctx.type.stackOffset = (short)offset;

					CompileInitList(&rctx.type, el, &rctx.bc);

					// Put the object on the stack
					rctx.bc.InstrSHORT(BC_PSF, rctx.type.stackOffset);
				}

				// Compile the lvalue
				lctx.bc.InstrDWORD(BC_PshC4, index);
				if( var->isVariable )
					lctx.bc.InstrSHORT(BC_PSF, var->stackOffset);
				else
					lctx.bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(var->stackOffset));
				lctx.bc.Instr(BC_RDSPTR);
				lctx.bc.Call(BC_CALLSYS, funcId, 1+PTR_SIZE);

				if( !var->dataType.GetSubType().IsPrimitive() )
					lctx.bc.Instr(BC_PshRPtr);

				lctx.type.Set(var->dataType.GetSubType());

				if( !lctx.type.dataType.IsObject() || lctx.type.dataType.IsObjectHandle() )
					lctx.type.dataType.MakeReference(true);
				
				// If the element type is handles, then we're expected to do handle assignments
				if( lctx.type.dataType.IsObjectHandle() )
					lctx.type.isExplicitHandle = true;

				asSExprContext ctx(engine);
				DoAssignment(&ctx, &lctx, &rctx, el, el, ttAssignment, el);

				if( !lctx.type.dataType.IsPrimitive() )
					ctx.bc.Pop(PTR_SIZE);

				// Release temporary variables used by expression
				ReleaseTemporaryVariable(ctx.type, &ctx.bc);

				ProcessDeferredParams(&ctx);

				bc->AddCode(&ctx.bc);
			}

			el = el->next;
			index++;
		}
	}
	else
	{
		asCString str;
		str.Format(TXT_INIT_LIST_CANNOT_BE_USED_WITH_s, var->dataType.Format().AddressOf());
		Error(str.AddressOf(), node);
	}
}

void asCCompiler::CompileStatement(asCScriptNode *statement, bool *hasReturn, asCByteCode *bc)
{
	*hasReturn = false;

	if( statement->nodeType == snStatementBlock )
		CompileStatementBlock(statement, true, hasReturn, bc);
	else if( statement->nodeType == snIf )
		CompileIfStatement(statement, hasReturn, bc);
	else if( statement->nodeType == snFor )
		CompileForStatement(statement, bc);
	else if( statement->nodeType == snWhile )
		CompileWhileStatement(statement, bc);
	else if( statement->nodeType == snDoWhile )
		CompileDoWhileStatement(statement, bc);
	else if( statement->nodeType == snExpressionStatement )
		CompileExpressionStatement(statement, bc);
	else if( statement->nodeType == snBreak )
		CompileBreakStatement(statement, bc);
	else if( statement->nodeType == snContinue )
		CompileContinueStatement(statement, bc);
	else if( statement->nodeType == snSwitch )
		CompileSwitchStatement(statement, hasReturn, bc);
	else if( statement->nodeType == snReturn )
	{
		CompileReturnStatement(statement, bc);
		*hasReturn = true;
	}
}

void asCCompiler::CompileSwitchStatement(asCScriptNode *snode, bool *, asCByteCode *bc)
{
	// Reserve label for break statements
	int breakLabel = nextLabel++;
	breakLabels.PushLast(breakLabel);

	// Add a variable scope that will be used by CompileBreak
	// to know where to stop deallocating variables
	AddVariableScope(true, false);

	//---------------------------
	// Compile the switch expression
	//-------------------------------

	// Compile the switch expression
	asSExprContext expr(engine);
	CompileAssignment(snode->firstChild, &expr);

	// Verify that the expression is a primitive type
	if( !expr.type.dataType.IsIntegerType() && !expr.type.dataType.IsUnsignedType() )
	{
		Error(TXT_SWITCH_MUST_BE_INTEGRAL, snode->firstChild);
		return;
	}

	// TODO: Need to support 64bit 
	// Convert the expression to a 32bit variable
	asCDataType to;
	if( expr.type.dataType.IsIntegerType() )
		to.SetTokenType(ttInt);
	else if( expr.type.dataType.IsUnsignedType() )
		to.SetTokenType(ttUInt);
	ImplicitConversion(&expr, to, snode->firstChild, false, true);

	ConvertToVariable(&expr);
	int offset = expr.type.stackOffset;

	//-------------------------------
	// Determine case values and labels
	//--------------------------------

	// Remember the first label so that we can later pass the
	// correct label to each CompileCase()
	int firstCaseLabel = nextLabel;
	int defaultLabel = 0;

	asCArray<int> caseValues;
	asCArray<int> caseLabels;

	// Compile all case comparisons and make them jump to the right label
	asCScriptNode *cnode = snode->firstChild->next;
	while( cnode )
	{
		// Each case should have a constant expression
		if( cnode->firstChild && cnode->firstChild->nodeType == snExpression )
		{
			// Compile expression
			asSExprContext c(engine);
			CompileExpression(cnode->firstChild, &c);

			// Verify that the result is a constant
			if( !c.type.isConstant )
				Error(TXT_SWITCH_CASE_MUST_BE_CONSTANT, cnode->firstChild);

			// Verify that the result is an integral number
			if( !c.type.dataType.IsIntegerType() && !c.type.dataType.IsUnsignedType() )
				Error(TXT_SWITCH_MUST_BE_INTEGRAL, cnode->firstChild);

			// Store constant for later use
			int val;
			if( c.type.dataType.GetSizeInMemoryBytes() == 1 )
				val = (signed char)c.type.byteValue;
			else if( c.type.dataType.GetSizeInMemoryBytes() == 2 )
				val = (short)c.type.wordValue;
			else
				val = c.type.intValue;

			caseValues.PushLast(val);

			// Reserve label for this case
			caseLabels.PushLast(nextLabel++);
		}
		else
		{
			// Is default the last case?
			if( cnode->next )
			{
				Error(TXT_DEFAULT_MUST_BE_LAST, cnode);
				break;
			}

			// Reserve label for this case
			defaultLabel = nextLabel++;
		}

		cnode = cnode->next;
	}

    // check for empty switch
	if (caseValues.GetLength() == 0)
	{
		Error(TXT_EMPTY_SWITCH, snode);
		return;
	}

	if( defaultLabel == 0 )
		defaultLabel = breakLabel;

	//---------------------------------
    // Output the optimized case comparisons
	// with jumps to the case code
	//------------------------------------

	// Sort the case values by increasing value. Do the sort together with the labels
	// A simple bubble sort is sufficient since we don't expect a huge number of values
	for( asUINT fwd = 1; fwd < caseValues.GetLength(); fwd++ )
	{
		for( int bck = fwd - 1; bck >= 0; bck-- )
		{
			int bckp = bck + 1;
			if( caseValues[bck] > caseValues[bckp] )
			{
				// Swap the values in both arrays
				int swap = caseValues[bckp];
				caseValues[bckp] = caseValues[bck];
				caseValues[bck] = swap;

				swap = caseLabels[bckp];
				caseLabels[bckp] = caseLabels[bck];
				caseLabels[bck] = swap;
			}
			else
				break;
		}
	}

	// Find ranges of consecutive numbers
	asCArray<int> ranges;
	ranges.PushLast(0);
	asUINT n;
	for( n = 1; n < caseValues.GetLength(); ++n )
	{
		// We can join numbers that are less than 5 numbers
		// apart since the output code will still be smaller
		if( caseValues[n] > caseValues[n-1] + 5 )
			ranges.PushLast(n);
	}

	// If the value is larger than the largest case value, jump to default
	int tmpOffset = AllocateVariable(asCDataType::CreatePrimitive(ttInt, false), true);
	expr.bc.InstrSHORT_DW(BC_SetV4, (short)tmpOffset, caseValues[caseValues.GetLength()-1]);
	expr.bc.InstrW_W(BC_CMPi, offset, tmpOffset);
	expr.bc.InstrDWORD(BC_JP, defaultLabel);
	ReleaseTemporaryVariable(tmpOffset, &expr.bc);

	// TODO: We could possible optimize this even more by doing a
	// binary search instead of a linear search through the ranges

	// For each range
	int range;
	for( range = 0; range < (int)ranges.GetLength(); range++ )
	{
		// Find the largest value in this range
		int maxRange = caseValues[ranges[range]];
		int index = ranges[range];
		for( ; (index < (int)caseValues.GetLength()) && (caseValues[index] <= maxRange + 5); index++ )
			maxRange = caseValues[index];

		// If there are only 2 numbers then it is better to compare them directly
		if( index - ranges[range] > 2 )
		{
			// If the value is smaller than the smallest case value in the range, jump to default
			tmpOffset = AllocateVariable(asCDataType::CreatePrimitive(ttInt, false), true);
			expr.bc.InstrSHORT_DW(BC_SetV4, (short)tmpOffset, caseValues[ranges[range]]);
			expr.bc.InstrW_W(BC_CMPi, offset, tmpOffset);
			expr.bc.InstrDWORD(BC_JS, defaultLabel);
			ReleaseTemporaryVariable(tmpOffset, &expr.bc);

			int nextRangeLabel = nextLabel++;
			// If this is the last range we don't have to make this test
			if( range < (int)ranges.GetLength() - 1 )
			{
				// If the value is larger than the largest case value in the range, jump to the next range
				tmpOffset = AllocateVariable(asCDataType::CreatePrimitive(ttInt, false), true);
				expr.bc.InstrSHORT_DW(BC_SetV4, (short)tmpOffset, maxRange);
				expr.bc.InstrW_W(BC_CMPi, offset, tmpOffset);
				expr.bc.InstrDWORD(BC_JP, nextRangeLabel);
				ReleaseTemporaryVariable(tmpOffset, &expr.bc);
			}

			// Jump forward according to the value
			tmpOffset = AllocateVariable(asCDataType::CreatePrimitive(ttInt, false), true);
			expr.bc.InstrSHORT_DW(BC_SetV4, (short)tmpOffset, caseValues[ranges[range]]);
			expr.bc.InstrW_W_W(BC_SUBi, tmpOffset, offset, tmpOffset);
			ReleaseTemporaryVariable(tmpOffset, &expr.bc);
			expr.bc.JmpP(tmpOffset, maxRange - caseValues[ranges[range]]);

			// Add the list of jumps to the correct labels (any holes, jump to default)
			index = ranges[range];
			for( int n = caseValues[index]; n <= maxRange; n++ )
			{
				if( caseValues[index] == n )
					expr.bc.InstrINT(BC_JMP, caseLabels[index++]);
				else
					expr.bc.InstrINT(BC_JMP, defaultLabel);
			}

			expr.bc.Label((short)nextRangeLabel);
		}
		else
		{
			// Simply make a comparison with each value
			int n;
			for( n = ranges[range]; n < index; ++n )
			{
				tmpOffset = AllocateVariable(asCDataType::CreatePrimitive(ttInt, false), true);
				expr.bc.InstrSHORT_DW(BC_SetV4, (short)tmpOffset, caseValues[n]);
				expr.bc.InstrW_W(BC_CMPi, offset, tmpOffset);
				expr.bc.InstrDWORD(BC_JZ, caseLabels[n]);
				ReleaseTemporaryVariable(tmpOffset, &expr.bc);
			}
		}
	}

	// Catch any value that falls trough
	expr.bc.InstrINT(BC_JMP, defaultLabel);

	// Release the temporary variable previously stored
	ReleaseTemporaryVariable(expr.type, &expr.bc);

	//----------------------------------
    // Output case implementations
	//----------------------------------

	// Compile case implementations, each one with the label before it
	cnode = snode->firstChild->next;
	while( cnode )
	{
		// Each case should have a constant expression
		if( cnode->firstChild && cnode->firstChild->nodeType == snExpression )
		{
			expr.bc.Label((short)firstCaseLabel++);

			CompileCase(cnode->firstChild->next, &expr.bc);
		}
		else
		{
			expr.bc.Label((short)defaultLabel);

			// Is default the last case?
			if( cnode->next )
			{
				// We've already reported this error
				break;
			}

			CompileCase(cnode->firstChild, &expr.bc);
		}

		cnode = cnode->next;
	}

	//--------------------------------

	bc->AddCode(&expr.bc);

	// Add break label
	bc->Label((short)breakLabel);

	breakLabels.PopLast();
	RemoveVariableScope();
}

void asCCompiler::CompileCase(asCScriptNode *node, asCByteCode *bc)
{
	bool isFinished = false;
	bool hasReturn = false;
	while( node )
	{
		if( hasReturn || isFinished )
		{
			Warning(TXT_UNREACHABLE_CODE, node);
			break;
		}

		if( node->nodeType == snBreak || node->nodeType == snContinue )
			isFinished = true;

		asCByteCode statement(engine);
		CompileStatement(node, &hasReturn, &statement);

		LineInstr(bc, node->tokenPos);
		bc->AddCode(&statement);

		if( !hasCompileErrors )
			assert( tempVariables.GetLength() == 0 );

		node = node->next;
	}
}

void asCCompiler::CompileIfStatement(asCScriptNode *inode, bool *hasReturn, asCByteCode *bc)
{
	// We will use one label for the if statement
	// and possibly another for the else statement
	int afterLabel = nextLabel++;

	// Compile the expression
	asSExprContext expr(engine);
	CompileAssignment(inode->firstChild, &expr);
	if( !expr.type.dataType.IsEqualExceptRefAndConst(asCDataType::CreatePrimitive(ttBool, true)) )
		Error(TXT_EXPR_MUST_BE_BOOL, inode->firstChild);

	if( expr.type.dataType.IsReference() ) ConvertToVariable(&expr);
	ProcessDeferredParams(&expr);

	if( !expr.type.isConstant )
	{
		ConvertToVariable(&expr);

		// Add byte code from the expression
		bc->AddCode(&expr.bc);

		// Add a test
		bc->InstrSHORT(BC_CpyVtoR4, expr.type.stackOffset);
		bc->Instr(BC_ClrHi);
		bc->InstrDWORD(BC_JZ, afterLabel);
		ReleaseTemporaryVariable(expr.type, bc);
	}
	else if( expr.type.dwordValue == 0 )
	{
		// Jump to the else case
		bc->InstrINT(BC_JMP, afterLabel);

		// TODO: Should we warn?
	}

	// Compile the if statement
	bool hasReturn1;
	asCByteCode ifBC(engine);
	CompileStatement(inode->firstChild->next, &hasReturn1, &ifBC);

	// Add the byte code
	LineInstr(bc, inode->firstChild->next->tokenPos);
	bc->AddCode(&ifBC);

	// Do we have an else statement?
	if( inode->firstChild->next != inode->lastChild )
	{
		int afterElse = 0;
		if( !hasReturn1 )
		{
			afterElse = nextLabel++;

			// Add jump to after the else statement
			bc->InstrINT(BC_JMP, afterElse);
		}

		// Add label for the else statement
		bc->Label((short)afterLabel);

		bool hasReturn2;
		asCByteCode elseBC(engine);
		CompileStatement(inode->lastChild, &hasReturn2, &elseBC);

		// Add byte code for the else statement
		LineInstr(bc, inode->lastChild->tokenPos);
		bc->AddCode(&elseBC);

		if( !hasReturn1 )
		{
			// Add label for the end of else statement
			bc->Label((short)afterElse);
		}

		// The if statement only has return if both alternatives have
		*hasReturn = hasReturn1 && hasReturn2;
	}
	else
	{
		// Add label for the end of if statement
		bc->Label((short)afterLabel);
		*hasReturn = false;
	}
}

void asCCompiler::CompileForStatement(asCScriptNode *fnode, asCByteCode *bc)
{
	// Add a variable scope that will be used by CompileBreak/Continue to know where to stop deallocating variables
	AddVariableScope(true, true);

	// We will use three labels for the for loop
	int beforeLabel = nextLabel++;
	int afterLabel = nextLabel++;
	int continueLabel = nextLabel++;

	continueLabels.PushLast(continueLabel);
	breakLabels.PushLast(afterLabel);

	//---------------------------------------
	// Compile the initialization statement
	asCByteCode initBC(engine);
	asSExprContext expr(engine);

	if( fnode->firstChild->nodeType == snDeclaration )
		CompileDeclaration(fnode->firstChild, &initBC);
	else
		CompileExpressionStatement(fnode->firstChild, &initBC);

	//-----------------------------------
	// Compile the condition statement
	asCScriptNode *second = fnode->firstChild->next;
	if( second->firstChild )
	{
		CompileAssignment(second->firstChild, &expr);
		if( !expr.type.dataType.IsEqualExceptRefAndConst(asCDataType::CreatePrimitive(ttBool, true)) )
			Error(TXT_EXPR_MUST_BE_BOOL, second);

		if( expr.type.dataType.IsReference() ) ConvertToVariable(&expr);
		ProcessDeferredParams(&expr);

		// If expression is false exit the loop
		ConvertToVariable(&expr);
		expr.bc.InstrSHORT(BC_CpyVtoR4, expr.type.stackOffset);
		expr.bc.Instr(BC_ClrHi);
		expr.bc.InstrDWORD(BC_JZ, afterLabel);
		ReleaseTemporaryVariable(expr.type, &expr.bc);
	}

	//---------------------------
	// Compile the increment statement
	asSExprContext nextBC(engine);
	asCScriptNode *third = second->next;
	if( third->nodeType == snAssignment )
	{
		CompileAssignment(third, &nextBC);
		expr.type = nextBC.type;

		// Release temporary variables used by expression
		ReleaseTemporaryVariable(expr.type, &nextBC.bc);
	}

	//------------------------------
	// Compile loop statement
	bool hasReturn;
	asCByteCode forBC(engine);
	CompileStatement(fnode->lastChild, &hasReturn, &forBC);

	//-------------------------------
	// Join the code pieces
	bc->AddCode(&initBC);
	bc->Label((short)beforeLabel);

	// Add a suspend bytecode inside the loop to guarantee
	// that the application can suspend the execution
	bc->Instr(BC_SUSPEND);

	bc->AddCode(&expr.bc);
	LineInstr(bc, fnode->lastChild->tokenPos);
	bc->AddCode(&forBC);
	bc->Label((short)continueLabel);
	bc->AddCode(&nextBC.bc);
	bc->InstrINT(BC_JMP, beforeLabel);
	bc->Label((short)afterLabel);

	continueLabels.PopLast();
	breakLabels.PopLast();

	// Deallocate variables in this block, in reverse order
	for( int n = (int)variables->variables.GetLength() - 1; n >= 0; n-- )
	{
		sVariable *v = variables->variables[n];

		// Call variable destructors here, for variables not yet destroyed
		CompileDestructor(v->type, v->stackOffset, bc);

		// Don't deallocate function parameters
		if( v->stackOffset > 0 )
			DeallocateVariable(v->stackOffset);
	}

	RemoveVariableScope();
}

void asCCompiler::CompileWhileStatement(asCScriptNode *wnode, asCByteCode *bc)
{
	// Add a variable scope that will be used by CompileBreak/Continue to know where to stop deallocating variables
	AddVariableScope(true, true);

	// We will use two labels for the while loop
	int beforeLabel = nextLabel++;
	int afterLabel = nextLabel++;

	continueLabels.PushLast(beforeLabel);
	breakLabels.PushLast(afterLabel);

	// Add label before the expression
	bc->Label((short)beforeLabel);

	// Compile expression
	asSExprContext expr(engine);
	CompileAssignment(wnode->firstChild, &expr);
	if( !expr.type.dataType.IsEqualExceptRefAndConst(asCDataType::CreatePrimitive(ttBool, true)) )
		Error(TXT_EXPR_MUST_BE_BOOL, wnode->firstChild);

	if( expr.type.dataType.IsReference() ) ConvertToVariable(&expr);
	ProcessDeferredParams(&expr);

	// Add byte code for the expression
	ConvertToVariable(&expr);
	bc->AddCode(&expr.bc);

	// Jump to end of statement if expression is false
	bc->InstrSHORT(BC_CpyVtoR4, expr.type.stackOffset);
	bc->Instr(BC_ClrHi);
	bc->InstrDWORD(BC_JZ, afterLabel);
	ReleaseTemporaryVariable(expr.type, bc);

	// Add a suspend bytecode inside the loop to guarantee
	// that the application can suspend the execution
	bc->Instr(BC_SUSPEND);

	// Compile statement
	bool hasReturn;
	asCByteCode whileBC(engine);
	CompileStatement(wnode->lastChild, &hasReturn, &whileBC);

	// Add byte code for the statement
	LineInstr(bc, wnode->lastChild->tokenPos);
	bc->AddCode(&whileBC);

	// Jump to the expression
	bc->InstrINT(BC_JMP, beforeLabel);

	// Add label after the statement
	bc->Label((short)afterLabel);

	continueLabels.PopLast();
	breakLabels.PopLast();

	RemoveVariableScope();
}

void asCCompiler::CompileDoWhileStatement(asCScriptNode *wnode, asCByteCode *bc)
{
	// Add a variable scope that will be used by CompileBreak/Continue to know where to stop deallocating variables
	AddVariableScope(true, true);

	// We will use two labels for the while loop
	int beforeLabel = nextLabel++;
	int beforeTest = nextLabel++;
	int afterLabel = nextLabel++;

	continueLabels.PushLast(beforeTest);
	breakLabels.PushLast(afterLabel);

	// Add label before the statement
	bc->Label((short)beforeLabel);

	// Compile statement
	bool hasReturn;
	asCByteCode whileBC(engine);
	CompileStatement(wnode->firstChild, &hasReturn, &whileBC);

	// Add byte code for the statement
	LineInstr(bc, wnode->firstChild->tokenPos);
	bc->AddCode(&whileBC);

	// Add label before the expression
	bc->Label((short)beforeTest);

	// Add a suspend bytecode inside the loop to guarantee
	// that the application can suspend the execution
	bc->Instr(BC_SUSPEND);

	// Add a line instruction
	LineInstr(bc, wnode->lastChild->tokenPos);

	// Compile expression
	asSExprContext expr(engine);
	CompileAssignment(wnode->lastChild, &expr);
	if( !expr.type.dataType.IsEqualExceptRefAndConst(asCDataType::CreatePrimitive(ttBool, true)) )
		Error(TXT_EXPR_MUST_BE_BOOL, wnode->firstChild);

	if( expr.type.dataType.IsReference() ) ConvertToVariable(&expr);
	ProcessDeferredParams(&expr);

	// Add byte code for the expression
	ConvertToVariable(&expr);
	bc->AddCode(&expr.bc);

	// Jump to next iteration if expression is true
	bc->InstrSHORT(BC_CpyVtoR4, expr.type.stackOffset);
	bc->Instr(BC_ClrHi);
	bc->InstrDWORD(BC_JNZ, beforeLabel);
	ReleaseTemporaryVariable(expr.type, bc);

	// Add label after the statement
	bc->Label((short)afterLabel);

	continueLabels.PopLast();
	breakLabels.PopLast();

	RemoveVariableScope();
}

void asCCompiler::CompileBreakStatement(asCScriptNode *node, asCByteCode *bc)
{
	if( breakLabels.GetLength() == 0 )
	{
		Error(TXT_INVALID_BREAK, node);
		return;
	}

	// Add destructor calls for all variables that will go out of scope
	asCVariableScope *vs = variables;
	while( !vs->isBreakScope )
	{
		for( int n = (int)vs->variables.GetLength() - 1; n >= 0; n-- )
			CompileDestructor(vs->variables[n]->type, vs->variables[n]->stackOffset, bc);

		vs = vs->parent;
	}

	bc->InstrINT(BC_JMP, breakLabels[breakLabels.GetLength()-1]);
}

void asCCompiler::CompileContinueStatement(asCScriptNode *node, asCByteCode *bc)
{
	if( continueLabels.GetLength() == 0 )
	{
		Error(TXT_INVALID_CONTINUE, node);
		return;
	}

	// Add destructor calls for all variables that will go out of scope
	asCVariableScope *vs = variables;
	while( !vs->isContinueScope )
	{
		for( int n = (int)vs->variables.GetLength() - 1; n >= 0; n-- )
			CompileDestructor(vs->variables[n]->type, vs->variables[n]->stackOffset, bc);

		vs = vs->parent;
	}

	bc->InstrINT(BC_JMP, continueLabels[continueLabels.GetLength()-1]);
}

void asCCompiler::CompileExpressionStatement(asCScriptNode *enode, asCByteCode *bc)
{
	if( enode->firstChild )
	{
		// Compile the expression
		asSExprContext expr(engine);
		CompileAssignment(enode->firstChild, &expr);

		// Pop the value from the stack
		if( !expr.type.dataType.IsPrimitive() )
			expr.bc.Pop(expr.type.dataType.GetSizeOnStackDWords());

		// Release temporary variables used by expression
		ReleaseTemporaryVariable(expr.type, &expr.bc);

		ProcessDeferredParams(&expr);

		bc->AddCode(&expr.bc);
	}
}

void asCCompiler::PrepareTemporaryObject(asCScriptNode *node, asSExprContext *ctx)
{
	// If the object already is stored in temporary variable then nothing needs to be done
	if( ctx->type.isTemporary ) return;

	// Allocate temporary variable
	asCDataType dt = ctx->type.dataType;
	dt.MakeReference(false);
	dt.MakeReadOnly(false);

	int offset = AllocateVariable(dt, true);

	// Allocate and construct the temporary object
	CompileConstructor(dt, offset, &ctx->bc);

	// Assign the object to the temporary variable
	asCTypeInfo lvalue;
	dt.MakeReference(true);
	lvalue.Set(dt);
	lvalue.isTemporary = true;
	lvalue.stackOffset = (short)offset;
	lvalue.isVariable = true;
	lvalue.isExplicitHandle = ctx->type.isExplicitHandle;

	PrepareForAssignment(&lvalue.dataType, ctx, node);

	ctx->bc.InstrSHORT(BC_PSF, (short)offset);
	PerformAssignment(&lvalue, &ctx->type, &ctx->bc, node);

	// Pop the original reference
	ctx->bc.Pop(PTR_SIZE);

	// Push the reference to the temporary variable on the stack
	ctx->bc.InstrSHORT(BC_PSF, (short)offset);
	lvalue.dataType.MakeReference(true);

	ctx->type = lvalue;
}

void asCCompiler::CompileReturnStatement(asCScriptNode *rnode, asCByteCode *bc)
{
	// Get return type and location
	sVariable *v = variables->GetVariable("return");
	if( v->type.GetSizeOnStackDWords() > 0 )
	{
		// Is there an expression?
		if( rnode->firstChild )
		{
			// Compile the expression
			asSExprContext expr(engine);
			CompileAssignment(rnode->firstChild, &expr);

			// Prepare the value for assignment
			IsVariableInitialized(&expr.type, rnode->firstChild);

			if( v->type.IsPrimitive() )
			{
				if( expr.type.dataType.IsReference() ) ConvertToVariable(&expr);

				// Implicitly convert the value to the return type
				ImplicitConversion(&expr, v->type, rnode->firstChild, false);

				// Verify that the conversion was successful
				if( expr.type.dataType != v->type )
				{
					asCString str;
					str.Format(TXT_NO_CONVERSION_s_TO_s, expr.type.dataType.Format().AddressOf(), v->type.Format().AddressOf());
					Error(str.AddressOf(), rnode);
				}
				else
				{
					ConvertToVariable(&expr);
					ReleaseTemporaryVariable(expr.type, &expr.bc);

					// Load the variable in the register
					if( v->type.GetSizeOnStackDWords() == 1 )
						expr.bc.InstrSHORT(BC_CpyVtoR4, expr.type.stackOffset);
					else
						expr.bc.InstrSHORT(BC_CpyVtoR8, expr.type.stackOffset);
				}
			}
			else if( v->type.IsObject() )
			{
				PrepareArgument(&v->type, &expr, rnode->firstChild);

				// Pop the reference to the temporary variable again
				expr.bc.Pop(PTR_SIZE);

				// Load the object pointer into the object register
				expr.bc.InstrSHORT(BC_LOADOBJ, expr.type.stackOffset);

				// The temporary variable must not be freed as it will no longer hold an object
				DeallocateVariable(expr.type.stackOffset);
				expr.type.isTemporary = false;
			}

			// Release temporary variables used by expression
			ReleaseTemporaryVariable(expr.type, &expr.bc);

			bc->AddCode(&expr.bc);
		}
		else
			Error(TXT_MUST_RETURN_VALUE, rnode);
	}
	else
		if( rnode->firstChild )
			Error(TXT_CANT_RETURN_VALUE, rnode);

	// Call destructor on all variables except for the function parameters
	asCVariableScope *vs = variables;
	while( vs )
	{
		for( int n = (int)vs->variables.GetLength() - 1; n >= 0; n-- )
			if( vs->variables[n]->stackOffset > 0 )
				CompileDestructor(vs->variables[n]->type, vs->variables[n]->stackOffset, bc);

		vs = vs->parent;
	}

	// Jump to the end of the function
	bc->InstrINT(BC_JMP, 0);
}

void asCCompiler::AddVariableScope(bool isBreakScope, bool isContinueScope)
{
	variables = NEW(asCVariableScope)(variables);
	variables->isBreakScope    = isBreakScope;
	variables->isContinueScope = isContinueScope;
}

void asCCompiler::RemoveVariableScope()
{
	if( variables )
	{
		asCVariableScope *var = variables;
		variables = variables->parent;
		DELETE(var,asCVariableScope);
	}
}

void asCCompiler::Error(const char *msg, asCScriptNode *node)
{
	asCString str;

	int r, c;
	script->ConvertPosToRowCol(node->tokenPos, &r, &c);

	builder->WriteError(script->name.AddressOf(), msg, r, c);

	hasCompileErrors = true;
}

void asCCompiler::Warning(const char *msg, asCScriptNode *node)
{
	asCString str;

	int r, c;
	script->ConvertPosToRowCol(node->tokenPos, &r, &c);

	builder->WriteWarning(script->name.AddressOf(), msg, r, c);
}

int asCCompiler::AllocateVariable(const asCDataType &type, bool isTemporary)
{
	asCDataType t(type);

	if( t.IsPrimitive() && t.GetSizeOnStackDWords() == 1 )
		t.SetTokenType(ttInt);

	if( t.IsPrimitive() && t.GetSizeOnStackDWords() == 2 )
		t.SetTokenType(ttDouble);

	// Find a free location with the same type
	for( asUINT n = 0; n < freeVariables.GetLength(); n++ )
	{
		int slot = freeVariables[n];
		if( variableAllocations[slot].IsEqualExceptConst(t) && variableIsTemporary[slot] == isTemporary )
		{
			if( n != freeVariables.GetLength() - 1 )
				freeVariables[n] = freeVariables.PopLast();
			else
				freeVariables.PopLast();

			// We can't return by slot, must count variable sizes
			int offset = GetVariableOffset(slot);

			if( isTemporary )
				tempVariables.PushLast(offset);

			return offset;
		}
	}

	variableAllocations.PushLast(t);
	variableIsTemporary.PushLast(isTemporary);

	int offset = GetVariableOffset((int)variableAllocations.GetLength()-1);

	if( isTemporary )
		tempVariables.PushLast(offset);

	return offset;
}

// TODO: This function should share code with AllocateVariable
int asCCompiler::AllocateVariableNotIn(const asCDataType &type, bool isTemporary, asCArray<int> *vars)
{
	asCDataType t(type);

	if( t.IsPrimitive() && t.GetSizeOnStackDWords() == 1 )
		t.SetTokenType(ttInt);

	if( t.IsPrimitive() && t.GetSizeOnStackDWords() == 2 )
		t.SetTokenType(ttDouble);

	// Find a free location with the same type
	for( asUINT n = 0; n < freeVariables.GetLength(); n++ )
	{
		int slot = freeVariables[n];
		if( variableAllocations[slot].IsEqualExceptConst(t) && variableIsTemporary[slot] == isTemporary )
		{
			// We can't return by slot, must count variable sizes
			int offset = GetVariableOffset(slot);

			// Verify that it is not in the list of used variables
			bool isUsed = false;
			if( vars )
			{
				for( asUINT m = 0; m < vars->GetLength(); m++ )
				{
					if( offset == (*vars)[m] )
					{
						isUsed = true;
						break;
					}
				}
			}

			if( !isUsed )
			{
				if( n != freeVariables.GetLength() - 1 )
					freeVariables[n] = freeVariables.PopLast();
				else
					freeVariables.PopLast();

				if( isTemporary )
					tempVariables.PushLast(offset);

				return offset;
			}
		}
	}

	variableAllocations.PushLast(t);
	variableIsTemporary.PushLast(isTemporary);

	int offset = GetVariableOffset((int)variableAllocations.GetLength()-1);

	if( isTemporary )
		tempVariables.PushLast(offset);

	return offset;
}

int asCCompiler::GetVariableOffset(int varIndex)
{
	// Return offset to the last dword on the stack
	int varOffset = 1;
	for( int n = 0; n < varIndex; n++ )
		varOffset += variableAllocations[n].GetSizeOnStackDWords();

	if( varIndex < (int)variableAllocations.GetLength() )
	{
		int size = variableAllocations[varIndex].GetSizeOnStackDWords();
		if( size > 1 )
			varOffset += size-1;
	}

	return varOffset;
}

int asCCompiler::GetVariableSlot(int offset)
{
	int varOffset = 1;
	for( asUINT n = 0; n < variableAllocations.GetLength(); n++ )
	{
		varOffset += -1 + variableAllocations[n].GetSizeOnStackDWords();
		if( varOffset == offset )
		{
			return n;
		}
		varOffset++;
	}

	return -1;
}

void asCCompiler::DeallocateVariable(int offset)
{
	// Remove temporary variable
	int n;
	for( n = 0; n < (int)tempVariables.GetLength(); n++ )
	{
		if( offset == tempVariables[n] )
		{
			if( n == (int)tempVariables.GetLength()-1 )
				tempVariables.PopLast();
			else
				tempVariables[n] = tempVariables.PopLast();
			break;
		}
	}

	n = GetVariableSlot(offset);
	if( n != -1 )
	{
		freeVariables.PushLast(n);
		return;
	}

	// We might get here if the variable was implicitly declared
	// because it was use before a formal declaration, in this case
	// the offset is 0x7FFF

	assert(offset == 0x7FFF);
}

void asCCompiler::ReleaseTemporaryVariable(asCTypeInfo &t, asCByteCode *bc)
{
	if( t.isTemporary )
	{
		if( bc )
		{
			// We need to call the destructor on the true variable type
			int n = GetVariableSlot(t.stackOffset);
			asCDataType dt = variableAllocations[n];

			// Call destructor
			CompileDestructor(dt, t.stackOffset, bc);
		}

		DeallocateVariable(t.stackOffset);
		t.isTemporary = false;
	}
}

void asCCompiler::ReleaseTemporaryVariable(int offset, asCByteCode *bc)
{
	if( bc )
	{
		// We need to call the destructor on the true variable type
		int n = GetVariableSlot(offset);
		asCDataType dt = variableAllocations[n];

		// Call destructor
		CompileDestructor(dt, offset, bc);
	}

	DeallocateVariable(offset);
}

void asCCompiler::Dereference(asSExprContext *ctx, bool generateCode)
{
	if( ctx->type.dataType.IsReference() )
	{
		if( ctx->type.dataType.IsObject() )
		{
			ctx->type.dataType.MakeReference(false);
			if( generateCode ) 
			{
				ctx->bc.Instr(BC_CHKREF);
				ctx->bc.Instr(BC_RDSPTR);
			}
		}
		else
		{
			// This should never happen as primitives are treated differently
			assert(false);
		}
	}
}


bool asCCompiler::IsVariableInitialized(asCTypeInfo *type, asCScriptNode *node)
{
	// Temporary variables are assumed to be initialized
	if( type->isTemporary ) return true;

	// Verify that it is a variable
	if( !type->isVariable ) return true;

	// Find the variable
	sVariable *v = variables->GetVariableByOffset(type->stackOffset);

	// The variable isn't found if it is a constant, in which case it is guaranteed to be initialized
	if( v == 0 ) return true;

	if( v->isInitialized ) return true;

	// Complex types don't need this test
	if( v->type.IsObject() ) return true;

	// Mark as initialized so that the user will not be bothered again
	v->isInitialized = true;

	// Write warning
	asCString str;
	str.Format(TXT_s_NOT_INITIALIZED, (const char *)v->name.AddressOf());
	Warning(str.AddressOf(), node);

	return false;
}

void asCCompiler::PrepareOperand(asSExprContext *ctx, asCScriptNode *node)
{
	// Check if the variable is initialized (if it indeed is a variable)
	IsVariableInitialized(&ctx->type, node);

	asCDataType to = ctx->type.dataType;
	to.MakeReference(false);

	ImplicitConversion(ctx, to, node, false);

	ProcessDeferredParams(ctx);
}

void asCCompiler::PrepareForAssignment(asCDataType *lvalue, asSExprContext *rctx, asCScriptNode *node, asSExprContext *lvalueExpr)
{
	// Make sure the rvalue is initialized if it is a variable
	IsVariableInitialized(&rctx->type, node);

	if( lvalue->IsPrimitive() )
	{
		if( rctx->type.dataType.IsPrimitive() )
		{
			if( rctx->type.dataType.IsReference() )
			{
				// Cannot do implicit conversion of references so we first convert the reference to a variable
				ConvertToVariableNotIn(rctx, lvalueExpr);
			}

			// Implicitly convert the value to the right type
			asCArray<int> usedVars;
			if( lvalueExpr ) lvalueExpr->bc.GetVarsUsed(usedVars);
			ImplicitConversion(rctx, *lvalue, node, false, true, &usedVars);
		}

		// Check data type
		if( !lvalue->IsEqualExceptRefAndConst(rctx->type.dataType) )
		{
			asCString str;
			str.Format(TXT_CANT_IMPLICITLY_CONVERT_s_TO_s, rctx->type.dataType.Format().AddressOf(), lvalue->Format().AddressOf());
			Error(str.AddressOf(), node);
			
			rctx->type.SetDummy();
		}

		// Make sure the rvalue is a variable
		if( !rctx->type.isVariable )
			ConvertToVariableNotIn(rctx, lvalueExpr);
	}
	else
	{
		asCDataType to = *lvalue;
		to.MakeReference(false);

		// Don't allow the implicit conversion to create an object
		ImplicitConversion(rctx, to, node, false, true, 0, false);

		// Check data type
		if( lvalue->IsObjectHandle() &&
			!rctx->type.dataType.IsObjectHandle() )
		{
			Error(TXT_NEED_TO_BE_A_HANDLE, node);
		}
		else if( !lvalue->IsEqualExceptRefAndConst(rctx->type.dataType) )
		{
			asCString str;
			str.Format(TXT_CANT_IMPLICITLY_CONVERT_s_TO_s, rctx->type.dataType.Format().AddressOf(), lvalue->Format().AddressOf());
			Error(str.AddressOf(), node);
		}
		else
		{
			// If the assignment will be made with the copy behaviour then the rvalue must not be a reference
			if( lvalue->IsObject() )
				assert(!rctx->type.dataType.IsReference());
		}
	}
}

bool asCCompiler::IsLValue(asCTypeInfo &type)
{
	if( type.dataType.IsReadOnly() ) return false;
	if( !type.dataType.IsObject() && !type.isVariable && !type.dataType.IsReference() ) return false;
	if( type.isTemporary ) return false;
	return true;
}

void asCCompiler::PerformAssignment(asCTypeInfo *lvalue, asCTypeInfo *rvalue, asCByteCode *bc, asCScriptNode *node)
{
	if( lvalue->dataType.IsReadOnly() )
		Error(TXT_REF_IS_READ_ONLY, node);

	if( lvalue->dataType.IsPrimitive() )
	{
		if( lvalue->isVariable )
		{
			// Copy the value between the variables directly
			if( lvalue->dataType.GetSizeInMemoryDWords() == 1 )
				bc->InstrW_W(BC_CpyVtoV4, lvalue->stackOffset, rvalue->stackOffset);
			else
				bc->InstrW_W(BC_CpyVtoV8, lvalue->stackOffset, rvalue->stackOffset);

			// Mark variable as initialized
			sVariable *v = variables->GetVariableByOffset(lvalue->stackOffset);
			if( v ) v->isInitialized = true;
		}
		else if( lvalue->dataType.IsReference() )
		{
			// Copy the value of the variable to the reference in the register
			int s = lvalue->dataType.GetSizeInMemoryBytes();
			if( s == 1 )
				bc->InstrSHORT(BC_WRTV1, rvalue->stackOffset);
			else if( s == 2 )
				bc->InstrSHORT(BC_WRTV2, rvalue->stackOffset);
			else if( s == 4 )
				bc->InstrSHORT(BC_WRTV4, rvalue->stackOffset);
			else if( s == 8 )
				bc->InstrSHORT(BC_WRTV8, rvalue->stackOffset);
		}
		else
		{
			Error(TXT_NOT_VALID_LVALUE, node);
			return;
		}
	}
	else if( !lvalue->isExplicitHandle )
	{
		// TODO: Call the assignment operator, or do a BC_COPY if none exist

		// TODO: Convert
		asSExprContext ctx(engine);
		ctx.type = *lvalue;
		Dereference(&ctx, true);
		*lvalue = ctx.type;
		bc->AddCode(&ctx.bc);

		// TODO: Can't this leave deferred output params unhandled?

		asSTypeBehaviour *beh = lvalue->dataType.GetBehaviour();
		if( beh->copy )
		{
			// Call the copy operator
			bc->Call(BC_CALLSYS, (asDWORD)beh->copy, 2*PTR_SIZE);
			bc->Instr(BC_PshRPtr);
		}
		else
		{
			// Default copy operator
			if( lvalue->dataType.GetSizeInMemoryDWords() == 0 )
			{
				Error(TXT_NO_DEFAULT_COPY_OP, node);
			}

			// Copy larger data types from a reference
			bc->InstrWORD(BC_COPY, (asWORD)lvalue->dataType.GetSizeInMemoryDWords());
		}
	}
	else
	{
		// TODO: Convert
		// TODO: The object handle can be stored in a variable as well
		if( !lvalue->dataType.IsReference() )
		{
			Error(TXT_NOT_VALID_REFERENCE, node);
			return;
		}

		// TODO: Convert to register based
		bc->InstrPTR(BC_REFCPY, builder->module->RefObjectType(lvalue->dataType.GetObjectType()));

		// Mark variable as initialized
		if( variables )
		{
			sVariable *v = variables->GetVariableByOffset(lvalue->stackOffset);
			if( v ) v->isInitialized = true;
		}
	}
}

void asCCompiler::ImplicitConversion(asSExprContext *ctx, const asCDataType &to, asCScriptNode *node, bool isExplicit, bool generateCode, asCArray<int> *reservedVars, bool allowObjectConstruct)
{
	// No conversion from void to any other type
	if( ctx->type.dataType.GetTokenType() == ttVoid ) return;
	
	// Do we want a var type?
	if( to.GetTokenType() == ttQuestion )
	{
		// Any type can be converted to a var type, but only when not generating code
		assert( !generateCode );

		ctx->type.dataType = to;

		return;
	}
	// Do we want a primitive?
	else if( to.IsPrimitive() )
	{
		// No conversion from objects to primitives yet
		if( !ctx->type.dataType.IsPrimitive() ) return;

		// Start by implicitly converting constant values
		if( ctx->type.isConstant )
			ImplicitConversionConstant(ctx, to, node, isExplicit);

		if( to == ctx->type.dataType )
			return;

		// After the constant value has been converted we have the following possibilities

		// Allow implicit conversion between numbers
		if( generateCode )
		{
			// Convert smaller types to 32bit first
			int s = ctx->type.dataType.GetSizeInMemoryBytes();
			if( s < 4 )
			{
				ConvertToTempVariableNotIn(ctx, reservedVars);
				if( ctx->type.dataType.IsIntegerType() )
				{
					if( s == 1 )
						ctx->bc.InstrSHORT(BC_sbTOi, ctx->type.stackOffset);
					else if( s == 2 )
						ctx->bc.InstrSHORT(BC_swTOi, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(ttInt);
				}
				else if( ctx->type.dataType.IsUnsignedType() )
				{
					if( s == 1 )
						ctx->bc.InstrSHORT(BC_ubTOi, ctx->type.stackOffset);
					else if( s == 2 )
						ctx->bc.InstrSHORT(BC_uwTOi, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(ttUInt);
				}
			}

			if( to.IsIntegerType() && to.GetSizeInMemoryDWords() == 1 )
			{
				if( ctx->type.dataType.IsIntegerType() || 
					ctx->type.dataType.IsUnsignedType() )
				{
					if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
						ctx->type.dataType.SetTokenType(to.GetTokenType());
					else
					{
						ConvertToTempVariableNotIn(ctx, reservedVars);
						ReleaseTemporaryVariable(ctx->type, &ctx->bc);
						int offset = AllocateVariableNotIn(to, true, reservedVars);
						ctx->bc.InstrW_W(BC_i64TOi, offset, ctx->type.stackOffset);
						ctx->type.SetVariable(to, offset, true);
					}
				}
				else if( ctx->type.dataType.IsFloatType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_fTOi, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
				else if( ctx->type.dataType.IsDoubleType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_dTOi, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}

				// Convert to smaller integer if necessary
				int s = to.GetSizeInMemoryBytes();
				if( s < 4 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					if( s == 1 )
						ctx->bc.InstrSHORT(BC_iTOb, ctx->type.stackOffset);
					else if( s == 2 )
						ctx->bc.InstrSHORT(BC_iTOw, ctx->type.stackOffset);
				}
			}
			if( to.IsIntegerType() && to.GetSizeInMemoryDWords() == 2 )
			{
				if( ctx->type.dataType.IsIntegerType() || 
					ctx->type.dataType.IsUnsignedType() )
				{
					if( ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
						ctx->type.dataType.SetTokenType(to.GetTokenType());
					else
					{
						ConvertToTempVariableNotIn(ctx, reservedVars);
						ReleaseTemporaryVariable(ctx->type, &ctx->bc);
						int offset = AllocateVariableNotIn(to, true, reservedVars);
						if( ctx->type.dataType.IsUnsignedType() )
							ctx->bc.InstrW_W(BC_uTOi64, offset, ctx->type.stackOffset);
						else
							ctx->bc.InstrW_W(BC_iTOi64, offset, ctx->type.stackOffset);
						ctx->type.SetVariable(to, offset, true);
					}
				}
				else if( ctx->type.dataType.IsFloatType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_fTOi64, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
				else if( ctx->type.dataType.IsDoubleType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_dTOi64, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
			}
			else if( to.IsUnsignedType() && to.GetSizeInMemoryDWords() == 1  )
			{
				if( ctx->type.dataType.IsIntegerType() || 
					ctx->type.dataType.IsUnsignedType() )
				{
					if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
						ctx->type.dataType.SetTokenType(to.GetTokenType());
					else
					{
						ConvertToTempVariableNotIn(ctx, reservedVars);
						ReleaseTemporaryVariable(ctx->type, &ctx->bc);
						int offset = AllocateVariableNotIn(to, true, reservedVars);
						ctx->bc.InstrW_W(BC_i64TOi, offset, ctx->type.stackOffset);
						ctx->type.SetVariable(to, offset, true);
					}
				}
				else if( ctx->type.dataType.IsFloatType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_fTOu, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
				else if( ctx->type.dataType.IsDoubleType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_dTOu, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}

				// Convert to smaller integer if necessary
				int s = to.GetSizeInMemoryBytes();
				if( s < 4 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					if( s == 1 )
						ctx->bc.InstrSHORT(BC_iTOb, ctx->type.stackOffset);
					else if( s == 2 )
						ctx->bc.InstrSHORT(BC_iTOw, ctx->type.stackOffset);
				}
			}
			if( to.IsUnsignedType() && to.GetSizeInMemoryDWords() == 2 )
			{
				if( ctx->type.dataType.IsIntegerType() || 
					ctx->type.dataType.IsUnsignedType() )
				{
					if( ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
						ctx->type.dataType.SetTokenType(to.GetTokenType());
					else
					{
						ConvertToTempVariableNotIn(ctx, reservedVars);
						ReleaseTemporaryVariable(ctx->type, &ctx->bc);
						int offset = AllocateVariableNotIn(to, true, reservedVars);
						if( ctx->type.dataType.IsUnsignedType() )
							ctx->bc.InstrW_W(BC_uTOi64, offset, ctx->type.stackOffset);
						else
							ctx->bc.InstrW_W(BC_iTOi64, offset, ctx->type.stackOffset);
						ctx->type.SetVariable(to, offset, true);
					}
				}
				else if( ctx->type.dataType.IsFloatType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_fTOu64, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
				else if( ctx->type.dataType.IsDoubleType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_dTOu64, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
			}
			else if( to.IsFloatType() )
			{
				if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_iTOf, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
				else if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_i64TOf, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
				else if( ctx->type.dataType.IsUnsignedType() && ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_uTOf, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
				else if( ctx->type.dataType.IsUnsignedType() && ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_u64TOf, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
				else if( ctx->type.dataType.IsDoubleType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_dTOf, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
			}
			else if( to.IsDoubleType() )
			{
				if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_iTOd, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
				else if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_i64TOd, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
				else if( ctx->type.dataType.IsUnsignedType() && ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_uTOd, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
				else if( ctx->type.dataType.IsUnsignedType() && ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ctx->bc.InstrSHORT(BC_u64TOd, ctx->type.stackOffset);
					ctx->type.dataType.SetTokenType(to.GetTokenType());
				}
				else if( ctx->type.dataType.IsFloatType() )
				{
					ConvertToTempVariableNotIn(ctx, reservedVars);
					ReleaseTemporaryVariable(ctx->type, &ctx->bc);
					int offset = AllocateVariableNotIn(to, true, reservedVars);
					ctx->bc.InstrW_W(BC_fTOd, offset, ctx->type.stackOffset);
					ctx->type.SetVariable(to, offset, true);
				}
			}
		}
		else
		{
			if( (to.IsIntegerType() || to.IsUnsignedType() ||
				 to.IsFloatType()   || to.IsDoubleType()) &&
				(ctx->type.dataType.IsIntegerType() || ctx->type.dataType.IsUnsignedType() ||
				 ctx->type.dataType.IsFloatType()   || ctx->type.dataType.IsDoubleType()) )
				ctx->type.dataType.SetTokenType(to.GetTokenType());
		}

		// Primitive types on the stack, can be const or non-const
		ctx->type.dataType.MakeReadOnly(to.IsReadOnly());
	}
	else // The target is a complex type
	{
		// Convert null to any object type handle
		if( ctx->type.IsNullConstant() && to.IsObjectHandle() )
			ctx->type.dataType = to;

		if( to.IsObjectHandle() )
		{
			// An object type can be directly converted to a handle of the same type
			if( ctx->type.dataType.SupportHandles() )
			{
				ctx->type.dataType.MakeHandle(true);
			}

			if( ctx->type.dataType.IsObjectHandle() )
				ctx->type.dataType.MakeReadOnly(to.IsReadOnly());

			if( to.IsHandleToConst() && ctx->type.dataType.IsObjectHandle() )
				ctx->type.dataType.MakeHandleToConst(true);
		}

		if( !to.IsReference() )
		{
			if( to.IsObjectHandle() )
			{
				if( ctx->type.dataType.IsReference() )
				{
					Dereference(ctx, generateCode);

					// TODO: Can't this leave unhandled deferred output params?
				}

				// TODO: If the type is handle, then we can't use IsReadOnly to determine the constness of the basetype

				// If the rvalue is a handle to a const object, then
				// the lvalue must also be a handle to a const object
				if( ctx->type.dataType.IsReadOnly() && !to.IsReadOnly() )
				{
					if( isExplicit )
					{
						assert(node);
						asCString str;
						str.Format(TXT_CANT_IMPLICITLY_CONVERT_s_TO_s, ctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
						Error(str.AddressOf(), node);
					}
				}

				// If the to type is an interface and the from type implements it, then we can convert it immediately
				if( ctx->type.dataType.GetObjectType() && 
					ctx->type.dataType.GetObjectType()->Implements(to.GetObjectType()) )
				{
					assert(ctx->type.dataType.IsObjectHandle());
					ctx->type.dataType.SetObjectType(to.GetObjectType());
				}
			}
			else
			{
				ImplicitConversionToObject(ctx, to, node, isExplicit, generateCode, reservedVars, allowObjectConstruct);
			}
		}
		else // to.IsReference()
		{
			if( ctx->type.dataType.IsReference() )
			{
				// A reference to a handle can be converted to a reference to an object
				// by first reading the address, then verifying that it is not null, then putting the address back on the stack
				if( !to.IsObjectHandle() && ctx->type.dataType.IsObjectHandle() && !ctx->type.isExplicitHandle )
				{
					ctx->type.dataType.MakeHandle(false);
					if( generateCode )
						ctx->bc.Instr(BC_ChkRefS);
				}

				// A reference to a non-const can be converted to a reference to a const
				if( to.IsReadOnly() )
					ctx->type.dataType.MakeReadOnly(true);
				else if( ctx->type.dataType.IsReadOnly() )
				{
					// A reference to a const can be converted to a reference to a non-const by putting the object in a temporary variable
					ctx->type.dataType.MakeReadOnly(false);

					if( generateCode )
					{
						assert(!ctx->type.isTemporary);

						// Allocate a temporary variable
						int offset = AllocateVariableNotIn(ctx->type.dataType, true, reservedVars);
						ctx->type.isTemporary = true;
						ctx->type.stackOffset = (short)offset;

						CompileConstructor(ctx->type.dataType, offset, &ctx->bc);

						asSExprContext rctx(engine);
						rctx.type = ctx->type;
						rctx.bc.AddCode(&ctx->bc);
						asSExprContext lctx(engine);
						lctx.type = ctx->type;
						lctx.bc.InstrSHORT(BC_PSF, (short)offset);
						DoAssignment(ctx, &lctx, &rctx, node, node, ttAssignment, node);
					}
				}
			}
			else
			{
				if( generateCode )
				{
					asCTypeInfo type;
					type.Set(ctx->type.dataType);

					// Allocate a temporary variable
					int offset = AllocateVariableNotIn(type.dataType, true, reservedVars);
					type.isTemporary = true;
					type.stackOffset = (short)offset;
					if( type.dataType.IsObjectHandle() )
						type.isExplicitHandle = true;

					CompileConstructor(type.dataType, offset, &ctx->bc);
					type.dataType.MakeReference(true);

					PrepareForAssignment(&type.dataType, ctx, node);

					ctx->bc.InstrSHORT(BC_PSF, type.stackOffset);

					// If the input type is read-only we'll need to temporarily remove this constness, otherwise the assignment will fail
					bool typeIsReadOnly = type.dataType.IsReadOnly();
					type.dataType.MakeReadOnly(false);
					PerformAssignment(&type, &ctx->type, &ctx->bc, node);
					type.dataType.MakeReadOnly(typeIsReadOnly);

					ctx->bc.Pop(ctx->type.dataType.GetSizeOnStackDWords());

					ReleaseTemporaryVariable(ctx->type, &ctx->bc);

					ctx->bc.InstrSHORT(BC_PSF, type.stackOffset);

					ctx->type = type;
				}

				// A non-reference can be converted to a reference,
				// by putting the value in a temporary variable
				ctx->type.dataType.MakeReference(true);

				// Since it is a new temporary variable it doesn't have to be const
				ctx->type.dataType.MakeReadOnly(to.IsReadOnly());
			}
		}
	}
}

void asCCompiler::ImplicitConversionToObject(asSExprContext *ctx, const asCDataType &to, asCScriptNode *node, bool /*isExplicit*/, bool generateCode, asCArray<int> * /*reservedVars*/, bool allowObjectConstruct)
{
	if( ctx->type.dataType.IsReference() && !ctx->type.dataType.IsPrimitive() )
	{
		Dereference(ctx, generateCode);

		// TODO: Can't this leave unhandled deferred output params?
	}

	if( ctx->type.dataType.IsObject() && ctx->type.dataType.GetObjectType() == to.GetObjectType() )
	{
		if( ctx->type.dataType.IsObjectHandle() && !ctx->type.isExplicitHandle )
		{
			if( ctx->type.dataType.IsReference() )
			{
				if( generateCode ) ctx->bc.Instr(BC_RDSPTR);
				ctx->type.dataType.MakeReference(false);
			}

			if( generateCode )
				ctx->bc.Instr(BC_CHKREF);

			ctx->type.dataType.MakeHandle(false);
		}
	}
	else if( allowObjectConstruct )
	{
		// Since the expression is not of the same object type we need to check if there
		// is any constructor that can be used to create an object of the correct type.

		asCArray<int> funcs;
		asSTypeBehaviour *beh = to.GetBehaviour();
		if( beh )
			funcs = beh->constructors;

		// Compile the arguments
		asCArray<asSExprContext *> args;
		asCArray<asCTypeInfo> temporaryVariables;

		args.PushLast(ctx);

		MatchFunctions(funcs, args, node, to.GetObjectType()->name.AddressOf(), false, true, false);

		// Verify that we found 1 matching function
		if( funcs.GetLength() == 1 )
		{
			asCTypeInfo tempObj;
			tempObj.dataType = to;
			tempObj.dataType.MakeReference(true);
			tempObj.isTemporary = true;
			tempObj.isVariable = true;

			if( generateCode )
			{
				tempObj.stackOffset = (short)AllocateVariable(to, true);

				asSExprContext tmp(engine);

				// Push the address of the object on the stack
				tmp.bc.InstrSHORT(BC_VAR, tempObj.stackOffset);

				PrepareFunctionCall(funcs[0], &tmp.bc, args);

				MoveArgsToStack(funcs[0], &tmp.bc, args, false);

				int offset = 0;
				for( asUINT n = 0; n < args.GetLength(); n++ )
					offset += args[n]->type.dataType.GetSizeOnStackDWords();

				tmp.bc.InstrWORD(BC_GETREF, (asWORD)offset);

				PerformFunctionCall(funcs[0], &tmp, true, &args, tempObj.dataType.GetObjectType());

				// The constructor doesn't return anything,
				// so we have to manually inform the type of
				// the return value
				tmp.type = tempObj;

				// Push the address of the object on the stack again
				tmp.bc.InstrSHORT(BC_PSF, tempObj.stackOffset);

				// Copy the newly generated code to the input context
				// ctx is already empty, since it was merged as part of argument expression
				assert(ctx->bc.GetLastInstr() == -1);
				MergeExprContexts(ctx, &tmp);
			}

			ctx->type = tempObj;
		}
	}
}

void asCCompiler::ImplicitConversionConstant(asSExprContext *from, const asCDataType &to, asCScriptNode *node, bool isExplicit)
{
	assert(from->type.isConstant);

	// TODO: node should be the node of the value that is
	// converted (not the operator that provokes the implicit
	// conversion)

	// If the base type is correct there is no more to do
	if( to.IsEqualExceptRefAndConst(from->type.dataType) ) return;

	// References cannot be constants
	if( from->type.dataType.IsReference() ) return;

	// Arrays can't be constants
	if( to.IsArrayType() ) return;

	if( to.IsIntegerType() && to.GetSizeInMemoryDWords() == 1 )
	{
		// Float constants can be implicitly converted to int
		if( from->type.dataType.IsFloatType() )
		{
			float fc = from->type.floatValue;
			int ic = int(fc);

			if( float(ic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);
			from->type.intValue = ic;

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		// Double constants can be implicitly converted to int
		else if( from->type.dataType.IsDoubleType() )
		{
			double fc = from->type.doubleValue;
			int ic = int(fc);

			if( double(ic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);
			from->type.intValue = ic;

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsUnsignedType() && from->type.dataType.GetSizeInMemoryDWords() == 1 )
		{
			// Verify that it is possible to convert to signed without getting negative
			if( from->type.intValue < 0 )
			{
				if( !isExplicit && node ) Warning(TXT_CHANGE_SIGN, node);
			}

			// Convert to 32bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.intValue = from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.intValue = from->type.wordValue;

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsUnsignedType() && from->type.dataType.GetSizeInMemoryDWords() == 2 )
		{
			// Convert to 32bit
			from->type.intValue = int(from->type.qwordValue);

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsIntegerType() && 
		         from->type.dataType.GetSizeInMemoryBytes() < 4 )
		{
			// Convert to 32bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.intValue = (signed char)from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.intValue = (short)from->type.wordValue;

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsIntegerType() &&
		         from->type.dataType.GetSizeInMemoryBytes() > to.GetSizeInMemoryBytes() )
		{
			// Verify if it is possible
			if( to.GetSizeInMemoryBytes() == 1 )
			{
				if( char(from->type.intValue) != from->type.intValue )
					if( !isExplicit && node ) Warning(TXT_VALUE_TOO_LARGE_FOR_TYPE, node);

				from->type.byteValue = char(from->type.intValue);
			}
			else if( to.GetSizeInMemoryBytes() == 2 )
			{
				if( short(from->type.intValue) != from->type.intValue )
					if( !isExplicit && node ) Warning(TXT_VALUE_TOO_LARGE_FOR_TYPE, node);

				from->type.wordValue = short(from->type.intValue);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
		}
	}
	else if( to.IsIntegerType() && to.GetSizeInMemoryDWords() == 2 )
	{
		// Float constants can be implicitly converted to int
		if( from->type.dataType.IsFloatType() )
		{
			float fc = from->type.floatValue;
			asINT64 ic = asINT64(fc);

			if( float(ic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt64, true);
			from->type.qwordValue = ic;
		}
		// Double constants can be implicitly converted to int
		else if( from->type.dataType.IsDoubleType() )
		{
			double fc = from->type.doubleValue;
			asINT64 ic = asINT64(fc);

			if( double(ic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt64, true);
			from->type.qwordValue = ic;
		}
		else if( from->type.dataType.IsUnsignedType() )
		{
			// Convert to 64bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.qwordValue = from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.qwordValue = from->type.wordValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 4 )
				from->type.qwordValue = from->type.dwordValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 8 )
			{
				if( asINT64(from->type.qwordValue) < 0 )
				{
					if( !isExplicit && node ) Warning(TXT_CHANGE_SIGN, node);
				}
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt64, true);
		}
		else if( from->type.dataType.IsIntegerType() )
		{
			// Convert to 64bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.qwordValue = (signed char)from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.qwordValue = (short)from->type.wordValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 4 )
				from->type.qwordValue = from->type.intValue;

			from->type.dataType = asCDataType::CreatePrimitive(ttInt64, true);
		}
	}
	else if( to.IsUnsignedType() && to.GetSizeInMemoryDWords() == 1 )
	{
		if( from->type.dataType.IsFloatType() )
		{
			float fc = from->type.floatValue;
			int uic = int(fc);

			if( float(uic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}
			else if( uic < 0 )
			{
				if( !isExplicit && node ) Warning(TXT_CHANGE_SIGN, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);
			from->type.intValue = uic;

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsDoubleType() )
		{
			double fc = from->type.doubleValue;
			int uic = int(fc);

			if( double(uic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttInt, true);
			from->type.intValue = uic;

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsIntegerType() )
		{
			// Verify that it is possible to convert to unsigned without loosing negative
			if( from->type.intValue < 0 )
			{
				if( !isExplicit && node ) Warning(TXT_CHANGE_SIGN, node);
			}

			// Convert to 32bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.intValue = (signed char)from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.intValue = (short)from->type.wordValue;

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt, true);

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsUnsignedType() && 
		         from->type.dataType.GetSizeInMemoryBytes() < 4 )
		{
			// Convert to 32bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.dwordValue = from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.dwordValue = from->type.wordValue;

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt, true);

			// Try once more, in case of a smaller type
			ImplicitConversionConstant(from, to, node, isExplicit);
		}
		else if( from->type.dataType.IsUnsignedType() &&
		         from->type.dataType.GetSizeInMemoryBytes() > to.GetSizeInMemoryBytes() )
		{
			// Verify if it is possible
			if( to.GetSizeInMemoryBytes() == 1 )
			{
				if( asBYTE(from->type.dwordValue) != from->type.dwordValue )
					if( !isExplicit && node ) Warning(TXT_VALUE_TOO_LARGE_FOR_TYPE, node);

				from->type.byteValue = asBYTE(from->type.dwordValue);
			}
			else if( to.GetSizeInMemoryBytes() == 2 )
			{
				if( asWORD(from->type.dwordValue) != from->type.dwordValue )
					if( !isExplicit && node ) Warning(TXT_VALUE_TOO_LARGE_FOR_TYPE, node);

				from->type.wordValue = asWORD(from->type.dwordValue);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
		}
	}
	else if( to.IsUnsignedType() && to.GetSizeInMemoryDWords() == 2 )
	{
		if( from->type.dataType.IsFloatType() )
		{
			float fc = from->type.floatValue;
			// Convert first to int64 then to uint64 to avoid negative float becoming 0 on gnuc base compilers
			asQWORD uic = asQWORD(asINT64(fc));

			// TODO: MSVC6 doesn't permit UINT64 to double
			if( float((signed)uic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt64, true);
			from->type.qwordValue = uic;
		}
		else if( from->type.dataType.IsDoubleType() )
		{
			double fc = from->type.doubleValue;
			// Convert first to int64 then to uint64 to avoid negative float becoming 0 on gnuc base compilers
			asQWORD uic = asQWORD(asINT64(fc));

			// TODO: MSVC6 doesn't permit UINT64 to double
			if( double((signed)uic) != fc )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt64, true);
			from->type.qwordValue = uic;
		}
		else if( from->type.dataType.IsIntegerType() && from->type.dataType.GetSizeInMemoryDWords() == 1 )
		{
			// Convert to 64bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.qwordValue = (asINT64)(signed char)from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.qwordValue = (asINT64)(short)from->type.wordValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 4 )
				from->type.qwordValue = (asINT64)from->type.intValue;

			// Verify that it is possible to convert to unsigned without loosing negative
			if( asINT64(from->type.qwordValue) < 0 )
			{
				if( !isExplicit && node ) Warning(TXT_CHANGE_SIGN, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt64, true);
		}
		else if( from->type.dataType.IsIntegerType() && from->type.dataType.GetSizeInMemoryDWords() == 2 )
		{
			// Verify that it is possible to convert to unsigned without loosing negative
			if( asINT64(from->type.qwordValue) < 0 )
			{
				if( !isExplicit && node ) Warning(TXT_CHANGE_SIGN, node);
			}

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt64, true);
		}
		else if( from->type.dataType.IsUnsignedType() )
		{
			// Convert to 64bit
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				from->type.qwordValue = from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				from->type.qwordValue = from->type.wordValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 4 )
				from->type.qwordValue = from->type.dwordValue;

			from->type.dataType = asCDataType::CreatePrimitive(ttUInt64, true);
		}
	}
	else if( to.IsFloatType() )
	{
		if( from->type.dataType.IsDoubleType() )
		{
			double ic = from->type.doubleValue;
			float fc = float(ic);

			if( double(fc) != ic )
			{
				asCString str;
				str.Format(TXT_POSSIBLE_LOSS_OF_PRECISION);
				if( !isExplicit && node ) Warning(str.AddressOf(), node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.floatValue = fc;
		}
		else if( from->type.dataType.IsIntegerType() && from->type.dataType.GetSizeInMemoryDWords() == 1 )
		{
			// Must properly convert value in case the from value is smaller
			int ic;
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				ic = (signed char)from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				ic = (short)from->type.wordValue;
			else
				ic = from->type.intValue;
			float fc = float(ic);

			if( int(fc) != ic )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.floatValue = fc;
		}
		else if( from->type.dataType.IsIntegerType() && from->type.dataType.GetSizeInMemoryDWords() == 2 )
		{
			float fc = float(asINT64(from->type.qwordValue));
			if( asINT64(fc) != asINT64(from->type.qwordValue) )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.floatValue = fc;
		}
		else if( from->type.dataType.IsUnsignedType() && from->type.dataType.GetSizeInMemoryDWords() == 1 )
		{
			// Must properly convert value in case the from value is smaller
			unsigned int uic;
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				uic = from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				uic = from->type.wordValue;
			else
				uic = from->type.dwordValue;
			float fc = float(uic);

			if( (unsigned int)(fc) != uic )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.floatValue = fc;
		}
		else if( from->type.dataType.IsUnsignedType() && from->type.dataType.GetSizeInMemoryDWords() == 2 )
		{
			// TODO: MSVC6 doesn't permit UINT64 to double
			float fc = float((signed)from->type.qwordValue);

			if( asQWORD(fc) != from->type.qwordValue )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.floatValue = fc;
		}
	}
	else if( to.IsDoubleType() )
	{
		if( from->type.dataType.IsFloatType() )
		{
			float ic = from->type.floatValue;
			double fc = double(ic);

			// Don't check for float->double
		//	if( float(fc) != ic )
		//	{
		//		acCString str;
		//		str.Format(TXT_NOT_EXACT_g_g_g, ic, fc, float(fc));
		//		if( !isExplicit ) Warning(str, node);
		//	}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.doubleValue = fc;
		}
		else if( from->type.dataType.IsIntegerType() && from->type.dataType.GetSizeInMemoryDWords() == 1 )
		{
			// Must properly convert value in case the from value is smaller
			int ic;
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				ic = (signed char)from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				ic = (short)from->type.wordValue;
			else
				ic = from->type.intValue;
			double fc = double(ic);

			if( int(fc) != ic )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.doubleValue = fc;
		}
		else if( from->type.dataType.IsIntegerType() && from->type.dataType.GetSizeInMemoryDWords() == 2 )
		{
			double fc = double(asINT64(from->type.qwordValue));

			if( asINT64(fc) != asINT64(from->type.qwordValue) )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.doubleValue = fc;
		}
		else if( from->type.dataType.IsUnsignedType() && from->type.dataType.GetSizeInMemoryDWords() == 1 )
		{
			// Must properly convert value in case the from value is smaller
			unsigned int uic;
			if( from->type.dataType.GetSizeInMemoryBytes() == 1 )
				uic = from->type.byteValue;
			else if( from->type.dataType.GetSizeInMemoryBytes() == 2 )
				uic = from->type.wordValue;
			else
				uic = from->type.dwordValue;
			double fc = double(uic);

			if( (unsigned int)(fc) != uic )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.doubleValue = fc;
		}
		else if( from->type.dataType.IsUnsignedType() && from->type.dataType.GetSizeInMemoryDWords() == 2 )
		{
			// TODO: MSVC6 doesn't permit UINT64 to double
			double fc = double((signed)from->type.qwordValue);

			if( asQWORD(fc) != from->type.qwordValue )
			{
				if( !isExplicit && node ) Warning(TXT_NOT_EXACT, node);
			}

			from->type.dataType.SetTokenType(to.GetTokenType());
			from->type.doubleValue = fc;
		}
	}
}

void asCCompiler::DoAssignment(asSExprContext *ctx, asSExprContext *lctx, asSExprContext *rctx, asCScriptNode *lexpr, asCScriptNode *rexpr, int op, asCScriptNode *opNode)
{
	if( lctx->type.dataType.IsPrimitive() )
	{
		if( op != ttAssignment )
		{
			// Compute the operator before the assignment
			asCTypeInfo lvalue = lctx->type;
			asSExprContext o(engine);
			CompileOperator(opNode, lctx, rctx, &o);
			MergeExprContexts(rctx, &o);
			rctx->type = o.type;

			// Convert the rvalue to the right type and validate it
			PrepareForAssignment(&lvalue.dataType, rctx, rexpr);

			MergeExprContexts(ctx, rctx);
			lctx->type = lvalue;

			// The lvalue continues the same, either it was a variable, or a reference in the register
		}
		else
		{
			// Convert the rvalue to the right type and validate it
			PrepareForAssignment(&lctx->type.dataType, rctx, rexpr, lctx);

			MergeExprContexts(ctx, rctx);
			MergeExprContexts(ctx, lctx);
		}

		ReleaseTemporaryVariable(rctx->type, &ctx->bc);

		PerformAssignment(&lctx->type, &rctx->type, &ctx->bc, opNode);

		ctx->type = lctx->type;
	}
	else if( lctx->type.isExplicitHandle )
	{
		// Verify that the left hand value isn't a temporary variable
		if( lctx->type.isTemporary )
		{
			Error(TXT_REF_IS_TEMP, lexpr);
			return;
		}

		// Object handles don't have any compound assignment operators
		if( op != ttAssignment )
		{
			asCString str;
			str.Format(TXT_ILLEGAL_OPERATION_ON_s, lctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), lexpr);
			return;
		}

		asCDataType dt = lctx->type.dataType;
		dt.MakeReference(false);

		PrepareArgument(&dt, rctx, rexpr, true, 1);
		if( !dt.IsEqualExceptRefAndConst(rctx->type.dataType) )
		{
			asCString str;
			str.Format(TXT_CANT_IMPLICITLY_CONVERT_s_TO_s, rctx->type.dataType.Format().AddressOf(), lctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), rexpr);
		}

		MergeExprContexts(ctx, rctx);
		MergeExprContexts(ctx, lctx);

		ctx->bc.InstrWORD(BC_GETOBJREF, PTR_SIZE);

		PerformAssignment(&lctx->type, &rctx->type, &ctx->bc, opNode);

		ReleaseTemporaryVariable(rctx->type, &ctx->bc);

		ctx->type = rctx->type;
	}
	else // if( lctx->type.dataType.IsObject() )
	{
		// Verify that the left hand value isn't a temporary variable
		if( lctx->type.isTemporary )
		{
			Error(TXT_REF_IS_TEMP, lexpr);
			return;
		}

		if( lctx->type.dataType.IsObjectHandle() && !lctx->type.isExplicitHandle )
		{
			// Convert the handle to a object reference
			asCDataType to;
			to = lctx->type.dataType;
			to.MakeHandle(false);
			ImplicitConversion(lctx, to, lexpr, false);
		}

		// If left expression resolves into a registered type
		// check if the assignment operator is overloaded, and check
		// the type of the right hand expression. If none is found
		// the default action is a direct copy if it is the same type
		// and a simple assignment.
		asSTypeBehaviour *beh = lctx->type.dataType.GetBehaviour();
		assert(beh);

		// Find the matching overloaded operators
		asCArray<int> ops;
		asUINT n;
		for( n = 0; n < beh->operators.GetLength(); n += 2 )
		{
			if( op == beh->operators[n] )
				ops.PushLast(beh->operators[n+1]);
		}

		asCArray<int> match;
		MatchArgument(ops, match, &rctx->type, 0);

		if( match.GetLength() == 1 )
		{
			// If it is an array, both sides must have the same subtype
			if( lctx->type.dataType.IsArrayType() )
				if( !lctx->type.dataType.IsEqualExceptRefAndConst(rctx->type.dataType) )
					Error(TXT_BOTH_MUST_BE_SAME, opNode);

			// We must verify that the lvalue isn't const
			if( lctx->type.dataType.IsReadOnly() )
				Error(TXT_REF_IS_READ_ONLY, lexpr);

			// Prepare the rvalue
			asCScriptFunction *descr = engine->scriptFunctions[match[0]];
			PrepareArgument(&descr->parameterTypes[0], rctx, rexpr, true, descr->inOutFlags[0]);

			if( rctx->type.isTemporary && lctx->bc.IsVarUsed(rctx->type.stackOffset) )
			{
				// Release the current temporary variable
				ReleaseTemporaryVariable(rctx->type, 0);

				asCArray<int> usedVars;
				lctx->bc.GetVarsUsed(usedVars);
				rctx->bc.GetVarsUsed(usedVars);

				asCDataType dt = rctx->type.dataType;
				dt.MakeReference(false);
				int newOffset = AllocateVariableNotIn(dt, true, &usedVars);

				rctx->bc.ExchangeVar(rctx->type.stackOffset, newOffset);
				rctx->type.stackOffset = (short)newOffset;
				rctx->type.isTemporary = true;
				rctx->type.isVariable = true;
			}

			// Add code for arguments
			MergeExprContexts(ctx, rctx);

			// Add the code for the object
			Dereference(lctx, true);
			MergeExprContexts(ctx, lctx);

			asCArray<asSExprContext*> args;
			args.PushLast(rctx);
			MoveArgsToStack(match[0], &ctx->bc, args, true);

			PerformFunctionCall(match[0], ctx, false, &args);

			return;
		}
		else if( match.GetLength() > 1 )
		{
			Error(TXT_MORE_THAN_ONE_MATCHING_OP, opNode);

			ctx->type.Set(lctx->type.dataType);

			return;
		}

		// No registered operator was found. In case the operation is a direct
		// assignment and the rvalue is the same type as the lvalue, then we can
		// still use the byte-for-byte copy to do the assignment

		if( op != ttAssignment )
		{
			asCString str;
			str.Format(TXT_ILLEGAL_OPERATION_ON_s, lctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), lexpr);
			return;
		}

		// Implicitly convert the rvalue to the type of the lvalue
		asCDataType dt = lctx->type.dataType;
		PrepareArgument(&dt, rctx, rexpr, true, 1);
		if( !dt.IsEqualExceptRefAndConst(rctx->type.dataType) )
		{
			asCString str;
			str.Format(TXT_CANT_IMPLICITLY_CONVERT_s_TO_s, rctx->type.dataType.Format().AddressOf(), lctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), rexpr);
		}

		MergeExprContexts(ctx, rctx);
		MergeExprContexts(ctx, lctx);

		ctx->bc.InstrWORD(BC_GETOBJREF, PTR_SIZE);

		PerformAssignment(&lctx->type, &rctx->type, &ctx->bc, opNode);

		ReleaseTemporaryVariable(rctx->type, &ctx->bc);

		ctx->type = rctx->type;
	}
}

void asCCompiler::CompileAssignment(asCScriptNode *expr, asSExprContext *ctx)
{
	asCScriptNode *lexpr = expr->firstChild;
	if( lexpr->next )
	{
		if( globalExpression )
		{
			Error(TXT_ASSIGN_IN_GLOBAL_EXPR, expr);
			ctx->type.SetDummy();
			return;
		}

		int currNumErrors = builder->numErrors;

		// Compile the two expression terms
		asSExprContext lctx(engine), rctx(engine);
		CompileAssignment(lexpr->next->next, &rctx);
		CompileCondition(lexpr, &lctx);

		if( currNumErrors == builder->numErrors )
			DoAssignment(ctx, &lctx, &rctx, lexpr, lexpr->next->next, lexpr->next->tokenType, lexpr->next);
		else
		{
			// Since the operands failed, the assignment was not computed
			ctx->type.SetDummy();
			return;
		}
	}
	else
		CompileCondition(lexpr, ctx);
}

// TODO: Convert
void asCCompiler::CompileCondition(asCScriptNode *expr, asSExprContext *ctx)
{
	asCTypeInfo ctype;

	// Compile the conditional expression
	asCScriptNode *cexpr = expr->firstChild;
	if( cexpr->next )
	{
		//-------------------------------
		// Compile the expressions
		asSExprContext e(engine);
		CompileExpression(cexpr, &e);
		ctype = e.type;
		if( !ctype.dataType.IsEqualExceptRefAndConst(asCDataType::CreatePrimitive(ttBool, true)) )
			Error(TXT_EXPR_MUST_BE_BOOL, cexpr);

		if( e.type.dataType.IsReference() ) ConvertToVariable(&e);
		ProcessDeferredParams(&e);

		asSExprContext le(engine);
		CompileAssignment(cexpr->next, &le);

		asSExprContext re(engine);
		CompileAssignment(cexpr->next->next, &re);

		bool isExplicitHandle = le.type.isExplicitHandle || re.type.isExplicitHandle;

		// Allow a 0 in the first case to be implicitly converted to the second type
		if( le.type.isConstant && le.type.intValue == 0 && le.type.dataType.IsUnsignedType() )
		{
			asCDataType to = re.type.dataType;
			to.MakeReference(false);
			to.MakeReadOnly(true);
			ImplicitConversionConstant(&le, to, cexpr->next, false);
		}

		//---------------------------------
		// Output the byte code
		int afterLabel = nextLabel++;
		int elseLabel = nextLabel++;

		// If left expression is void, then we don't need to store the result
		if( le.type.dataType.IsEqualExceptConst(asCDataType::CreatePrimitive(ttVoid, false)) )
		{
			// Put the code for the condition expression on the output
			MergeExprContexts(ctx, &e);
			
			// Added the branch decision
			ctx->type = e.type;
			ConvertToVariable(ctx);
			ctx->bc.InstrSHORT(BC_CpyVtoR4, ctx->type.stackOffset);
			ctx->bc.Instr(BC_ClrHi);
			ctx->bc.InstrDWORD(BC_JZ, elseLabel);
			ReleaseTemporaryVariable(ctx->type, &ctx->bc);

			// Add the left expression
			MergeExprContexts(ctx, &le);
			ctx->bc.InstrINT(BC_JMP, afterLabel);

			// Add the right expression
			ctx->bc.Label((short)elseLabel);
			MergeExprContexts(ctx, &re);
			ctx->bc.Label((short)afterLabel);

			// Make sure both expressions have the same type
			if( le.type.dataType != re.type.dataType )
				Error(TXT_BOTH_MUST_BE_SAME, expr);

			// Set the type of the result
			ctx->type = le.type;
		}
		else
		{
			// Allocate temporary variable and copy the result to that one
			asCTypeInfo temp;
			temp = le.type;
			temp.dataType.MakeReference(false);
			temp.dataType.MakeReadOnly(false);
			// Make sure the variable isn't used in the initial expression
			asCArray<int> vars;
			e.bc.GetVarsUsed(vars);
			int offset = AllocateVariableNotIn(temp.dataType, true, &vars);
			temp.SetVariable(temp.dataType, offset, true);

			CompileConstructor(temp.dataType, offset, &ctx->bc);

			// Put the code for the condition expression on the output
			MergeExprContexts(ctx, &e);

			// Added the branch decision
			ctx->type = e.type;
			ConvertToVariable(ctx);
			ctx->bc.InstrSHORT(BC_CpyVtoR4, ctx->type.stackOffset);
			ctx->bc.Instr(BC_ClrHi);
			ctx->bc.InstrDWORD(BC_JZ, elseLabel);
			ReleaseTemporaryVariable(ctx->type, &ctx->bc);

			// Assign the result of the left expression to the temporary variable
			asCTypeInfo rtemp;
			rtemp = temp;
			if( rtemp.dataType.IsObjectHandle() )
				rtemp.isExplicitHandle = true;

			PrepareForAssignment(&rtemp.dataType, &le, cexpr->next);
			MergeExprContexts(ctx, &le);

			if( !rtemp.dataType.IsPrimitive() )
			{
				ctx->bc.InstrSHORT(BC_PSF, (short)offset);
				rtemp.dataType.MakeReference(true);
			}
			PerformAssignment(&rtemp, &le.type, &ctx->bc, cexpr->next);
			if( !rtemp.dataType.IsPrimitive() )
				ctx->bc.Pop(le.type.dataType.GetSizeOnStackDWords()); // Pop the original value

			// Release the old temporary variable
			ReleaseTemporaryVariable(le.type, &ctx->bc);

			ctx->bc.InstrINT(BC_JMP, afterLabel);

			// Start of the right expression
			ctx->bc.Label((short)elseLabel);

			// Copy the result to the same temporary variable
			PrepareForAssignment(&rtemp.dataType, &re, cexpr->next);
			MergeExprContexts(ctx, &re);

			if( !rtemp.dataType.IsPrimitive() )
			{
				ctx->bc.InstrSHORT(BC_PSF, (short)offset);
				rtemp.dataType.MakeReference(true);
			}
			PerformAssignment(&rtemp, &re.type, &ctx->bc, cexpr->next);
			if( !rtemp.dataType.IsPrimitive() )
				ctx->bc.Pop(le.type.dataType.GetSizeOnStackDWords()); // Pop the original value

			// Release the old temporary variable
			ReleaseTemporaryVariable(re.type, &ctx->bc);

			ctx->bc.Label((short)afterLabel);

			// Make sure both expressions have the same type
			if( le.type.dataType != re.type.dataType )
				Error(TXT_BOTH_MUST_BE_SAME, expr);

			// Set the temporary variable as output
			ctx->type = rtemp;
			ctx->type.isExplicitHandle = isExplicitHandle;

			if( !ctx->type.dataType.IsPrimitive() )
			{
				ctx->bc.InstrSHORT(BC_PSF, (short)offset);
				ctx->type.dataType.MakeReference(true);
			}

			// Make sure the output isn't marked as being a literal constant
			ctx->type.isConstant = false;
		}
	}
	else
		CompileExpression(cexpr, ctx);
}

void asCCompiler::CompileExpression(asCScriptNode *expr, asSExprContext *ctx)
{
	assert(expr->nodeType == snExpression);

	// Count the nodes
	int count = 0;
	asCScriptNode *node = expr->firstChild;
	while( node )
	{
		count++;
		node = node->next;
	}

	// Convert to polish post fix, i.e: a+b => ab+
	asCArray<asCScriptNode *> stack(count);
	asCArray<asCScriptNode *> stack2(count);
	asCArray<asCScriptNode *> postfix(count);

	node = expr->firstChild;
	while( node )
	{
		int precedence = GetPrecedence(node);

		while( stack.GetLength() > 0 &&
			   precedence <= GetPrecedence(stack[stack.GetLength()-1]) )
			stack2.PushLast(stack.PopLast());

		stack.PushLast(node);

		node = node->next;
	}

	while( stack.GetLength() > 0 )
		stack2.PushLast(stack.PopLast());

	// We need to swap operands so that the left
	// operand is always computed before the right
	SwapPostFixOperands(stack2, postfix);

	// Compile the postfix formatted expression
	CompilePostFixExpression(&postfix, ctx);
}

void asCCompiler::SwapPostFixOperands(asCArray<asCScriptNode *> &postfix, asCArray<asCScriptNode *> &target)
{
	if( postfix.GetLength() == 0 ) return;

	asCScriptNode *node = postfix.PopLast();
	if( node->nodeType == snExprTerm )
	{
		target.PushLast(node);
		return;
	}

	SwapPostFixOperands(postfix, target);
	SwapPostFixOperands(postfix, target);

	target.PushLast(node);
}

void asCCompiler::CompilePostFixExpression(asCArray<asCScriptNode *> *postfix, asSExprContext *ctx)
{
	// Shouldn't send any byte code
	assert(ctx->bc.GetLastInstr() == -1);

	// Pop the last node
	asCScriptNode *node = postfix->PopLast();
	ctx->exprNode = node;

	// If term, compile the term
	if( node->nodeType == snExprTerm )
	{
		CompileExpressionTerm(node, ctx);
		return;
	}

	// Compile the two expression terms
	asSExprContext r(engine), l(engine);

	CompilePostFixExpression(postfix, &l);
	CompilePostFixExpression(postfix, &r);

	// Compile the operation
	CompileOperator(node, &l, &r, ctx);
}

void asCCompiler::CompileExpressionTerm(asCScriptNode *node, asSExprContext *ctx)
{
	// Shouldn't send any byte code
	assert(ctx->bc.GetLastInstr() == -1);

	// Compile the value node
	asCScriptNode *vnode = node->firstChild;
	while( vnode->nodeType != snExprValue )
		vnode = vnode->next;

	asSExprContext v(engine);
	CompileExpressionValue(vnode, &v);

	// Compile post fix operators
	asCScriptNode *pnode = vnode->next;
	while( pnode )
	{
		CompileExpressionPostOp(pnode, &v);
		pnode = pnode->next;
	}

	// Compile pre fix operators
	pnode = vnode->prev;
	while( pnode )
	{
		CompileExpressionPreOp(pnode, &v);
		pnode = pnode->prev;
	}

	// Return the byte code and final type description
	MergeExprContexts(ctx, &v);

    ctx->type = v.type;
}

void asCCompiler::CompileExpressionValue(asCScriptNode *node, asSExprContext *ctx)
{
	// Shouldn't receive any byte code
	assert(ctx->bc.GetLastInstr() == -1);

	asCScriptNode *vnode = node->firstChild;
	if( vnode->nodeType == snIdentifier )
	{
		GETSTRING(name, &script->code[vnode->tokenPos], vnode->tokenLength);

		sVariable *v = variables->GetVariable(name.AddressOf());
		if( v == 0 )
		{
			bool found = false;
			if( outFunc && outFunc->objectType )
			{
				if( name == THIS_TOKEN )
				{
					asCDataType dt = asCDataType::CreateObject(outFunc->objectType, false);

					// The object pointer is located at stack position 0
					ctx->bc.InstrSHORT(BC_PSF, 0);
					ctx->type.SetVariable(dt, 0, false);
					ctx->type.dataType.MakeReference(true);

					found = true;
				}
	
				if( !found )
				{
					asCDataType dt = asCDataType::CreateObject(outFunc->objectType, false);
					asCProperty *prop = builder->GetObjectProperty(dt, name.AddressOf());
					if( prop )
					{
						// The object pointer is located at stack position 0
						ctx->bc.InstrSHORT(BC_PSF, 0);
						ctx->type.SetVariable(dt, 0, false);
						ctx->type.dataType.MakeReference(true);

						Dereference(ctx, true);

						// TODO: This is the same as what is in CompileExpressionPostOp
						// Put the offset on the stack
						ctx->bc.InstrINT(BC_ADDSi, prop->byteOffset);

						if( prop->type.IsReference() )
							ctx->bc.Instr(BC_RDSPTR);

						// Reference to primitive must be stored in the temp register
						if( prop->type.IsPrimitive() )
						{
							// The ADD offset command should store the reference in the register directly
							ctx->bc.Instr(BC_PopRPtr);
						}

						// Set the new type (keeping info about temp variable)
						ctx->type.dataType = prop->type;
						ctx->type.dataType.MakeReference(true);
						ctx->type.isVariable = false;

						if( ctx->type.dataType.IsObject() && !ctx->type.dataType.IsObjectHandle() )
						{
							// Objects that are members are not references
							ctx->type.dataType.MakeReference(false);
						}

						// TODO: Check constness
						//ctx->type.dataType.MakeReadOnly(isConst ? true : prop->type.IsReadOnly());

						found = true;
					}
				}
			}

			if( !found )
			{
				// Is it a global property?
				bool isCompiled = true;
				bool isPureConstant = false;
				asQWORD constantValue;
				asCProperty *prop = builder->GetGlobalProperty(name.AddressOf(), &isCompiled, &isPureConstant, &constantValue);
				if( prop )
				{
					// Verify that the global property has been compiled already
					if( isCompiled )
					{
						// If the global property is a pure constant
						// we can allow the compiler to optimize it. Pure
						// constants are global constant variables that were
						// initialized by literal constants.
						if( isPureConstant )
							ctx->type.SetConstantQW(prop->type, constantValue);
						else
						{
							ctx->type.Set(prop->type);
							ctx->type.dataType.MakeReference(true);

							// Push the address of the variable on the stack
							if( ctx->type.dataType.IsPrimitive() )
								ctx->bc.InstrWORD(BC_LDG, (asWORD)builder->module->GetGlobalVarIndex(prop->index));
							else
								ctx->bc.InstrWORD(BC_PGA, (asWORD)builder->module->GetGlobalVarIndex(prop->index));

							// If the global property is registered by the application
							// then module must keep a reference to the config group
							builder->module->RefConfigGroupForGlobalVar(prop->index);
						}
					}
					else
					{
						asCString str;
						str.Format(TXT_UNINITIALIZED_GLOBAL_VAR_s, prop->name.AddressOf());
						Error(str.AddressOf(), vnode);
					}
				}
				else
				{
					asCString str;
					str.Format(TXT_s_NOT_DECLARED, name.AddressOf());
					Error(str.AddressOf(), vnode);

					// Give dummy value
					ctx->type.SetDummy();

					// Declare the variable now so that it will not be reported again
					variables->DeclareVariable(name.AddressOf(), asCDataType::CreatePrimitive(ttInt, false), 0x7FFF);

					// Mark the variable as initialized so that the user will not be bother by it again
					sVariable *v = variables->GetVariable(name.AddressOf());
					assert(v);
					if( v ) v->isInitialized = true;
				}
			}
		}
		else
		{
			if( v->isPureConstant )
				ctx->type.SetConstantQW(v->type, v->constantValue);
			else
			{
				if( v->type.IsPrimitive() )
				{
					if( v->type.IsReference() )
					{
						ctx->bc.InstrSHORT(BC_CpyVtoR4, (short)v->stackOffset);
						ctx->type.Set(v->type);
					}
					else
						ctx->type.SetVariable(v->type, v->stackOffset, false);
				}
				else
				{
					ctx->bc.InstrSHORT(BC_PSF, (short)v->stackOffset);
					ctx->type.SetVariable(v->type, v->stackOffset, false);
					ctx->type.dataType.MakeReference(true);

					// Implicitly dereference primitive parameters sent by reference
					if( v->type.IsReference() && (!v->type.IsObject() || v->type.IsObjectHandle()) )
						ctx->bc.Instr(BC_RDSPTR);
				}
			}
		}
	}
	else if( vnode->nodeType == snConstant )
	{
		if( vnode->tokenType == ttIntConstant )
		{
			GETSTRING(value, &script->code[vnode->tokenPos], vnode->tokenLength);

			asQWORD val = 0;

			// Process either character literal or numeric constant
			if( value[0] == '\'' )
			{
				asCString cat;
				cat = value.SubString(1, value.GetLength()-2);
				ProcessStringConstant(cat);
				if( cat.GetLength() == 0 )
					Error(TXT_EMPTY_CHAR_LITERAL, vnode);
				else
					val = (unsigned)cat[0];
			}
			else
				val = asStringScanUInt64(value.AddressOf(), 10, 0);

			// Do we need 64 bits?
			if( val>>32 )
				ctx->type.SetConstantQW(asCDataType::CreatePrimitive(ttUInt64, true), val);
			else
				ctx->type.SetConstantDW(asCDataType::CreatePrimitive(ttUInt, true), asDWORD(val));
		}
		else if( vnode->tokenType == ttBitsConstant )
		{
			GETSTRING(value, &script->code[vnode->tokenPos+2], vnode->tokenLength-2);

			// TODO: Check for overflow
			asQWORD val = asStringScanUInt64(value.AddressOf(), 16, 0);

			// Do we need 64 bits?
			if( val>>32 )
				ctx->type.SetConstantQW(asCDataType::CreatePrimitive(ttUInt64, true), val);
			else
				ctx->type.SetConstantDW(asCDataType::CreatePrimitive(ttUInt, true), asDWORD(val));
		}
		else if( vnode->tokenType == ttFloatConstant )
		{
			GETSTRING(value, &script->code[vnode->tokenPos], vnode->tokenLength);

			// TODO: Check for overflow
			float v = float(asStringScanDouble(value.AddressOf(), 0));
			ctx->type.SetConstantF(asCDataType::CreatePrimitive(ttFloat, true), v);
		}
		else if( vnode->tokenType == ttDoubleConstant )
		{
			GETSTRING(value, &script->code[vnode->tokenPos], vnode->tokenLength);

			// TODO: Check for overflow
			double v = asStringScanDouble(value.AddressOf(), 0);
			ctx->type.SetConstantD(asCDataType::CreatePrimitive(ttDouble, true), v);
		}
		else if( vnode->tokenType == ttTrue ||
			     vnode->tokenType == ttFalse )
		{
#if AS_SIZEOF_BOOL == 1 
			ctx->type.SetConstantB(asCDataType::CreatePrimitive(ttBool, true), vnode->tokenType == ttTrue ? VALUE_OF_BOOLEAN_TRUE : 0);
#else
			ctx->type.SetConstantDW(asCDataType::CreatePrimitive(ttBool, true), vnode->tokenType == ttTrue ? VALUE_OF_BOOLEAN_TRUE : 0);
#endif
		}
		else if( vnode->tokenType == ttStringConstant || vnode->tokenType == ttHeredocStringConstant )
		{
			asCString str;

			asCScriptNode *snode = vnode->firstChild;
			while( snode )
			{
				asCString cat;
				if( snode->tokenType == ttStringConstant )
				{
					cat.Assign(&script->code[snode->tokenPos+1], snode->tokenLength-2);
					ProcessStringConstant(cat);
				}
				else if( snode->tokenType == ttHeredocStringConstant )
				{
					cat.Assign(&script->code[snode->tokenPos+3], snode->tokenLength-6);
					ProcessHeredocStringConstant(cat);
				}

				str += cat;

				snode = snode->next;
			}

			// Call the string factory function to create a string object
			asCScriptFunction *descr = engine->stringFactory;
			if( descr == 0 )
			{
				// Error
				Error(TXT_STRINGS_NOT_RECOGNIZED, vnode);

				// Give dummy value
				ctx->type.SetDummy();
			}
			else
			{
				// Register the constant string with the engine
				int id = builder->module->AddConstantString(str.AddressOf(), str.GetLength());
				ctx->bc.InstrWORD(BC_STR, (asWORD)id);
				PerformFunctionCall(descr->id, ctx);
			}
		}
		else if( vnode->tokenType == ttNull )
		{
#ifndef AS_64BIT_PTR
			ctx->bc.InstrDWORD(BC_PshC4, 0);
#else
			ctx->bc.InstrQWORD(BC_SET8, 0);
#endif
			ctx->type.SetNullConstant();
		}
		else
			assert(false);
	}
	else if( vnode->nodeType == snFunctionCall )
	{
		bool found = false;
		if( outFunc && outFunc->objectType )
		{
			// Check if a class method is being called
			asCScriptNode *nm = vnode->firstChild;
			asCString name;
			name.Assign(&script->code[nm->tokenPos], nm->tokenLength);
			asCArray<int> funcs;
			// TODO: Verify the constness
			builder->GetObjectMethodDescriptions(name.AddressOf(), outFunc->objectType, funcs, false);
			if( funcs.GetLength() )
			{
				asCDataType dt = asCDataType::CreateObject(outFunc->objectType, false);

				// The object pointer is located at stack position 0
				ctx->bc.InstrSHORT(BC_PSF, 0);
				ctx->type.SetVariable(dt, 0, false);
				ctx->type.dataType.MakeReference(true);

				Dereference(ctx, true);

				CompileFunctionCall(vnode, ctx, outFunc->objectType, false);
				found = true;
			}
		}

		if( !found )
			CompileFunctionCall(vnode, ctx, 0, false);
	}
	else if( vnode->nodeType == snConstructCall )
	{
		CompileConstructCall(vnode, ctx);
	}
	else if( vnode->nodeType == snAssignment )
	{
		asSExprContext e(engine);
		CompileAssignment(vnode, &e);
		MergeExprContexts(ctx, &e);
		ctx->type = e.type;
	}
	else if( vnode->nodeType == snCast )
	{
		// Implement the cast operator
		CompileConversion(vnode, ctx);
	}
	else
		assert(false);
}

void asCCompiler::ProcessStringConstant(asCString &cstr)
{
	// Process escape sequences
	asCArray<char> str((int)cstr.GetLength());

	for( asUINT n = 0; n < cstr.GetLength(); n++ )
	{
		if( cstr[n] == '\\' )
		{
			++n;
			if( n == cstr.GetLength() ) return;

			if( cstr[n] == '"' )
				str.PushLast('"');
			else if( cstr[n] == '\'' )
				str.PushLast('\'');
			else if( cstr[n] == 'n' )
				str.PushLast('\n');
			else if( cstr[n] == 'r' )
				str.PushLast('\r');
			else if( cstr[n] == 't' )
				str.PushLast('\t');
			else if( cstr[n] == '0' )
				str.PushLast('\0');
			else if( cstr[n] == '\\' )
				str.PushLast('\\');
			else if( cstr[n] == 'x' || cstr[n] == 'X' )
			{
				++n;
				if( n == cstr.GetLength() ) break;

				int val = 0;
				if( cstr[n] >= '0' && cstr[n] <= '9' )
					val = cstr[n] - '0';
				else if( cstr[n] >= 'a' && cstr[n] <= 'f' )
					val = cstr[n] - 'a' + 10;
				else if( cstr[n] >= 'A' && cstr[n] <= 'F' )
					val = cstr[n] - 'A' + 10;
				else
					continue;

				++n;
				if( n == cstr.GetLength() )
				{
					str.PushLast((char)val);
					break;
				}

				if( cstr[n] >= '0' && cstr[n] <= '9' )
					val = val*16 + cstr[n] - '0';
				else if( cstr[n] >= 'a' && cstr[n] <= 'f' )
					val = val*16 + cstr[n] - 'a' + 10;
				else if( cstr[n] >= 'A' && cstr[n] <= 'F' )
					val = val*16 + cstr[n] - 'A' + 10;
				else
				{
					str.PushLast((char)val);
					continue;
				}

				str.PushLast((char)val);
			}
			else
				continue;
		}
		else
			str.PushLast(cstr[n]);
	}

	cstr.Assign(str.AddressOf(), str.GetLength());
}

void asCCompiler::ProcessHeredocStringConstant(asCString &str)
{
	// Remove first line if it only contains whitespace
	asUINT start;
	for( start = 0; start < str.GetLength(); start++ )
	{
		if( str[start] == '\n' )
		{
			// Remove the linebreak as well
			start++;
			break;
		}

		if( str[start] != ' '  &&
			str[start] != '\t' &&
			str[start] != '\r' )
		{
			// Don't remove anything
			start = 0;
			break;
		}
	}

	// Remove last line break and the line after that if it only contains whitespaces
	int end;
	for( end = (int)str.GetLength() - 1; end >= 0; end-- )
	{
		if( str[end] == '\n' )
			break;

		if( str[end] != ' '  &&
			str[end] != '\t' &&
			str[end] != '\r' )
		{
			// Don't remove anything
			end = (int)str.GetLength();
			break;
		}
	}

	if( end < 0 ) end = 0;

	asCString tmp;
	tmp.Assign(&str[start], end-start);

	str = tmp;
}

void asCCompiler::CompileConversion(asCScriptNode *node, asSExprContext *ctx)
{
	asSExprContext expr(engine);
	asCDataType to;
	if( node->nodeType == snConstructCall )
	{
		// Verify that there is only one argument
		if( node->lastChild->firstChild != node->lastChild->lastChild )
		{
			Error(TXT_ONLY_ONE_ARGUMENT_IN_CAST, node->lastChild);
		}

		// Compile the expression
		CompileAssignment(node->lastChild->firstChild, &expr);

		// Determine the requested type
		to = builder->CreateDataTypeFromNode(node->firstChild, script);
		to.MakeReadOnly(true); // Default to const
		assert(to.IsPrimitive());
	}
	else
	{
		// Compile the expression
		CompileAssignment(node->lastChild, &expr);

		// Determine the requested type
		to = builder->CreateDataTypeFromNode(node->firstChild, script);
		to = builder->ModifyDataTypeFromNode(to, node->firstChild->next, script, 0, 0);

		// If the type support object handles, then use it
		if( to.SupportHandles() )
		{
			to.MakeHandle(true);
			to.MakeReference(true);
		}
	}

	// Convert any reference to a variable
	if( expr.type.dataType.IsReference() ) ConvertToVariable(&expr);

	ImplicitConversion(&expr, to, node, true);

	IsVariableInitialized(&expr.type, node);

	// If no type conversion is really tried ignore it
	if( to == expr.type.dataType )
	{
		// This will keep information about constant type
		MergeExprContexts(ctx, &expr);
		ctx->type = expr.type;
		return;
	}

	if( to.IsEqualExceptConst(expr.type.dataType) && to.IsPrimitive() )
	{
		MergeExprContexts(ctx, &expr);
		ctx->type = expr.type;
		ctx->type.dataType.MakeReadOnly(true);
		return;
	}

	// The implicit conversion already does most of the conversions permitted,
	// here we'll only treat those conversions that require an explicit cast.

	bool conversionOK = false;
	if( !expr.type.isConstant )
	{
		if( !expr.type.dataType.IsObject() )
			ConvertToTempVariable(&expr);
		MergeExprContexts(ctx, &expr);

		if( to.IsObjectHandle() &&
			expr.type.dataType.IsObjectHandle() &&
			(expr.type.dataType.GetObjectType()->flags & asOBJ_SCRIPT_STRUCT) &&
			!(!to.IsHandleToConst() && expr.type.dataType.IsHandleToConst()) )
		{
			// Allow dynamic cast between object handles (only for script objects).
			// At run time this may result in a null handle,   
			// which when used will throw an exception
			conversionOK = true;
			ctx->bc.InstrDWORD(BC_Cast, engine->GetTypeIdFromDataType(to));

			ctx->type = expr.type;
			ctx->type.dataType = to;
		}
	}

	if( conversionOK )
		return;

	// Conversion not available
	ctx->type.SetDummy();

	asCString strTo, strFrom;

	strTo = to.Format();
	strFrom = expr.type.dataType.Format();

	asCString msg;
	msg.Format(TXT_NO_CONVERSION_s_TO_s, strFrom.AddressOf(), strTo.AddressOf());

	Error(msg.AddressOf(), node);
}

void asCCompiler::AfterFunctionCall(int funcID, asCArray<asSExprContext*> &args, asSExprContext *ctx, bool deferAll)
{
	asCScriptFunction *descr = builder->GetFunctionDescription(funcID);

	// Parameters that are sent by reference should be assigned
	// to the evaluated expression if it is an lvalue

	// Evaluate the arguments from last to first
	int n = (int)descr->parameterTypes.GetLength() - 1;
	for( ; n >= 0; n-- )
	{
		if( (descr->parameterTypes[n].IsReference() && (descr->inOutFlags[n] & 2)) ||
		    (descr->parameterTypes[n].IsObject() && deferAll) )
		{
			assert( !(descr->parameterTypes[n].IsReference() && (descr->inOutFlags[n] == 2)) || args[n]->origExpr );

			// For &inout, only store the argument if it is for a temporary variable
			if( engine->allowUnsafeReferences || 
				descr->inOutFlags[n] != 3 || args[n]->type.isTemporary )
			{
				// Store the argument for later processing
				asSDeferredParam outParam;
				outParam.argNode = args[n]->exprNode;
				outParam.argType = args[n]->type;
				outParam.argInOutFlags = descr->inOutFlags[n];
				outParam.origExpr = args[n]->origExpr;

				ctx->deferredParams.PushLast(outParam);
			}
		}
		else
		{
			// Release the temporary variable now
			ReleaseTemporaryVariable(args[n]->type, &ctx->bc);
		}
	}
}

void asCCompiler::ProcessDeferredParams(asSExprContext *ctx)
{
	if( isProcessingDeferredParams ) return;

	isProcessingDeferredParams = true;

	for( asUINT n = 0; n < ctx->deferredParams.GetLength(); n++ )
	{
		asSDeferredParam outParam = ctx->deferredParams[n];
		if( outParam.argInOutFlags < 2 ) // &in, or not reference
		{
			// Just release the variable
			ReleaseTemporaryVariable(outParam.argType, &ctx->bc);
		}
		else if( outParam.argInOutFlags == 2 ) // &out
		{
			asSExprContext *expr = outParam.origExpr;

			if( outParam.argType.dataType.IsObjectHandle() )
			{
				// Implicitly convert the value to a handle
				if( expr->type.dataType.IsObjectHandle() )
					expr->type.isExplicitHandle = true;
			}

			// Verify that the expression result in a lvalue
			if( IsLValue(expr->type) )
			{
				asSExprContext rctx(engine);
				rctx.type = outParam.argType;
				if( rctx.type.dataType.IsPrimitive() )
					rctx.type.dataType.MakeReference(false);
				else
				{
					rctx.bc.InstrSHORT(BC_PSF, outParam.argType.stackOffset);
					rctx.type.dataType.MakeReference(true);
				}

				asSExprContext o(engine);
				DoAssignment(&o, expr, &rctx, outParam.argNode, outParam.argNode, ttAssignment, outParam.argNode);

				if( !o.type.dataType.IsPrimitive() ) o.bc.Pop(PTR_SIZE);

				MergeExprContexts(ctx, &o);
			}
			else
			{
				// We must still evaluate the expression
				MergeExprContexts(ctx, expr);
				ctx->bc.Pop(expr->type.dataType.GetSizeOnStackDWords());

				// Give a warning
				Warning(TXT_ARG_NOT_LVALUE, outParam.argNode);

				ReleaseTemporaryVariable(outParam.argType, &ctx->bc);
			}

			ReleaseTemporaryVariable(expr->type, &ctx->bc);

			// Delete the original expression context
			DELETE(expr,asSExprContext);
		}
		else // &inout
		{
			if( outParam.argType.isTemporary )
				ReleaseTemporaryVariable(outParam.argType, &ctx->bc);
			else if( !outParam.argType.isVariable )
			{
				if( outParam.argType.dataType.IsObject() &&
					outParam.argType.dataType.GetBehaviour()->addref &&
					outParam.argType.dataType.GetBehaviour()->release )
				{
					// Release the object handle that was taken to guarantee the reference
					ReleaseTemporaryVariable(outParam.argType, &ctx->bc);
				}
			}
		}
	}

	ctx->deferredParams.SetLength(0);
	isProcessingDeferredParams = false;
}


void asCCompiler::CompileConstructCall(asCScriptNode *node, asSExprContext *ctx)
{
	// The first node is a datatype node
	asCString name;
	asCTypeInfo tempObj;
	asCArray<int> funcs;

	// It is possible that the name is really a constructor
	asCDataType dt;
	dt = builder->CreateDataTypeFromNode(node->firstChild, script);
	if( dt.GetObjectType() == 0 )
	{
		// This is a cast to a primitive type
		CompileConversion(node, ctx);
		return;
	}
	else
	{
		asSTypeBehaviour *beh = dt.GetBehaviour();
		if( beh )
			funcs = beh->constructors;

		tempObj.dataType = dt;
		tempObj.stackOffset = (short)AllocateVariable(dt, true);
		tempObj.dataType.MakeReference(true);
		tempObj.isTemporary = true;
		tempObj.isVariable = true;

		// Push the address of the object on the stack
		ctx->bc.InstrSHORT(BC_VAR, tempObj.stackOffset);
	}

	if( globalExpression )
	{
		Error(TXT_FUNCTION_IN_GLOBAL_EXPR, node);

		// Output dummy code
		ctx->type.SetDummy();
		return;
	}

	// Compile the arguments
	asCArray<asSExprContext *> args;
	asCArray<asCTypeInfo> temporaryVariables;

	CompileArgumentList(node->lastChild, args, &tempObj.dataType);

	// Special case: Allow calling func(void) with a void expression.
	if( args.GetLength() == 1 && args[0]->type.dataType == asCDataType::CreatePrimitive(ttVoid, false) )
	{
		// Evaluate the expression before the function call
		MergeExprContexts(ctx, args[0]);
		DELETE(args[0],asSExprContext);
		args.SetLength(0);
	}

	// Special case: If this is an object constructor and there are no arguments use the default constructor.
	// If none has been registered, just allocate the variable and push it on the stack.
	if( args.GetLength() == 0 )
	{
		asSTypeBehaviour *beh = tempObj.dataType.GetBehaviour();
		if( beh && beh->construct == 0 )
		{
			// Call the default constructor
			ctx->type = tempObj;
			ctx->bc.InstrWORD(BC_GETREF, 0);
			DefaultConstructor(&ctx->bc, tempObj.dataType);

			// Push the reference on the stack
			ctx->bc.InstrSHORT(BC_PSF, tempObj.stackOffset);
			return;
		}
	}

	MatchFunctions(funcs, args, node, name.AddressOf(), false);

	if( funcs.GetLength() != 1 )
	{
		// The error was reported by MatchFunctions()

		// Dummy value
		ctx->type.SetDummy();
	}
	else
	{
		asCByteCode objBC(engine);

		PrepareFunctionCall(funcs[0], &ctx->bc, args);

		MoveArgsToStack(funcs[0], &ctx->bc, args, false);

		int offset = 0;
		asCScriptFunction *descr = builder->GetFunctionDescription(funcs[0]);
		for( asUINT n = 0; n < args.GetLength(); n++ )
			offset += descr->parameterTypes[n].GetSizeOnStackDWords();

		ctx->bc.InstrWORD(BC_GETREF, (asWORD)offset);

		PerformFunctionCall(funcs[0], ctx, true, &args, tempObj.dataType.GetObjectType());

		// The constructor doesn't return anything,
		// so we have to manually inform the type of
		// the return value
		ctx->type = tempObj;

		// Push the address of the object on the stack again
		ctx->bc.InstrSHORT(BC_PSF, tempObj.stackOffset);
	}

	// Cleanup
	for( asUINT n = 0; n < args.GetLength(); n++ )
		if( args[n] )
		{
			DELETE(args[n],asSExprContext);
		}
}


void asCCompiler::CompileFunctionCall(asCScriptNode *node, asSExprContext *ctx, asCObjectType *objectType, bool objIsConst)
{
	asCString name;
	asCTypeInfo tempObj;
	asCArray<int> funcs;
	asCScriptNode *nm = node->firstChild;
	name.Assign(&script->code[nm->tokenPos], nm->tokenLength);
	if( !builder->GetObjectType(name.AddressOf()) )
	{
		if( objectType )
			builder->GetObjectMethodDescriptions(name.AddressOf(), objectType, funcs, objIsConst);
		else
			builder->GetFunctionDescriptions(name.AddressOf(), funcs);
	}

	if( globalExpression )
	{
		Error(TXT_FUNCTION_IN_GLOBAL_EXPR, node);

		// Output dummy code
		ctx->type.SetDummy();
		return;
	}

	// Compile the arguments
	asCArray<asSExprContext *> args;
	asCArray<asCTypeInfo> temporaryVariables;

	CompileArgumentList(node->lastChild, args, 0);

	// Special case: Allow calling func(void) with a void expression.
	if( args.GetLength() == 1 && args[0]->type.dataType == asCDataType::CreatePrimitive(ttVoid, false) )
	{
		// Evaluate the expression before the function call
		MergeExprContexts(ctx, args[0]);
		DELETE(args[0],asSExprContext);
		args.SetLength(0);
	}

	MatchFunctions(funcs, args, node, name.AddressOf(), objIsConst);

	if( funcs.GetLength() != 1 )
	{
		// The error was reported by MatchFunctions()

		// Dummy value
		ctx->type.SetDummy();
	}
	else
	{
		asCByteCode objBC(engine);

		objBC.AddCode(&ctx->bc);

		PrepareFunctionCall(funcs[0], &ctx->bc, args);

		// Verify if any of the args variable offsets are used in the other code.
		// If they are exchange the offset for a new one
		asUINT n;
		for( n = 0; n < args.GetLength(); n++ )
		{
			if( args[n]->type.isTemporary && objBC.IsVarUsed(args[n]->type.stackOffset) )
			{
				// Release the current temporary variable
				ReleaseTemporaryVariable(args[n]->type, 0);

				asCArray<int> usedVars;
				objBC.GetVarsUsed(usedVars);
				ctx->bc.GetVarsUsed(usedVars);

				asCDataType dt = args[n]->type.dataType;
				dt.MakeReference(false);
				int newOffset = AllocateVariableNotIn(dt, true, &usedVars);

				ctx->bc.ExchangeVar(args[n]->type.stackOffset, newOffset);
				args[n]->type.stackOffset = (short)newOffset;
				args[n]->type.isTemporary = true;
				args[n]->type.isVariable = true;
			}
		}

		ctx->bc.AddCode(&objBC);

		MoveArgsToStack(funcs[0], &ctx->bc, args, objectType ? true : false);

		int offset = 0;
		for( n = 0; n < args.GetLength(); n++ )
			offset += args[n]->type.dataType.GetSizeOnStackDWords();

		PerformFunctionCall(funcs[0], ctx, false, &args, 0);
	}

	// Cleanup
	for( asUINT n = 0; n < args.GetLength(); n++ )
		if( args[n] )
		{
			DELETE(args[n],asSExprContext);
		}
}

void asCCompiler::CompileExpressionPreOp(asCScriptNode *node, asSExprContext *ctx)
{
	int op = node->tokenType;

	IsVariableInitialized(&ctx->type, node);

	if( op == ttHandle )
	{
		// Verify that the type allow its handle to be taken
		if( ctx->type.isExplicitHandle || !ctx->type.dataType.IsObject() || !ctx->type.dataType.GetObjectType()->beh.addref || !ctx->type.dataType.GetObjectType()->beh.release )
		{
			Error(TXT_OBJECT_HANDLE_NOT_SUPPORTED, node);
		}

		// Objects that are not local variables are not references
		if( !ctx->type.dataType.IsReference() && !(ctx->type.dataType.IsObject() && !ctx->type.isVariable) )
		{
			Error(TXT_NOT_VALID_REFERENCE, node);
		}

		// If this is really an object then the handle created is a const handle
		bool makeConst = !ctx->type.dataType.IsObjectHandle();

		// Mark the type as an object handle
		ctx->type.dataType.MakeHandle(true);
		ctx->type.isExplicitHandle = true;
		if( makeConst )
			ctx->type.dataType.MakeReadOnly(true);
	}
	else if( op == ttMinus && ctx->type.dataType.IsObject() )
	{
		asCTypeInfo objType = ctx->type;

		Dereference(ctx, true);

		// Check if the variable is initialized (if it indeed is a variable)
		IsVariableInitialized(&ctx->type, node);

		// Now find a matching function for the object type
		asSTypeBehaviour *beh = ctx->type.dataType.GetBehaviour();
		if( beh == 0 )
		{
			asCString str;
			str.Format(TXT_OBJECT_DOESNT_SUPPORT_NEGATE_OP);
			Error(str.AddressOf(), node);
		}
		else
		{
			// Find the negate operator
			int opNegate = 0;
			bool found = false;
			asUINT n;
			for( n = 0; n < beh->operators.GetLength(); n+= 2 )
			{
				// Only accept the negate operator
				if( ttMinus == beh->operators[n] &&
					engine->scriptFunctions[beh->operators[n+1]]->parameterTypes.GetLength() == 0 )
				{
					found = true;
					opNegate = beh->operators[n+1];
					break;
				}
			}

			// Did we find a suitable function?
			if( found )
			{
				PerformFunctionCall(opNegate, ctx);
			}
			else
			{
				asCString str;
				str.Format(TXT_OBJECT_DOESNT_SUPPORT_NEGATE_OP);
				Error(str.AddressOf(), node);
			}
		}

		// Release the potentially temporary object
		ReleaseTemporaryVariable(objType, &ctx->bc);
	}
	else if( op == ttPlus || op == ttMinus )
	{
		asCDataType to = ctx->type.dataType;

		// TODO: The case -2147483648 gives an unecessary warning of changed sign for implicit conversion

		if( ctx->type.dataType.IsUnsignedType() )
		{
			if( ctx->type.dataType.GetSizeInMemoryBytes() == 1 )
				to.SetTokenType(ttInt8);
			else if( ctx->type.dataType.GetSizeInMemoryBytes() == 2 )
				to.SetTokenType(ttInt16);
			else if( ctx->type.dataType.GetSizeInMemoryBytes() == 4 )
				to.SetTokenType(ttInt);
			else if( ctx->type.dataType.GetSizeInMemoryBytes() == 8 )
				to.SetTokenType(ttInt64);
			else
				Error(TXT_INVALID_TYPE, node);
		}

		if( ctx->type.dataType.IsReference() ) ConvertToVariable(ctx);
		ImplicitConversion(ctx, to, node, false);

		if( !ctx->type.isConstant )
		{
			ConvertToTempVariable(ctx);

			if( op == ttMinus )
			{
				if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
					ctx->bc.InstrSHORT(BC_NEGi, ctx->type.stackOffset);
				else if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
					ctx->bc.InstrSHORT(BC_NEGi64, ctx->type.stackOffset);
				else if( ctx->type.dataType.IsFloatType() )
					ctx->bc.InstrSHORT(BC_NEGf, ctx->type.stackOffset);
				else if( ctx->type.dataType.IsDoubleType() )
					ctx->bc.InstrSHORT(BC_NEGd, ctx->type.stackOffset);
				else
					Error(TXT_ILLEGAL_OPERATION, node);

				return;
			}
		}
		else
		{
			if( op == ttMinus )
			{
				if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
					ctx->type.intValue = -ctx->type.intValue;
				else if( ctx->type.dataType.IsIntegerType() && ctx->type.dataType.GetSizeInMemoryDWords() == 2 )
					ctx->type.qwordValue = -(asINT64)ctx->type.qwordValue;
				else if( ctx->type.dataType.IsFloatType() )
					ctx->type.floatValue = -ctx->type.floatValue;
				else if( ctx->type.dataType.IsDoubleType() )
					ctx->type.doubleValue = -ctx->type.doubleValue;
				else
					Error(TXT_ILLEGAL_OPERATION, node);

				return;
			}
		}

		if( op == ttPlus )
		{
			if( !ctx->type.dataType.IsIntegerType() &&
				!ctx->type.dataType.IsFloatType() &&
				!ctx->type.dataType.IsDoubleType() )
				Error(TXT_ILLEGAL_OPERATION, node);
		}
	}
	else if( op == ttNot )
	{
		if( ctx->type.dataType.IsEqualExceptRefAndConst(asCDataType::CreatePrimitive(ttBool, true)) )
		{
			if( ctx->type.isConstant )
			{
				ctx->type.dwordValue = (ctx->type.dwordValue == 0 ? VALUE_OF_BOOLEAN_TRUE : 0);
				return;
			}

			ConvertToTempVariable(ctx);

			ctx->bc.InstrSHORT(BC_NOT, ctx->type.stackOffset);
		}
		else
			Error(TXT_ILLEGAL_OPERATION, node);
	}
	else if( op == ttBitNot )
	{
		asCDataType to = ctx->type.dataType;

		if( ctx->type.dataType.IsIntegerType() )
		{
			if( ctx->type.dataType.GetSizeInMemoryBytes() == 1 )
				to.SetTokenType(ttUInt8);
			else if( ctx->type.dataType.GetSizeInMemoryBytes() == 2 )
				to.SetTokenType(ttUInt16);
			else if( ctx->type.dataType.GetSizeInMemoryBytes() == 4 )
				to.SetTokenType(ttUInt);
			else if( ctx->type.dataType.GetSizeInMemoryBytes() == 8 )
				to.SetTokenType(ttUInt64);
			else
				Error(TXT_INVALID_TYPE, node);
		}

		if( ctx->type.dataType.IsReference() ) ConvertToVariable(ctx);
		ImplicitConversion(ctx, to, node, false);

		if( ctx->type.dataType.IsUnsignedType() )
		{
			if( ctx->type.isConstant )
			{
				ctx->type.qwordValue = ~ctx->type.qwordValue;
				return;
			}

			ConvertToTempVariable(ctx);
			if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				ctx->bc.InstrSHORT(BC_BNOT, ctx->type.stackOffset);
			else
				ctx->bc.InstrSHORT(BC_BNOT64, ctx->type.stackOffset);
		}
		else
			Error(TXT_ILLEGAL_OPERATION, node);
	}
	else if( op == ttInc || op == ttDec )
	{
		// Need a reference to the primitive that will be updated
		// The result of this expression is the same reference as before
		if( globalExpression )
			Error(TXT_INC_OP_IN_GLOBAL_EXPR, node);

		// Make sure the reference isn't a temporary variable
		if( ctx->type.isTemporary )
			Error(TXT_REF_IS_TEMP, node);
		if( ctx->type.dataType.IsReadOnly() )
			Error(TXT_REF_IS_READ_ONLY, node);

		if( ctx->type.isVariable )
			ConvertToReference(ctx);
		else if( !ctx->type.dataType.IsReference() )
			Error(TXT_NOT_VALID_REFERENCE, node);

		if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt64, false)) ||
			ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt64, false)) )
		{
			if( op == ttInc )
				ctx->bc.Instr(BC_INCi64);
			else
				ctx->bc.Instr(BC_DECi64);
		}
		else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt, false)) ||
			ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt, false)) )
		{
			if( op == ttInc )
				ctx->bc.Instr(BC_INCi);
			else
				ctx->bc.Instr(BC_DECi);
		}
		else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt16, false)) ||
			     ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt16, false)) )
		{
			if( op == ttInc )
				ctx->bc.Instr(BC_INCi16);
			else
				ctx->bc.Instr(BC_DECi16);
		}
		else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt8, false)) ||
			     ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt8, false)) )
		{
			if( op == ttInc )
				ctx->bc.Instr(BC_INCi8);
			else
				ctx->bc.Instr(BC_DECi8);
		}
		else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttFloat, false)) )
		{
			if( op == ttInc )
				ctx->bc.Instr(BC_INCf);
			else
				ctx->bc.Instr(BC_DECf);
		}
		else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttDouble, false)) )
		{
			if( op == ttInc )
				ctx->bc.Instr(BC_INCd);
			else
				ctx->bc.Instr(BC_DECd);
		}
		else
			Error(TXT_ILLEGAL_OPERATION, node);
	}
	else
		// Unknown operator
		assert(false);
}

void asCCompiler::ConvertToReference(asSExprContext *ctx)
{
	if( ctx->type.isVariable )
	{
		ctx->bc.InstrSHORT(BC_LDV, ctx->type.stackOffset);
		ctx->type.dataType.MakeReference(true);
		ctx->type.Set(ctx->type.dataType);
	}
}

void asCCompiler::CompileExpressionPostOp(asCScriptNode *node, asSExprContext *ctx)
{
	int op = node->tokenType;

	if( op == ttInc || op == ttDec )
	{
		if( globalExpression )
			Error(TXT_INC_OP_IN_GLOBAL_EXPR, node);

		// Make sure the reference isn't a temporary variable
		if( ctx->type.isTemporary )
			Error(TXT_REF_IS_TEMP, node);
		if( ctx->type.dataType.IsReadOnly() )
			Error(TXT_REF_IS_READ_ONLY, node);

		if( ctx->type.isVariable )
			ConvertToReference(ctx);
		else if( !ctx->type.dataType.IsReference() )
			Error(TXT_NOT_VALID_REFERENCE, node);

		// Copy the value to a temp before changing it
		ConvertToTempVariable(ctx);

		// Increment the value pointed to by the reference still in the register
		bcInstr iInc = BC_INCi, iDec = BC_DECi;
		if( ctx->type.dataType.IsDoubleType() )
		{
			iInc = BC_INCd;
			iDec = BC_DECd;
		}
		else if( ctx->type.dataType.IsFloatType() )
		{
			iInc = BC_INCf;
			iDec = BC_DECf;
		}
		else if( ctx->type.dataType.IsIntegerType() || ctx->type.dataType.IsUnsignedType() )
		{
			if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt16, false)) ||
					 ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt16, false)) )
			{
				iInc = BC_INCi16;
				iDec = BC_DECi16;
			}
			else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt8, false)) ||
					 ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt8, false)) )
			{
				iInc = BC_INCi8;
				iDec = BC_DECi8;
			}
			else if( ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttInt64, false)) ||
					 ctx->type.dataType.IsEqualExceptRef(asCDataType::CreatePrimitive(ttUInt64, false)) )
			{
				iInc = BC_INCi64;
				iDec = BC_DECi64;
			}
		}
		else
		{
			Error(TXT_ILLEGAL_OPERATION, node);
		}

		if( op == ttInc ) ctx->bc.Instr(iInc); else ctx->bc.Instr(iDec);
	}
	else if( op == ttDot )
	{
		// Check if the variable is initialized (if it indeed is a variable)
		IsVariableInitialized(&ctx->type, node);

		if( node->firstChild->nodeType == snIdentifier )
		{
			// Get the property name
			GETSTRING(name, &script->code[node->firstChild->tokenPos], node->firstChild->tokenLength);

			if( !ctx->type.dataType.IsPrimitive() )
				Dereference(ctx, true);

			if( ctx->type.dataType.IsObjectHandle() )
			{
				// Convert the handle to a normal object
				asCDataType dt = ctx->type.dataType;
				dt.MakeHandle(false);

				ImplicitConversion(ctx, dt, node, false);
			}

			// Find the property offset and type
			if( ctx->type.dataType.IsObject() )
			{
				bool isConst = ctx->type.dataType.IsReadOnly();

				asCProperty *prop = builder->GetObjectProperty(ctx->type.dataType, name.AddressOf());
				if( prop )
				{
					// Put the offset on the stack
					ctx->bc.InstrINT(BC_ADDSi, prop->byteOffset);

					if( prop->type.IsReference() )
						ctx->bc.Instr(BC_RDSPTR);

					// Reference to primitive must be stored in the temp register
					if( prop->type.IsPrimitive() )
					{
						// The ADD offset command should store the reference in the register directly
						ctx->bc.Instr(BC_PopRPtr);
					}

					// Set the new type (keeping info about temp variable)
					ctx->type.dataType = prop->type;
					ctx->type.dataType.MakeReference(true);
					ctx->type.isVariable = false;

					if( ctx->type.dataType.IsObject() && !ctx->type.dataType.IsObjectHandle() )
					{
						// Objects that are members are not references
						ctx->type.dataType.MakeReference(false);
					}

					ctx->type.dataType.MakeReadOnly(isConst ? true : prop->type.IsReadOnly());
				}
				else
				{
					asCString str;
					str.Format(TXT_s_NOT_MEMBER_OF_s, name.AddressOf(), ctx->type.dataType.Format().AddressOf());
					Error(str.AddressOf(), node);
				}
			}
			else
			{
				asCString str;
				str.Format(TXT_s_NOT_MEMBER_OF_s, name.AddressOf(), ctx->type.dataType.Format().AddressOf());
				Error(str.AddressOf(), node);
			}
		}
		else
		{
			if( globalExpression )
				Error(TXT_METHOD_IN_GLOBAL_EXPR, node);

			// Make sure it is an object we are accessing
			if( !ctx->type.dataType.IsObject() )
			{
				asCString str;
				str.Format(TXT_ILLEGAL_OPERATION_ON_s, ctx->type.dataType.Format().AddressOf());
				Error(str.AddressOf(), node);
			}
			else
			{
				// We need the object pointer
				Dereference(ctx, true);

				bool isConst = ctx->type.dataType.IsReadOnly();

				asCObjectType *trueObj = ctx->type.dataType.GetObjectType();

				asCTypeInfo objType = ctx->type;

				// Compile function call
				CompileFunctionCall(node->firstChild, ctx, trueObj, isConst);

				// Release the potentially temporary object
				ReleaseTemporaryVariable(objType, &ctx->bc);
			}
		}
	}
	else if( op == ttOpenBracket )
	{
		Dereference(ctx, true);
		bool isConst = ctx->type.dataType.IsReadOnly();

		if( ctx->type.dataType.IsObjectHandle() )
		{
			// Convert the handle to a normal object
			asCDataType dt = ctx->type.dataType;
			dt.MakeHandle(false);

			ImplicitConversion(ctx, dt, node, false);
		}

		// Check if the variable is initialized (if it indeed is a variable)
		IsVariableInitialized(&ctx->type, node);

		// Compile the expression
		asSExprContext expr(engine);
		CompileAssignment(node->firstChild, &expr);

		asCTypeInfo objType = ctx->type;

		// Now find a matching function for the object type and indexing type
		asSTypeBehaviour *beh = ctx->type.dataType.GetBehaviour();
		if( beh == 0 )
		{
			asCString str;
			str.Format(TXT_OBJECT_DOESNT_SUPPORT_INDEX_OP, ctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), node);
		}
		else
		{
			asCArray<int> ops;
			asUINT n;
			if( isConst )
			{
				// Only list const behaviours
				for( n = 0; n < beh->operators.GetLength(); n += 2 )
				{
					if( ttOpenBracket == beh->operators[n] && engine->scriptFunctions[beh->operators[n+1]]->isReadOnly )
						ops.PushLast(beh->operators[n+1]);
				}
			}
			else
			{
				// TODO:
				// Prefer non-const over const
				for( n = 0; n < beh->operators.GetLength(); n += 2 )
				{
					if( ttOpenBracket == beh->operators[n] )
						ops.PushLast(beh->operators[n+1]);
				}
			}

			asCArray<int> ops1;
			MatchArgument(ops, ops1, &expr.type, 0);

			if( !isConst )
				FilterConst(ops1);

			// Did we find a suitable function?
			if( ops1.GetLength() == 1 )
			{
				asCScriptFunction *descr = engine->scriptFunctions[ops1[0]];

				// Store the code for the object
				asCByteCode objBC(engine);
				objBC.AddCode(&ctx->bc);

				// Add code for arguments

				PrepareArgument(&descr->parameterTypes[0], &expr, node->firstChild, true, descr->inOutFlags[0]);
				MergeExprContexts(ctx, &expr);

				if( descr->parameterTypes[0].IsReference() )
				{
					if( descr->parameterTypes[0].IsObject() && !descr->parameterTypes[0].IsObjectHandle() )
						ctx->bc.InstrWORD(BC_GETOBJREF, 0);
					else
						ctx->bc.InstrWORD(BC_GETREF, 0);
				}
				else if( descr->parameterTypes[0].IsObject() )
				{
					ctx->bc.InstrWORD(BC_GETOBJ, 0);

					// The temporary variable must not be freed as it will no longer hold an object
					DeallocateVariable(expr.type.stackOffset);
					expr.type.isTemporary = false;
				}

				// Add the code for the object again
				ctx->bc.AddCode(&objBC);

				asCArray<asSExprContext*> args;
				args.PushLast(&expr);
				PerformFunctionCall(descr->id, ctx, false, &args);

				// TODO: Ugly code
				// The default array returns a reference to the subtype
				if( objType.dataType.IsScriptArray() )
				{
					ctx->type.dataType = objType.dataType.GetSubType();
					if( !ctx->type.dataType.IsPrimitive() )
					{
						// The reference is currently stored in the register
						ctx->bc.Instr(BC_PshRPtr);
					}
					if( !ctx->type.dataType.IsObject() || ctx->type.dataType.IsObjectHandle() )
						ctx->type.dataType.MakeReference(true);
					if( isConst )
						ctx->type.dataType.MakeReadOnly(true);
				}
			}
			else if( ops.GetLength() > 1 )
			{
				Error(TXT_MORE_THAN_ONE_MATCHING_OP, node);
			}
			else
			{
				asCString str;
				str.Format(TXT_NO_MATCHING_OP_FOUND_FOR_TYPE_s, expr.type.dataType.Format().AddressOf());
				Error(str.AddressOf(), node);
			}
		}

		// Release the potentially temporary object
		ReleaseTemporaryVariable(objType, &ctx->bc);
	}
}

int asCCompiler::GetPrecedence(asCScriptNode *op)
{
	// x*y, x/y, x%y
	// x+y, x-y
	// x<=y, x<y, x>=y, x>y
	// x==y, x!=y
	// x and y
	// (x xor y)
	// x or y

	// The following are not used in this function,
	// but should have lower precedence than the above
	// x ? y : z
	// x = y

	// The expression term have the highest precedence
	if( op->nodeType == snExprTerm )
		return 1;

	// Evaluate operators by token
	int tokenType = op->tokenType;
	if( tokenType == ttStar || tokenType == ttSlash || tokenType == ttPercent )
		return 0;

	if( tokenType == ttPlus || tokenType == ttMinus )
		return -1;

	if( tokenType == ttBitShiftLeft ||
		tokenType == ttBitShiftRight ||
		tokenType == ttBitShiftRightArith )
		return -2;

	if( tokenType == ttAmp )
		return -3;

	if( tokenType == ttBitXor )
		return -4;

	if( tokenType == ttBitOr )
		return -5;

	if( tokenType == ttLessThanOrEqual ||
		tokenType == ttLessThan ||
		tokenType == ttGreaterThanOrEqual ||
		tokenType == ttGreaterThan )
		return -6;

	if( tokenType == ttEqual || tokenType == ttNotEqual || tokenType == ttXor )
		return -7;

	if( tokenType == ttAnd )
		return -8;

	if( tokenType == ttOr )
		return -9;

	// Unknown operator
	assert(false);

	return 0;
}

int asCCompiler::MatchArgument(asCArray<int> &funcs, asCArray<int> &matches, const asCTypeInfo *argType, int paramNum, bool allowObjectConstruct)
{
	bool isExactMatch        = false;
	bool isMatchExceptConst  = false;
	bool isMatchWithBaseType = false;
	bool isMatchExceptSign   = false;
	bool isMatchNotVarType   = false;

	asUINT n;

	matches.SetLength(0);

	for( n = 0; n < funcs.GetLength(); n++ )
	{
		asCScriptFunction *desc = builder->GetFunctionDescription(funcs[n]);

		// Does the function have arguments enough?
		if( (int)desc->parameterTypes.GetLength() <= paramNum )
			continue;

		// Can we make the match by implicit conversion?
		asSExprContext ti(engine);
		ti.type = *argType;
		if( argType->dataType.IsPrimitive() ) ti.type.dataType.MakeReference(false);
		ImplicitConversion(&ti, desc->parameterTypes[paramNum], 0, false, false, 0, allowObjectConstruct);
		if( desc->parameterTypes[paramNum].IsEqualExceptRef(ti.type.dataType) )
		{
			// Is it an exact match?
			if( argType->dataType.IsEqualExceptRef(ti.type.dataType) )
			{
				if( !isExactMatch ) matches.SetLength(0);

				isExactMatch = true;

				matches.PushLast(funcs[n]);
				continue;
			}

			if( !isExactMatch )
			{
				// Is it a match except const?
				if( argType->dataType.IsEqualExceptRefAndConst(ti.type.dataType) )
				{
					if( !isMatchExceptConst ) matches.SetLength(0);

					isMatchExceptConst = true;

					matches.PushLast(funcs[n]);
					continue;
				}

				if( !isMatchExceptConst )
				{
					// Is it a size promotion, e.g. int8 -> int?
					if( argType->dataType.IsSamePrimitiveBaseType(ti.type.dataType) )
					{
						if( !isMatchWithBaseType ) matches.SetLength(0);

						isMatchWithBaseType = true;

						matches.PushLast(funcs[n]);
						continue;
					}

					if( !isMatchWithBaseType )
					{
						// Conversion between signed and unsigned integer is better than between integer and float

						// Is it a match except for sign?
						if( argType->dataType.IsIntegerType() && ti.type.dataType.IsUnsignedType() ||
							argType->dataType.IsUnsignedType() && ti.type.dataType.IsIntegerType() )
						{
							if( !isMatchExceptSign ) matches.SetLength(0);

							isMatchExceptSign = true;

							matches.PushLast(funcs[n]);
							continue;
						}

						if( !isMatchExceptSign )
						{
							// If there was any match without a var type it has higher priority
							if( desc->parameterTypes[paramNum].GetTokenType() != ttQuestion )
							{
								if( !isMatchNotVarType ) matches.SetLength(0);

								isMatchNotVarType = true;

								matches.PushLast(funcs[n]);
								continue;
							}

							// Implicit conversion to ?& has the smallest priority
							if( !isMatchNotVarType )
								matches.PushLast(funcs[n]);
						}
					}
				}
			}
		}
	}

	return (int)matches.GetLength();
}

void asCCompiler::PrepareArgument2(asSExprContext *ctx, asSExprContext *arg, asCDataType *paramType, bool isFunction, int refType, asCArray<int> *reservedVars)
{
	asSExprContext e(engine);

	// Reference parameters whose value won't be used don't evaluate the expression
	if( !paramType->IsReference() || (refType & 1) )
	{
		MergeExprContexts(&e, arg);
	}
	else
	{
		// Store the original bytecode so that it can be reused when processing the deferred output parameter
		asSExprContext *orig = NEW(asSExprContext)(engine);
		MergeExprContexts(orig, arg);
		orig->exprNode = arg->exprNode;
		orig->type = arg->type;

		arg->origExpr = orig;
	}

	e.type = arg->type;
	PrepareArgument(paramType, &e, arg->exprNode, isFunction, refType, reservedVars);
	arg->type = e.type;
	ctx->bc.AddCode(&e.bc);
}

bool asCCompiler::CompileOverloadedOperator(asCScriptNode *node, asSExprContext *lctx, asSExprContext *rctx, asSExprContext *ctx)
{
	// TODO: An operator can be overloaded for an object or can be global

	// What type of operator is it?
	int token = node->tokenType;

	// TODO: Only search in config groups to which the module has access
	// What overloaded operators of this type do we have?
	asCArray<int> ops;
	asUINT n;
	for( n = 0; n < engine->globalBehaviours.operators.GetLength(); n += 2 )
	{
		if( token == engine->globalBehaviours.operators[n] )
		{
			int funcId = engine->globalBehaviours.operators[n+1];

			// Find the config group for the global function
			asCConfigGroup *group = engine->FindConfigGroupForFunction(funcId);
			if( !group || group->HasModuleAccess(builder->module->name.AddressOf()) )
				ops.PushLast(funcId);
		}
	}

	// Find the best matches for each argument
	asCArray<int> ops1;
	asCArray<int> ops2;
	MatchArgument(ops, ops1, &lctx->type, 0);
	MatchArgument(ops, ops2, &rctx->type, 1);

	// Find intersection of the two sets of matching operators
	ops.SetLength(0);
	for( n = 0; n < ops1.GetLength(); n++ )
	{
		for( asUINT m = 0; m < ops2.GetLength(); m++ )
		{
			if( ops1[n] == ops2[m] )
			{
				ops.PushLast(ops1[n]);
				break;
			}
		}
	}

	// Did we find an operator?
	if( ops.GetLength() == 1 )
	{
		asCScriptFunction *descr = engine->scriptFunctions[ops[0]];

		// Add code for arguments
		asCArray<int> reserved;
		rctx->bc.GetVarsUsed(reserved);

		PrepareArgument2(ctx, lctx, &descr->parameterTypes[0], true, descr->inOutFlags[0], &reserved);
		PrepareArgument2(ctx, rctx, &descr->parameterTypes[1], true, descr->inOutFlags[1]);

		// Swap the order of the arguments
		if( lctx->type.dataType.GetSizeOnStackDWords() == 2 && rctx->type.dataType.GetSizeOnStackDWords() == 2 )
			ctx->bc.Instr(BC_SWAP8);
		else if( lctx->type.dataType.GetSizeOnStackDWords() == 2 )
			ctx->bc.Instr(BC_SWAP48);
		else if( rctx->type.dataType.GetSizeOnStackDWords() == 2 )
			ctx->bc.Instr(BC_SWAP84);
		else
			ctx->bc.Instr(BC_SWAP4);

		asCArray<asSExprContext*> args(2);
		args.PushLast(lctx);
		args.PushLast(rctx);

		MoveArgsToStack(descr->id, &ctx->bc, args, false);

		PerformFunctionCall(descr->id, ctx, false, &args);

		// Don't continue
		return true;
	}
	else if( ops.GetLength() > 1 )
	{
		Error(TXT_MORE_THAN_ONE_MATCHING_OP, node);

		// Don't continue
		return true;
	}

	// No suitable operator was found
	return false;
}

void asCCompiler::CompileOperator(asCScriptNode *node, asSExprContext *lctx, asSExprContext *rctx, asSExprContext *ctx)
{
	IsVariableInitialized(&lctx->type, node);
	IsVariableInitialized(&rctx->type, node);

	if( lctx->type.isExplicitHandle || rctx->type.isExplicitHandle )
	{
		CompileOperatorOnHandles(node, lctx, rctx, ctx);

		return;
	}
	else
	{
		// Compile an overloaded operator for the two operands
		if( CompileOverloadedOperator(node, lctx, rctx, ctx) )
			return;

		// If either of the types are objects we shouldn't continue
		if( lctx->type.dataType.IsObject() || rctx->type.dataType.IsObject() )
		{
			asCString str;
			str.Format(TXT_NO_MATCHING_OP_FOUND_FOR_TYPES_s_AND_s, lctx->type.dataType.Format().AddressOf(), rctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), node);
			ctx->type.SetDummy();
			return;
		}

		// Make sure we have two variables or constants
		if( lctx->type.dataType.IsReference() ) ConvertToVariableNotIn(lctx, rctx);
		if( rctx->type.dataType.IsReference() ) ConvertToVariable(rctx);

		// Math operators
		// + - * / % += -= *= /= %=
		int op = node->tokenType;
		if( op == ttPlus    || op == ttAddAssign ||
			op == ttMinus   || op == ttSubAssign ||
			op == ttStar    || op == ttMulAssign ||
			op == ttSlash   || op == ttDivAssign ||
			op == ttPercent || op == ttModAssign )
		{
			CompileMathOperator(node, lctx, rctx, ctx);

			return;
		}

		// Bitwise operators
		// << >> >>> & | ^ <<= >>= >>>= &= |= ^=
		if( op == ttAmp                || op == ttAndAssign         ||
			op == ttBitOr              || op == ttOrAssign          ||
			op == ttBitXor             || op == ttXorAssign         ||
			op == ttBitShiftLeft       || op == ttShiftLeftAssign   ||
			op == ttBitShiftRight      || op == ttShiftRightLAssign ||
			op == ttBitShiftRightArith || op == ttShiftRightAAssign )
		{
			CompileBitwiseOperator(node, lctx, rctx, ctx);

			return;
		}

		// Comparison operators
		// == != < > <= >=
		if( op == ttEqual       || op == ttNotEqual           ||
			op == ttLessThan    || op == ttLessThanOrEqual    ||
			op == ttGreaterThan || op == ttGreaterThanOrEqual )
		{
			CompileComparisonOperator(node, lctx, rctx, ctx);

			return;
		}

		// Boolean operators
		// && || ^^
		if( op == ttAnd || op == ttOr || op == ttXor )
		{
			CompileBooleanOperator(node, lctx, rctx, ctx);

			return;
		}
	}

	assert(false);
}

void asCCompiler::ConvertToTempVariableNotIn(asSExprContext *ctx, asSExprContext *exclude)
{
	asCArray<int> excludeVars;
	if( exclude ) exclude->bc.GetVarsUsed(excludeVars);
	ConvertToTempVariableNotIn(ctx, &excludeVars);
}

void asCCompiler::ConvertToTempVariableNotIn(asSExprContext *ctx, asCArray<int> *reservedVars)
{
	ConvertToVariableNotIn(ctx, reservedVars);
	if( !ctx->type.isTemporary )
	{
		if( ctx->type.dataType.IsPrimitive() )
		{
			int offset = AllocateVariableNotIn(ctx->type.dataType, true, reservedVars);
			if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				ctx->bc.InstrW_W(BC_CpyVtoV4, offset, ctx->type.stackOffset);
			else
				ctx->bc.InstrW_W(BC_CpyVtoV8, offset, ctx->type.stackOffset);
			ctx->type.SetVariable(ctx->type.dataType, offset, true);
		}
		else if( ctx->type.dataType.IsObjectHandle() )
		{
			assert(false);
		}
		else // ctx->type.dataType.IsObject()
		{
			// Make sure the variable is not used in the expression
			asCArray<int> vars;
			ctx->bc.GetVarsUsed(vars);
			int offset = AllocateVariableNotIn(ctx->type.dataType, true, &vars);

			// Allocate and construct the temporary object
			asCByteCode tmpBC(engine);
			CompileConstructor(ctx->type.dataType, offset, &tmpBC);

			// Insert the code before the expression code
			tmpBC.AddCode(&ctx->bc);
			ctx->bc.AddCode(&tmpBC);

			// Assign the evaluated expression to the temporary variable
			PrepareForAssignment(&ctx->type.dataType, ctx, 0);

			asCTypeInfo type;
			type.SetVariable(ctx->type.dataType, offset, true);

			ctx->bc.InstrSHORT(BC_PSF, (short)offset);

			PerformAssignment(&type, &ctx->type, &ctx->bc, 0);

			ReleaseTemporaryVariable(ctx->type, &ctx->bc);

			ctx->type = type;
		}
	}
}

void asCCompiler::ConvertToTempVariable(asSExprContext *ctx)
{
	ConvertToVariable(ctx);
	if( !ctx->type.isTemporary )
	{
		if( ctx->type.dataType.IsPrimitive() )
		{
			int offset = AllocateVariable(ctx->type.dataType, true);
			if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				ctx->bc.InstrW_W(BC_CpyVtoV4, offset, ctx->type.stackOffset);
			else
				ctx->bc.InstrW_W(BC_CpyVtoV8, offset, ctx->type.stackOffset);
			ctx->type.SetVariable(ctx->type.dataType, offset, true);
		}
		else if( ctx->type.dataType.IsObjectHandle() )
		{
			assert(false);
		}
		else // ctx->type.dataType.IsObject()
		{
			// Make sure the variable is not used in the expression
			asCArray<int> vars;
			ctx->bc.GetVarsUsed(vars);
			int offset = AllocateVariableNotIn(ctx->type.dataType, true, &vars);

			// Allocate and construct the temporary object
			asCByteCode tmpBC(engine);
			CompileConstructor(ctx->type.dataType, offset, &tmpBC);

			// Insert the code before the expression code
			tmpBC.AddCode(&ctx->bc);
			ctx->bc.AddCode(&tmpBC);

			// Assign the evaluated expression to the temporary variable
			PrepareForAssignment(&ctx->type.dataType, ctx, 0);

			asCTypeInfo type;
			type.SetVariable(ctx->type.dataType, offset, true);

			ctx->bc.InstrSHORT(BC_PSF, (short)offset);

			PerformAssignment(&type, &ctx->type, &ctx->bc, 0);

			ReleaseTemporaryVariable(ctx->type, &ctx->bc);

			ctx->type = type;
		}
	}
}

void asCCompiler::ConvertToVariable(asSExprContext *ctx)
{
	ConvertToVariableNotIn(ctx, (asCArray<int>*)0);
}

void asCCompiler::ConvertToVariableNotIn(asSExprContext *ctx, asCArray<int> *reservedVars)
{
	if( !ctx->type.isVariable )
	{
		asCArray<int> excludeVars;
		if( reservedVars ) excludeVars.Concatenate(*reservedVars);
		int offset;
		if( ctx->type.dataType.IsObjectHandle() )
		{
			offset = AllocateVariableNotIn(ctx->type.dataType, true, &excludeVars);
			if( ctx->type.IsNullConstant() )
			{
				// TODO: Adapt pointer size
				ctx->bc.InstrSHORT_DW(BC_SetV4, (short)offset, 0);
			}
			else
			{
				// Copy the object handle to a variable
				ctx->bc.InstrSHORT(BC_PSF, (short)offset);
				ctx->bc.InstrPTR(BC_REFCPY, builder->module->RefObjectType(ctx->type.dataType.GetObjectType()));
				ctx->bc.Pop(PTR_SIZE);
			}

			ReleaseTemporaryVariable(ctx->type, &ctx->bc);
			ctx->type.SetVariable(ctx->type.dataType, offset, true);
		}
		else if( ctx->type.dataType.IsPrimitive() )
		{
			if( ctx->type.isConstant )
			{
				offset = AllocateVariableNotIn(ctx->type.dataType, true, &excludeVars);
				if( ctx->type.dataType.GetSizeInMemoryBytes() == 1 )
					ctx->bc.InstrSHORT_B(BC_SetV1, (short)offset, ctx->type.byteValue);
				else if( ctx->type.dataType.GetSizeInMemoryBytes() == 2 )
					ctx->bc.InstrSHORT_W(BC_SetV2, (short)offset, ctx->type.wordValue);
				else if( ctx->type.dataType.GetSizeInMemoryBytes() == 4 )
					ctx->bc.InstrSHORT_DW(BC_SetV4, (short)offset, ctx->type.dwordValue);
				else
					ctx->bc.InstrSHORT_QW(BC_SetV8, (short)offset, ctx->type.qwordValue);

				ctx->type.SetVariable(ctx->type.dataType, offset, true);
				return;
			}
			else
			{
				assert(ctx->type.dataType.IsPrimitive());
				assert(ctx->type.dataType.IsReference());

				ctx->type.dataType.MakeReference(false);
				offset = AllocateVariableNotIn(ctx->type.dataType, true, &excludeVars);

				// Read the value from the address in the register directly into the variable
				if( ctx->type.dataType.GetSizeInMemoryBytes() == 1 )
					ctx->bc.InstrSHORT(BC_RDR1, (short)offset);
				else if( ctx->type.dataType.GetSizeInMemoryBytes() == 2 )
					ctx->bc.InstrSHORT(BC_RDR2, (short)offset);
				else if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
					ctx->bc.InstrSHORT(BC_RDR4, (short)offset);
				else
					ctx->bc.InstrSHORT(BC_RDR8, (short)offset);
			}

			ReleaseTemporaryVariable(ctx->type, &ctx->bc);
			ctx->type.SetVariable(ctx->type.dataType, offset, true);
		}
	}
}

void asCCompiler::ConvertToVariableNotIn(asSExprContext *ctx, asSExprContext *exclude)
{
	asCArray<int> excludeVars;
	if( exclude ) exclude->bc.GetVarsUsed(excludeVars);
	ConvertToVariableNotIn(ctx, &excludeVars);
}


void asCCompiler::CompileMathOperator(asCScriptNode *node, asSExprContext *lctx, asSExprContext  *rctx, asSExprContext *ctx)
{
	// TODO: If a constant is only using 32bits, then a 32bit operation is preferred

	// Implicitly convert the operands to a number type
	asCDataType to;
	if( lctx->type.dataType.IsDoubleType() || rctx->type.dataType.IsDoubleType() )
		to.SetTokenType(ttDouble);
	else if( lctx->type.dataType.IsFloatType() || rctx->type.dataType.IsFloatType() )
		to.SetTokenType(ttFloat);
	else if( lctx->type.dataType.GetSizeInMemoryDWords() == 2 || rctx->type.dataType.GetSizeInMemoryDWords() == 2 )
	{
		if( lctx->type.dataType.IsIntegerType() || rctx->type.dataType.IsIntegerType() )
			to.SetTokenType(ttInt64);
		else if( lctx->type.dataType.IsUnsignedType() || rctx->type.dataType.IsUnsignedType() )
			to.SetTokenType(ttUInt64);
	}
	else
	{
		if( lctx->type.dataType.IsIntegerType() || rctx->type.dataType.IsIntegerType() )
			to.SetTokenType(ttInt);
		else if( lctx->type.dataType.IsUnsignedType() || rctx->type.dataType.IsUnsignedType() )
			to.SetTokenType(ttUInt);
	}

	// If doing an operation with double constant and float variable, the constant should be converted to float
	if( (lctx->type.isConstant && lctx->type.dataType.IsDoubleType() && !rctx->type.isConstant && rctx->type.dataType.IsFloatType()) ||
		(rctx->type.isConstant && rctx->type.dataType.IsDoubleType() && !lctx->type.isConstant && lctx->type.dataType.IsFloatType()) )
		to.SetTokenType(ttFloat);

	// Do the actual conversion
	asCArray<int> reservedVars;
	rctx->bc.GetVarsUsed(reservedVars);
	ImplicitConversion(lctx, to, node, false, true, &reservedVars);
	ImplicitConversion(rctx, to, node, false);

	// Verify that the conversion was successful
	if( !lctx->type.dataType.IsIntegerType() &&
		!lctx->type.dataType.IsUnsignedType() &&
		!lctx->type.dataType.IsFloatType() &&
		!lctx->type.dataType.IsDoubleType() )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_MATH_TYPE, lctx->type.dataType.Format().AddressOf());
		Error(str.AddressOf(), node);

		ctx->type.SetDummy();
		return;
	}

	if( !rctx->type.dataType.IsIntegerType() &&
		!rctx->type.dataType.IsUnsignedType() &&
		!rctx->type.dataType.IsFloatType() &&
		!rctx->type.dataType.IsDoubleType() )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_MATH_TYPE, rctx->type.dataType.Format().AddressOf());
		Error(str.AddressOf(), node);

		ctx->type.SetDummy();
		return;
	}

	bool isConstant = lctx->type.isConstant && rctx->type.isConstant;

	// Verify if we are dividing with a constant zero
	int op = node->tokenType;
	if( rctx->type.isConstant && rctx->type.qwordValue == 0 &&
		(op == ttSlash   || op == ttDivAssign ||
		 op == ttPercent || op == ttModAssign) )
	{
		Error(TXT_DIVIDE_BY_ZERO, node);
	}

	if( !isConstant )
	{
		ConvertToVariableNotIn(lctx, rctx);
		ConvertToVariable(rctx);
		ReleaseTemporaryVariable(lctx->type, &lctx->bc);
		ReleaseTemporaryVariable(rctx->type, &rctx->bc);

		if( op == ttAddAssign || op == ttSubAssign ||
			op == ttMulAssign || op == ttDivAssign ||
			op == ttModAssign )
		{
			// Merge the operands in the different order so that they are evaluated correctly
			MergeExprContexts(ctx, rctx);
			MergeExprContexts(ctx, lctx);
		}
		else
		{
			MergeExprContexts(ctx, lctx);
			MergeExprContexts(ctx, rctx);
		}

		bcInstr instruction = BC_ADDi;
		if( lctx->type.dataType.IsIntegerType() ||
			lctx->type.dataType.IsUnsignedType() )
		{
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				if( op == ttPlus || op == ttAddAssign )
					instruction = BC_ADDi;
				else if( op == ttMinus || op == ttSubAssign )
					instruction = BC_SUBi;
				else if( op == ttStar || op == ttMulAssign )
					instruction = BC_MULi;
				else if( op == ttSlash || op == ttDivAssign )
					instruction = BC_DIVi;
				else if( op == ttPercent || op == ttModAssign )
					instruction = BC_MODi;
			}
			else
			{
				if( op == ttPlus || op == ttAddAssign )
					instruction = BC_ADDi64;
				else if( op == ttMinus || op == ttSubAssign )
					instruction = BC_SUBi64;
				else if( op == ttStar || op == ttMulAssign )
					instruction = BC_MULi64;
				else if( op == ttSlash || op == ttDivAssign )
					instruction = BC_DIVi64;
				else if( op == ttPercent || op == ttModAssign )
					instruction = BC_MODi64;
			}
		}
		else if( lctx->type.dataType.IsFloatType() )
		{
			if( op == ttPlus || op == ttAddAssign )
				instruction = BC_ADDf;
			else if( op == ttMinus || op == ttSubAssign )
				instruction = BC_SUBf;
			else if( op == ttStar || op == ttMulAssign )
				instruction = BC_MULf;
			else if( op == ttSlash || op == ttDivAssign )
				instruction = BC_DIVf;
			else if( op == ttPercent || op == ttModAssign )
				instruction = BC_MODf;
		}
		else if( lctx->type.dataType.IsDoubleType() )
		{
			if( op == ttPlus || op == ttAddAssign )
				instruction = BC_ADDd;
			else if( op == ttMinus || op == ttSubAssign )
				instruction = BC_SUBd;
			else if( op == ttStar || op == ttMulAssign )
				instruction = BC_MULd;
			else if( op == ttSlash || op == ttDivAssign )
				instruction = BC_DIVd;
			else if( op == ttPercent || op == ttModAssign )
				instruction = BC_MODd;
		}
		else
		{
			// Shouldn't be possible
			assert(false);
		}

		// Do the operation
		int a = AllocateVariable(lctx->type.dataType, true);
		int b = lctx->type.stackOffset;
		int c = rctx->type.stackOffset;

		ctx->bc.InstrW_W_W(instruction, a, b, c);

		ctx->type.SetVariable(lctx->type.dataType, a, true);
	}
	else
	{
		// Both values are constants
		if( lctx->type.dataType.IsIntegerType() ||
			lctx->type.dataType.IsUnsignedType() )
		{
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				int v = 0;
				if( op == ttPlus )
					v = lctx->type.intValue + rctx->type.intValue;
				else if( op == ttMinus )
					v = lctx->type.intValue - rctx->type.intValue;
				else if( op == ttStar )
					v = lctx->type.intValue * rctx->type.intValue;
				else if( op == ttSlash )
				{
					if( rctx->type.intValue == 0 )
						v = 0;
					else
						v = lctx->type.intValue / rctx->type.intValue;
				}
				else if( op == ttPercent )
				{
					if( rctx->type.intValue == 0 )
						v = 0;
					else
						v = lctx->type.intValue % rctx->type.intValue;
				}

				ctx->type.SetConstantDW(lctx->type.dataType, v);
			}
			else
			{
				asQWORD v = 0;
				if( op == ttPlus )
					v = lctx->type.qwordValue + rctx->type.qwordValue;
				else if( op == ttMinus )
					v = lctx->type.qwordValue - rctx->type.qwordValue;
				else if( op == ttStar )
					v = lctx->type.qwordValue * rctx->type.qwordValue;
				else if( op == ttSlash )
				{
					if( rctx->type.qwordValue == 0 )
						v = 0;
					else
						v = lctx->type.qwordValue / rctx->type.qwordValue;
				}
				else if( op == ttPercent )
				{
					if( rctx->type.qwordValue == 0 )
						v = 0;
					else
						v = lctx->type.qwordValue % rctx->type.qwordValue;
				}

				ctx->type.SetConstantQW(lctx->type.dataType, v);
			}
		}
		else if( lctx->type.dataType.IsFloatType() )
		{
			float v = 0.0f;
			if( op == ttPlus )
				v = lctx->type.floatValue + rctx->type.floatValue;
			else if( op == ttMinus )
				v = lctx->type.floatValue - rctx->type.floatValue;
			else if( op == ttStar )
				v = lctx->type.floatValue * rctx->type.floatValue;
			else if( op == ttSlash )
			{
				if( rctx->type.floatValue == 0 )
					v = 0;
				else
					v = lctx->type.floatValue / rctx->type.floatValue;
			}
			else if( op == ttPercent )
			{
				if( rctx->type.floatValue == 0 )
					v = 0;
				else
					v = fmodf(lctx->type.floatValue, rctx->type.floatValue);
			}

			ctx->type.SetConstantF(lctx->type.dataType, v);
		}
		else if( lctx->type.dataType.IsDoubleType() )
		{
			double v = 0.0;
			if( op == ttPlus )
				v = lctx->type.doubleValue + rctx->type.doubleValue;
			else if( op == ttMinus )
				v = lctx->type.doubleValue - rctx->type.doubleValue;
			else if( op == ttStar )
				v = lctx->type.doubleValue * rctx->type.doubleValue;
			else if( op == ttSlash )
			{
				if( rctx->type.doubleValue == 0 )
					v = 0;
				else
					v = lctx->type.doubleValue / rctx->type.doubleValue;
			}
			else if( op == ttPercent )
			{
				if( rctx->type.doubleValue == 0 )
					v = 0;
				else
					v = fmod(lctx->type.doubleValue, rctx->type.doubleValue);
			}

			ctx->type.SetConstantD(lctx->type.dataType, v);
		}
		else
		{
			// Shouldn't be possible
			assert(false);
		}
	}
}

void asCCompiler::CompileBitwiseOperator(asCScriptNode *node, asSExprContext *lctx, asSExprContext *rctx, asSExprContext *ctx)
{
	// TODO: If a constant is only using 32bits, then a 32bit operation is preferred

	int op = node->tokenType;
	if( op == ttAmp    || op == ttAndAssign ||
		op == ttBitOr  || op == ttOrAssign  ||
		op == ttBitXor || op == ttXorAssign )
	{
		// Convert left hand operand to integer if it's not already one
		if( !(lctx->type.dataType.IsIntegerType() ||
			  lctx->type.dataType.IsUnsignedType()) )
		{
			asCDataType to;
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 2  )
				to.SetTokenType(ttInt64);
			else
				to.SetTokenType(ttInt);

			// Do the actual conversion
			asCArray<int> reservedVars;
			rctx->bc.GetVarsUsed(reservedVars);
			ImplicitConversion(lctx, to, node, false, true, &reservedVars);

			// Verify that the conversion was successful
			if( !lctx->type.dataType.IsUnsignedType() )
			{
				asCString str;
				str.Format(TXT_NO_CONVERSION_s_TO_s, lctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
				Error(str.AddressOf(), node);
			}
		}

		// Convert right hand operand to same type as left hand operand
		ImplicitConversion(rctx, lctx->type.dataType, node, false);
		if( !rctx->type.dataType.IsEqualExceptRef(lctx->type.dataType) )
		{
			asCString str;
			str.Format(TXT_NO_CONVERSION_s_TO_s, rctx->type.dataType.Format().AddressOf(), lctx->type.dataType.Format().AddressOf());
			Error(str.AddressOf(), node);
		}

		bool isConstant = lctx->type.isConstant && rctx->type.isConstant;

		if( !isConstant )
		{
			ConvertToVariableNotIn(lctx, rctx);
			ConvertToVariable(rctx);
			ReleaseTemporaryVariable(lctx->type, &lctx->bc);
			ReleaseTemporaryVariable(rctx->type, &rctx->bc);

			if( op == ttAndAssign || op == ttOrAssign || op == ttXorAssign )
			{
				// Compound assignments execute the right hand value first
				MergeExprContexts(ctx, rctx);
				MergeExprContexts(ctx, lctx);
			}
			else
			{
				MergeExprContexts(ctx, lctx);
				MergeExprContexts(ctx, rctx);
			}

			bcInstr instruction = BC_BAND;
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				if( op == ttAmp || op == ttAndAssign )
					instruction = BC_BAND;
				else if( op == ttBitOr || op == ttOrAssign )
					instruction = BC_BOR;
				else if( op == ttBitXor || op == ttXorAssign )
					instruction = BC_BXOR;
			}
			else
			{
				if( op == ttAmp || op == ttAndAssign )
					instruction = BC_BAND64;
				else if( op == ttBitOr || op == ttOrAssign )
					instruction = BC_BOR64;
				else if( op == ttBitXor || op == ttXorAssign )
					instruction = BC_BXOR64;
			}

			// Do the operation
			int a = AllocateVariable(lctx->type.dataType, true);
			int b = lctx->type.stackOffset;
			int c = rctx->type.stackOffset;

			ctx->bc.InstrW_W_W(instruction, a, b, c);

			ctx->type.SetVariable(lctx->type.dataType, a, true);
		}
		else
		{
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 2 )
			{
				asQWORD v = 0;
				if( op == ttAmp )
					v = lctx->type.qwordValue & rctx->type.qwordValue;
				else if( op == ttBitOr )
					v = lctx->type.qwordValue | rctx->type.qwordValue;
				else if( op == ttBitXor )
					v = lctx->type.qwordValue ^ rctx->type.qwordValue;

				// Remember the result
				ctx->type.SetConstantQW(lctx->type.dataType, v);
			}
			else
			{
				asDWORD v = 0;
				if( op == ttAmp )
					v = lctx->type.dwordValue & rctx->type.dwordValue;
				else if( op == ttBitOr )
					v = lctx->type.dwordValue | rctx->type.dwordValue;
				else if( op == ttBitXor )
					v = lctx->type.dwordValue ^ rctx->type.dwordValue;

				// Remember the result
				ctx->type.SetConstantDW(lctx->type.dataType, v);
			}
		}
	}
	else if( op == ttBitShiftLeft       || op == ttShiftLeftAssign   ||
		     op == ttBitShiftRight      || op == ttShiftRightLAssign ||
			 op == ttBitShiftRightArith || op == ttShiftRightAAssign )
	{
		// Convert left hand operand to integer if it's not already one
		if( !(lctx->type.dataType.IsIntegerType() ||
			  lctx->type.dataType.IsUnsignedType()) )
		{
			asCDataType to;
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 2  )
				to.SetTokenType(ttInt64);
			else
				to.SetTokenType(ttInt);

			// Do the actual conversion
			asCArray<int> reservedVars;
			rctx->bc.GetVarsUsed(reservedVars);
			ImplicitConversion(lctx, to, node, false, true, &reservedVars);

			// Verify that the conversion was successful
			if( !lctx->type.dataType.IsUnsignedType() )
			{
				asCString str;
				str.Format(TXT_NO_CONVERSION_s_TO_s, lctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
				Error(str.AddressOf(), node);
			}
		}

		// Right operand must be 32bit uint
		ImplicitConversion(rctx, asCDataType::CreatePrimitive(ttUInt, true), node, false);
		if( !rctx->type.dataType.IsUnsignedType() )
		{
			asCString str;
			str.Format(TXT_NO_CONVERSION_s_TO_s, rctx->type.dataType.Format().AddressOf(), "uint");
			Error(str.AddressOf(), node);
		}

		bool isConstant = lctx->type.isConstant && rctx->type.isConstant;

		if( !isConstant )
		{
			ConvertToVariableNotIn(lctx, rctx);
			ConvertToVariable(rctx);
			ReleaseTemporaryVariable(lctx->type, &lctx->bc);
			ReleaseTemporaryVariable(rctx->type, &rctx->bc);

			if( op == ttShiftLeftAssign || op == ttShiftRightLAssign || op == ttShiftRightAAssign )
			{
				// Compound assignments execute the right hand value first
				MergeExprContexts(ctx, rctx);
				MergeExprContexts(ctx, lctx);
			}
			else
			{
				MergeExprContexts(ctx, lctx);
				MergeExprContexts(ctx, rctx);
			}

			bcInstr instruction = BC_BSLL;
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				if( op == ttBitShiftLeft || op == ttShiftLeftAssign )
					instruction = BC_BSLL;
				else if( op == ttBitShiftRight || op == ttShiftRightLAssign )
					instruction = BC_BSRL;
				else if( op == ttBitShiftRightArith || op == ttShiftRightAAssign )
					instruction = BC_BSRA;
			}
			else
			{
				if( op == ttBitShiftLeft || op == ttShiftLeftAssign )
					instruction = BC_BSLL64;
				else if( op == ttBitShiftRight || op == ttShiftRightLAssign )
					instruction = BC_BSRL64;
				else if( op == ttBitShiftRightArith || op == ttShiftRightAAssign )
					instruction = BC_BSRA64;
			}

			// Do the operation
			int a = AllocateVariable(lctx->type.dataType, true);
			int b = lctx->type.stackOffset;
			int c = rctx->type.stackOffset;

			ctx->bc.InstrW_W_W(instruction, a, b, c);

			ctx->type.SetVariable(lctx->type.dataType, a, true);
		}
		else
		{
			if( lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				asDWORD v = 0;
				if( op == ttBitShiftLeft )
					v = lctx->type.dwordValue << rctx->type.dwordValue;
				else if( op == ttBitShiftRight )
					v = lctx->type.dwordValue >> rctx->type.dwordValue;
				else if( op == ttBitShiftRightArith )
					v = lctx->type.intValue >> rctx->type.dwordValue;

				ctx->type.SetConstantDW(lctx->type.dataType, v);
			}
			else
			{
				asQWORD v = 0;
				if( op == ttBitShiftLeft )
					v = lctx->type.qwordValue << rctx->type.dwordValue;
				else if( op == ttBitShiftRight )
					v = lctx->type.qwordValue >> rctx->type.dwordValue;
				else if( op == ttBitShiftRightArith )
					v = asINT64(lctx->type.qwordValue) >> rctx->type.dwordValue;

				ctx->type.SetConstantQW(lctx->type.dataType, v);
			}
		}
	}
}

void asCCompiler::CompileComparisonOperator(asCScriptNode *node, asSExprContext *lctx, asSExprContext *rctx, asSExprContext *ctx)
{
	// Both operands must be of the same type

	// Implicitly convert the operands to a number type
	asCDataType to;
	if( lctx->type.dataType.IsDoubleType() || rctx->type.dataType.IsDoubleType() )
		to.SetTokenType(ttDouble);
	else if( lctx->type.dataType.IsFloatType() || rctx->type.dataType.IsFloatType() )
		to.SetTokenType(ttFloat);
	else if( lctx->type.dataType.GetSizeInMemoryDWords() == 2 || rctx->type.dataType.GetSizeInMemoryDWords() == 2 )
	{
		if( lctx->type.dataType.IsIntegerType() || rctx->type.dataType.IsIntegerType() )
			to.SetTokenType(ttInt64);
		else if( lctx->type.dataType.IsUnsignedType() || rctx->type.dataType.IsUnsignedType() )
			to.SetTokenType(ttUInt64);
	}
	else
	{
		if( lctx->type.dataType.IsIntegerType() || rctx->type.dataType.IsIntegerType() )
			to.SetTokenType(ttInt);
		else if( lctx->type.dataType.IsUnsignedType() || rctx->type.dataType.IsUnsignedType() )
			to.SetTokenType(ttUInt);
		else if( lctx->type.dataType.IsBooleanType() || rctx->type.dataType.IsBooleanType() )
			to.SetTokenType(ttBool);
	}

	// If doing an operation with double constant and float variable, the constant should be converted to float
	if( (lctx->type.isConstant && lctx->type.dataType.IsDoubleType() && !rctx->type.isConstant && rctx->type.dataType.IsFloatType()) ||
		(rctx->type.isConstant && rctx->type.dataType.IsDoubleType() && !lctx->type.isConstant && lctx->type.dataType.IsFloatType()) )
		to.SetTokenType(ttFloat);

	// Is it an operation on signed values?
	bool signMismatch = false;
	if( !lctx->type.dataType.IsUnsignedType() || !rctx->type.dataType.IsUnsignedType() )
	{
		if( lctx->type.dataType.GetTokenType() == ttUInt64 )
		{
			if( !lctx->type.isConstant )
				signMismatch = true;
			else if( lctx->type.qwordValue & (I64(1)<<63) )
				signMismatch = true;
		}
		if( lctx->type.dataType.GetTokenType() == ttUInt )
		{
			if( !lctx->type.isConstant )
				signMismatch = true;
			else if( lctx->type.dwordValue & (1<<31) )
				signMismatch = true;
		}
		if( rctx->type.dataType.GetTokenType() == ttUInt64 )
		{
			if( !rctx->type.isConstant )
				signMismatch = true;
			else if( rctx->type.qwordValue & (I64(1)<<63) )
				signMismatch = true;
		}
		if( rctx->type.dataType.GetTokenType() == ttUInt )
		{
			if( !rctx->type.isConstant )
				signMismatch = true;
			else if( rctx->type.dwordValue & (1<<31) )
				signMismatch = true;
		}
	}

	// Check for signed/unsigned mismatch
	if( signMismatch )
		Warning(TXT_SIGNED_UNSIGNED_MISMATCH, node);

	// Do the actual conversion
	asCArray<int> reservedVars;
	rctx->bc.GetVarsUsed(reservedVars);
	ImplicitConversion(lctx, to, node, false, true, &reservedVars);
	ImplicitConversion(rctx, to, node, false);

	// Verify that the conversion was successful
	if( !lctx->type.dataType.IsEqualExceptConst(to) )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_s, lctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
		Error(str.AddressOf(), node);
	}

	if( !rctx->type.dataType.IsEqualExceptConst(to) )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_s, rctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
		Error(str.AddressOf(), node);
	}

	bool isConstant = lctx->type.isConstant && rctx->type.isConstant;
	int op = node->tokenType;

	if( !isConstant )
	{
		if( to.IsBooleanType() )
		{
			int op = node->tokenType;
			if( op == ttEqual || op == ttNotEqual )
			{
				// Must convert to temporary variable, because we are changing the value before comparison
				ConvertToTempVariableNotIn(lctx, rctx);
				ConvertToTempVariable(rctx);
				ReleaseTemporaryVariable(lctx->type, &lctx->bc);
				ReleaseTemporaryVariable(rctx->type, &rctx->bc);

				// Make sure they are equal if not false
				lctx->bc.InstrWORD(BC_NOT, lctx->type.stackOffset);
				rctx->bc.InstrWORD(BC_NOT, rctx->type.stackOffset);

				MergeExprContexts(ctx, lctx);
				MergeExprContexts(ctx, rctx);

				int a = AllocateVariable(ctx->type.dataType, true);
				int b = lctx->type.stackOffset;
				int c = rctx->type.stackOffset;

				if( op == ttEqual )
				{
					ctx->bc.InstrW_W(BC_CMPi,b,c);
					ctx->bc.Instr(BC_TZ);
					ctx->bc.InstrSHORT(BC_CpyRtoV4, (short)a);
				}
				else if( op == ttNotEqual )
				{
					ctx->bc.InstrW_W(BC_CMPi,b,c);
					ctx->bc.Instr(BC_TNZ);
					ctx->bc.InstrSHORT(BC_CpyRtoV4, (short)a);
				}

				ctx->type.SetVariable(asCDataType::CreatePrimitive(ttBool, true), a, true);
			}
			else
			{
				// TODO: Use TXT_ILLEGAL_OPERATION_ON
				Error(TXT_ILLEGAL_OPERATION, node);
			}
		}
		else
		{
			ConvertToVariableNotIn(lctx, rctx);
			ConvertToVariable(rctx);
			ReleaseTemporaryVariable(lctx->type, &lctx->bc);
			ReleaseTemporaryVariable(rctx->type, &rctx->bc);

			MergeExprContexts(ctx, lctx);
			MergeExprContexts(ctx, rctx);

			bcInstr iCmp = BC_CMPi, iT = BC_TZ;

			if( lctx->type.dataType.IsIntegerType() && lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				iCmp = BC_CMPi;
			else if( lctx->type.dataType.IsUnsignedType() && lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
				iCmp = BC_CMPu;
			else if( lctx->type.dataType.IsIntegerType() && lctx->type.dataType.GetSizeInMemoryDWords() == 2 )
				iCmp = BC_CMPi64;
			else if( lctx->type.dataType.IsUnsignedType() && lctx->type.dataType.GetSizeInMemoryDWords() == 2 )
				iCmp = BC_CMPu64;
			else if( lctx->type.dataType.IsFloatType() )
				iCmp = BC_CMPf;
			else if( lctx->type.dataType.IsDoubleType() )
				iCmp = BC_CMPd;
			else
				assert(false);

			if( op == ttEqual )
				iT = BC_TZ;
			else if( op == ttNotEqual )
				iT = BC_TNZ;
			else if( op == ttLessThan )
				iT = BC_TS;
			else if( op == ttLessThanOrEqual )
				iT = BC_TNP;
			else if( op == ttGreaterThan )
				iT = BC_TP;
			else if( op == ttGreaterThanOrEqual )
				iT = BC_TNS;

			int a = AllocateVariable(ctx->type.dataType, true);
			int b = lctx->type.stackOffset;
			int c = rctx->type.stackOffset;

			ctx->bc.InstrW_W(iCmp, b, c);
			ctx->bc.Instr(iT);
			ctx->bc.InstrSHORT(BC_CpyRtoV4, (short)a);

			ctx->type.SetVariable(asCDataType::CreatePrimitive(ttBool, true), a, true);
		}
	}
	else
	{
		if( to.IsBooleanType() )
		{
			int op = node->tokenType;
			if( op == ttEqual || op == ttNotEqual )
			{
				// Make sure they are equal if not false
				if( lctx->type.dwordValue != 0 ) lctx->type.dwordValue = VALUE_OF_BOOLEAN_TRUE;
				if( rctx->type.dwordValue != 0 ) rctx->type.dwordValue = VALUE_OF_BOOLEAN_TRUE;

				asDWORD v = 0;
				if( op == ttEqual )
				{
					v = lctx->type.intValue - rctx->type.intValue;
					if( v == 0 ) v = VALUE_OF_BOOLEAN_TRUE; else v = 0;
				}
				else if( op == ttNotEqual )
				{
					v = lctx->type.intValue - rctx->type.intValue;
					if( v != 0 ) v = VALUE_OF_BOOLEAN_TRUE; else v = 0;
				}

				ctx->type.SetConstantDW(asCDataType::CreatePrimitive(ttBool, true), v);
			}
			else
			{
				// TODO: Use TXT_ILLEGAL_OPERATION_ON
				Error(TXT_ILLEGAL_OPERATION, node);
			}
		}
		else
		{
			int i = 0;
			if( lctx->type.dataType.IsIntegerType() && lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				int v = lctx->type.intValue - rctx->type.intValue;
				if( v < 0 ) i = -1;
				if( v > 0 ) i = 1;
			}
			else if( lctx->type.dataType.IsUnsignedType() && lctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			{
				asDWORD v1 = lctx->type.dwordValue;
				asDWORD v2 = rctx->type.dwordValue;
				if( v1 < v2 ) i = -1;
				if( v1 > v2 ) i = 1;
			}
			else if( lctx->type.dataType.IsIntegerType() && lctx->type.dataType.GetSizeInMemoryDWords() == 2 )
			{
				asINT64 v = asINT64(lctx->type.qwordValue) - asINT64(rctx->type.qwordValue);
				if( v < 0 ) i = -1;
				if( v > 0 ) i = 1;
			}
			else if( lctx->type.dataType.IsUnsignedType() && lctx->type.dataType.GetSizeInMemoryDWords() == 2 )
			{
				asQWORD v1 = lctx->type.qwordValue;
				asQWORD v2 = rctx->type.qwordValue;
				if( v1 < v2 ) i = -1;
				if( v1 > v2 ) i = 1;
			}
			else if( lctx->type.dataType.IsFloatType() )
			{
				float v = lctx->type.floatValue - rctx->type.floatValue;
				if( v < 0 ) i = -1;
				if( v > 0 ) i = 1;
			}
			else if( lctx->type.dataType.IsDoubleType() )
			{
				double v = lctx->type.doubleValue - rctx->type.doubleValue;
				if( v < 0 ) i = -1;
				if( v > 0 ) i = 1;
			}


			if( op == ttEqual )
				i = (i == 0 ? VALUE_OF_BOOLEAN_TRUE : 0);
			else if( op == ttNotEqual )
				i = (i != 0 ? VALUE_OF_BOOLEAN_TRUE : 0);
			else if( op == ttLessThan )
				i = (i < 0 ? VALUE_OF_BOOLEAN_TRUE : 0);
			else if( op == ttLessThanOrEqual )
				i = (i <= 0 ? VALUE_OF_BOOLEAN_TRUE : 0);
			else if( op == ttGreaterThan )
				i = (i > 0 ? VALUE_OF_BOOLEAN_TRUE : 0);
			else if( op == ttGreaterThanOrEqual )
				i = (i >= 0 ? VALUE_OF_BOOLEAN_TRUE : 0);

			ctx->type.SetConstantDW(asCDataType::CreatePrimitive(ttBool, true), i);
		}
	}
}

void asCCompiler::PushVariableOnStack(asSExprContext *ctx, bool asReference)
{
	// Put the result on the stack
	ctx->bc.InstrSHORT(BC_PSF, ctx->type.stackOffset);
	if( asReference )
		ctx->type.dataType.MakeReference(true);
	else
	{
		if( ctx->type.dataType.GetSizeInMemoryDWords() == 1 )
			ctx->bc.Instr(BC_RDS4);
		else
			ctx->bc.Instr(BC_RDS8);
	}
}

void asCCompiler::CompileBooleanOperator(asCScriptNode *node, asSExprContext *lctx, asSExprContext *rctx, asSExprContext *ctx)
{
	// Both operands must be booleans
	asCDataType to;
	to.SetTokenType(ttBool);

	// Do the actual conversion
	asCArray<int> reservedVars;
	rctx->bc.GetVarsUsed(reservedVars);
	ImplicitConversion(lctx, to, node, false, true, &reservedVars);
	ImplicitConversion(rctx, to, node, false);

	// Verify that the conversion was successful
	if( !lctx->type.dataType.IsBooleanType() )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_s, lctx->type.dataType.Format().AddressOf(), "bool");
		Error(str.AddressOf(), node);
	}

	if( !rctx->type.dataType.IsBooleanType() )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_s, rctx->type.dataType.Format().AddressOf(), "bool");
		Error(str.AddressOf(), node);
	}

	bool isConstant = lctx->type.isConstant && rctx->type.isConstant;

	ctx->type.Set(asCDataType::CreatePrimitive(ttBool, true));

	// What kind of operator is it?
	int op = node->tokenType;
	if( op == ttXor )
	{
		if( !isConstant )
		{
			// Must convert to temporary variable, because we are changing the value before comparison
			ConvertToTempVariableNotIn(lctx, rctx);
			ConvertToTempVariable(rctx);
			ReleaseTemporaryVariable(lctx->type, &lctx->bc);
			ReleaseTemporaryVariable(rctx->type, &rctx->bc);

			// Make sure they are equal if not false
			lctx->bc.InstrWORD(BC_NOT, lctx->type.stackOffset);
			rctx->bc.InstrWORD(BC_NOT, rctx->type.stackOffset);

			MergeExprContexts(ctx, lctx);
			MergeExprContexts(ctx, rctx);

			int a = AllocateVariable(ctx->type.dataType, true);
			int b = lctx->type.stackOffset;
			int c = rctx->type.stackOffset;

			ctx->bc.InstrW_W_W(BC_BXOR,a,b,c);

			ctx->type.SetVariable(asCDataType::CreatePrimitive(ttBool, true), a, true);
		}
		else
		{
			// Make sure they are equal if not false
#if AS_SIZEOF_BOOL == 1
			if( lctx->type.byteValue != 0 ) lctx->type.byteValue = VALUE_OF_BOOLEAN_TRUE;
			if( rctx->type.byteValue != 0 ) rctx->type.byteValue = VALUE_OF_BOOLEAN_TRUE;

			asBYTE v = 0;
			v = lctx->type.byteValue - rctx->type.byteValue;
			if( v != 0 ) v = VALUE_OF_BOOLEAN_TRUE; else v = 0;

			ctx->type.isConstant = true;
			ctx->type.byteValue = v;
#else
			if( lctx->type.dwordValue != 0 ) lctx->type.dwordValue = VALUE_OF_BOOLEAN_TRUE;
			if( rctx->type.dwordValue != 0 ) rctx->type.dwordValue = VALUE_OF_BOOLEAN_TRUE;

			asDWORD v = 0;
			v = lctx->type.intValue - rctx->type.intValue;
			if( v != 0 ) v = VALUE_OF_BOOLEAN_TRUE; else v = 0;

			ctx->type.isConstant = true;
			ctx->type.dwordValue = v;
#endif
		}
	}
	else if( op == ttAnd ||
			 op == ttOr )
	{
		if( !isConstant )
		{
			// If or-operator and first value is 1 the second value shouldn't be calculated
			// if and-operator and first value is 0 the second value shouldn't be calculated
			ConvertToVariable(lctx);
			ReleaseTemporaryVariable(lctx->type, &lctx->bc);
			MergeExprContexts(ctx, lctx);

			int offset = AllocateVariable(asCDataType::CreatePrimitive(ttBool, false), true);

			int label1 = nextLabel++;
			int label2 = nextLabel++;
			if( op == ttAnd )
			{
				ctx->bc.InstrSHORT(BC_CpyVtoR4, lctx->type.stackOffset);
				ctx->bc.Instr(BC_ClrHi);
				ctx->bc.InstrDWORD(BC_JNZ, label1);
				ctx->bc.InstrW_DW(BC_SetV4, (asWORD)offset, 0);
				ctx->bc.InstrINT(BC_JMP, label2);
			}
			else if( op == ttOr )
			{
				ctx->bc.InstrSHORT(BC_CpyVtoR4, lctx->type.stackOffset);
				ctx->bc.Instr(BC_ClrHi);
				ctx->bc.InstrDWORD(BC_JZ, label1);
#if AS_SIZEOF_BOOL == 1
				ctx->bc.InstrSHORT_B(BC_SetV1, (short)offset, VALUE_OF_BOOLEAN_TRUE);
#else
				ctx->bc.InstrSHORT_DW(BC_SetV4, (short)offset, VALUE_OF_BOOLEAN_TRUE);
#endif
				ctx->bc.InstrINT(BC_JMP, label2);
			}

			ctx->bc.Label((short)label1);
			ConvertToVariable(rctx);
			ReleaseTemporaryVariable(rctx->type, &rctx->bc);
			rctx->bc.InstrW_W(BC_CpyVtoV4, offset, rctx->type.stackOffset);
			MergeExprContexts(ctx, rctx);
			ctx->bc.Label((short)label2);

			ctx->type.SetVariable(asCDataType::CreatePrimitive(ttBool, false), offset, true);
		}
		else
		{
#if AS_SIZEOF_BOOL == 1
			asBYTE v = 0;
			if( op == ttAnd )
				v = lctx->type.byteValue && rctx->type.byteValue;
			else if( op == ttOr )
				v = lctx->type.byteValue || rctx->type.byteValue;

			// Remember the result
			ctx->type.isConstant = true;
			ctx->type.byteValue = v;
#else
			asDWORD v = 0;
			if( op == ttAnd )
				v = lctx->type.dwordValue && rctx->type.dwordValue;
			else if( op == ttOr )
				v = lctx->type.dwordValue || rctx->type.dwordValue;

			// Remember the result
			ctx->type.isConstant = true;
			ctx->type.dwordValue = v;
#endif
		}
	}
}

void asCCompiler::CompileOperatorOnHandles(asCScriptNode *node, asSExprContext *lctx, asSExprContext *rctx, asSExprContext *ctx)
{
	// Implicitly convert null to the other type
	asCDataType to;
	if( lctx->type.IsNullConstant() )
		to = rctx->type.dataType;
	else if( rctx->type.IsNullConstant() )
		to = lctx->type.dataType;
	else
	{
		// TODO: Use the common base type
		to = lctx->type.dataType;
	}

	// Need to pop the value if it is a null constant
	if( lctx->type.IsNullConstant() )
		lctx->bc.Pop(PTR_SIZE);
	if( rctx->type.IsNullConstant() )
		rctx->bc.Pop(PTR_SIZE);

	// Convert both sides to explicit handles
	to.MakeHandle(true);
	to.MakeReference(false);

	// Do the conversion
	ImplicitConversion(lctx, to, node, false);
	ImplicitConversion(rctx, to, node, false);

	// Both operands must be of the same type

	// Verify that the conversion was successful
	if( !lctx->type.dataType.IsEqualExceptConst(to) )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_s, lctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
		Error(str.AddressOf(), node);
	}

	if( !rctx->type.dataType.IsEqualExceptConst(to) )
	{
		asCString str;
		str.Format(TXT_NO_CONVERSION_s_TO_s, rctx->type.dataType.Format().AddressOf(), to.Format().AddressOf());
		Error(str.AddressOf(), node);
	}

	ctx->type.Set(asCDataType::CreatePrimitive(ttBool, true));

	int op = node->tokenType;
	if( op == ttEqual || op == ttNotEqual )
	{
		// If the object handle already is in a variable we must manually pop it from the stack
		if( lctx->type.isVariable )
			lctx->bc.Pop(PTR_SIZE);
		if( rctx->type.isVariable )
			rctx->bc.Pop(PTR_SIZE);

		// TODO: Treat the object handles as two integers, i.e. don't do REFCPY
		ConvertToVariableNotIn(lctx, rctx);
		ConvertToVariable(rctx);

		MergeExprContexts(ctx, lctx);
		MergeExprContexts(ctx, rctx);

		int a = AllocateVariable(ctx->type.dataType, true);
		int b = lctx->type.stackOffset;
		int c = rctx->type.stackOffset;

		if( op == ttEqual )
		{
#ifdef AS_64BIT_PTR
			// TODO: Use a 64bit integer comparison instead of double
			ctx->bc.InstrW_W(BC_CMPd, b, c);
#else
			ctx->bc.InstrW_W(BC_CMPi, b, c);
#endif
			ctx->bc.Instr(BC_TZ);
			ctx->bc.InstrSHORT(BC_CpyRtoV4, (short)a);
		}
		else if( op == ttNotEqual )
		{
#ifdef AS_64BIT_PTR
			// TODO: Use a 64bit integer comparison instead of double
			ctx->bc.InstrW_W(BC_CMPd, b, c);
#else
			ctx->bc.InstrW_W(BC_CMPi, b, c);
#endif
			ctx->bc.Instr(BC_TNZ);
			ctx->bc.InstrSHORT(BC_CpyRtoV4, (short)a);
		}

		ctx->type.SetVariable(asCDataType::CreatePrimitive(ttBool, true), a, true);

		ReleaseTemporaryVariable(lctx->type, &ctx->bc);
		ReleaseTemporaryVariable(rctx->type, &ctx->bc);
	}
	else
	{
		// TODO: Use TXT_ILLEGAL_OPERATION_ON
		Error(TXT_ILLEGAL_OPERATION, node);
	}
}


void asCCompiler::PerformFunctionCall(int funcID, asSExprContext *ctx, bool isConstructor, asCArray<asSExprContext*> *args, asCObjectType *objType)
{
	asCScriptFunction *descr = builder->GetFunctionDescription(funcID);

	int argSize = descr->GetSpaceNeededForArguments();

	ctx->type.Set(descr->returnType);

	if( isConstructor )
		ctx->bc.Alloc(BC_ALLOC, builder->module->RefObjectType(objType), descr->id, argSize+PTR_SIZE);
	else if( descr->funcType == asFUNC_IMPORTED )
		ctx->bc.Call(BC_CALLBND , descr->id, argSize + (descr->objectType ? PTR_SIZE : 0));
	else if( descr->funcType == asFUNC_INTERFACE )
		ctx->bc.Call(BC_CALLINTF, descr->id, argSize + (descr->objectType ? PTR_SIZE : 0));
	else if( descr->funcType == asFUNC_SCRIPT )
		ctx->bc.Call(BC_CALL    , descr->id, argSize + (descr->objectType ? PTR_SIZE : 0));
	else // if( descr->funcType == asFUNC_SYSTEM )
	{
		ctx->bc.Call(BC_CALLSYS , descr->id, argSize + (descr->objectType ? PTR_SIZE : 0));
		// Add a reference for the configuration group here
		builder->module->RefConfigGroupForFunction(descr->id);
	}

	if( ctx->type.dataType.IsObject() && !descr->returnType.IsReference() )
	{
		// Allocate a temporary variable for the returned object
		int returnOffset = AllocateVariable(descr->returnType, true);
		ctx->type.SetVariable(descr->returnType, returnOffset, true);
		ctx->type.dataType.MakeReference(true);

		// Move the pointer from the object register to the temporary variable
		ctx->bc.InstrSHORT(BC_STOREOBJ, (short)returnOffset);

		// Clean up arguments
		if( args )
			AfterFunctionCall(funcID, *args, ctx, false);

		ProcessDeferredParams(ctx);

		ctx->bc.InstrSHORT(BC_PSF, (short)returnOffset);
	}
	else if( descr->returnType.IsReference() )
	{
		// We cannot clean up the arguments yet, because the
		// reference might be pointing to one of them.

		// Clean up arguments
		if( args )
			AfterFunctionCall(funcID, *args, ctx, true);

		// Do not process the output parameters yet, because it
		// might invalidate the returned reference

		if( descr->returnType.IsPrimitive() )
			ctx->type.Set(descr->returnType);
		else
		{
			ctx->bc.Instr(BC_PshRPtr);
			if( descr->returnType.IsObject() )
			{
				// We are getting the pointer to the object
				// not a pointer to a object variable
				ctx->type.dataType.MakeReference(false);
			}
		}
	}
	else
	{
		if( descr->returnType.GetSizeInMemoryBytes() )
		{
			int offset = AllocateVariable(descr->returnType, true);

			ctx->type.SetVariable(descr->returnType, offset, true);

			// Move the value from the return register to the variable
			if( descr->returnType.GetSizeOnStackDWords() == 1 )
				ctx->bc.InstrSHORT(BC_CpyRtoV4, (short)offset);
			else if( descr->returnType.GetSizeOnStackDWords() == 2 )
				ctx->bc.InstrSHORT(BC_CpyRtoV8, (short)offset);
		}
		else
			ctx->type.Set(descr->returnType);

		// Clean up arguments
		if( args )
			AfterFunctionCall(funcID, *args, ctx, false);

		ProcessDeferredParams(ctx);
	}
}


void asCCompiler::MergeExprContexts(asSExprContext *before, asSExprContext *after)
{
	before->bc.AddCode(&after->bc);

	for( asUINT n = 0; n < after->deferredParams.GetLength(); n++ )
		before->deferredParams.PushLast(after->deferredParams[n]);

	after->deferredParams.SetLength(0);

	assert( after->origExpr == 0 );
}

void asCCompiler::FilterConst(asCArray<int> &funcs)
{
	if( funcs.GetLength() == 0 ) return;

	// This is only done for object methods
	asCScriptFunction *desc = builder->GetFunctionDescription(funcs[0]);
	if( desc->objectType == 0 ) return;

	// Check if there are any non-const matches
	asUINT n;
	bool foundNonConst = false;
	for( n = 0; n < funcs.GetLength(); n++ )
	{
		desc = builder->GetFunctionDescription(funcs[n]);
		if( !desc->isReadOnly )
		{
			foundNonConst = true;
			break;
		}
	}

	if( foundNonConst )
	{
		// Remove all const methods
		for( n = 0; n < funcs.GetLength(); n++ )
		{
			desc = builder->GetFunctionDescription(funcs[n]);
			if( desc->isReadOnly )
			{
				if( n == funcs.GetLength() - 1 )
					funcs.PopLast();
				else
					funcs[n] = funcs.PopLast();

				n--;
			}
		}
	}
}

END_AS_NAMESPACE



