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
// as_config.h
//
// this file is used for configuring the compilation of the library
//

#ifndef AS_CONFIG_H
#define AS_CONFIG_H



//
// Features
//-----------------------------------------

// USE_THREADS
// Adds basic support for multithreading.
// This is currently only supported on Win32 platforms.

// BUILD_WITHOUT_LINE_CUES
// This flag makes the script compiler remove some extra bytecodes that is used
// to allow the VM to call the line callback after each statement. The performance
// is improved slightly but the scripts are only guaranteed to allow one suspension
// per loop iteration, not one per statement.

// AS_DEBUG
// This flag can be defined to make the library write some extra output when
// compiling and executing scripts.

// AS_DEPRECATED
// If this flag is defined then some backwards compatibility is maintained.
// There is no guarantee for how well deprecated functionality will work though
// so it is best to exchange it for the new functionality as soon as possible.

// AS_C_INTERFACE
// Make the C interface available.

// AS_NO_CLASS_METHODS
// Disables the possibility to add class methods. Can increase the
// portability of the library.

// AS_MAX_PORTABILITY
// Disables all platform specific code. Only the asCALL_GENERIC calling
// convention will be available in with this flag set.

// AS_NO_USER_ALLOC
// With this macro defined, the library will not use the overrideable memory
// allocation functions. Some compilers may not be compatible with this yet,
// so defining this macro may allow usage of the library on those. The macro
// will however criple the asSetGlobalMemoryFunctions so you won't get the same
// result.




//
// Library usage
//------------------------------------------

// ANGELSCRIPT_EXPORT
// This flag should be defined when compiling the library as a lib or dll.

// ANGELSCRIPT_DLL_LIBRARY_IMPORT
// This flag should be defined when using AngelScript as a dll with automatic
// library import.

// ANGELSCRIPT_DLL_MANUAL_IMPORT
// This flag should be defined when using AngelScript as a dll with manual
// loading of the library.




//
// Compiler differences
//-----------------------------------------

// vsnprintf()
// Some compilers use different names for this function. If your compiler
// doesn't use the name vsnprintf() then you need to write a macro to translate
// the function into its real name.

// ASM_AT_N_T or ASM_INTEL
// You should choose what inline assembly syntax to use when compiling.

// VALUE_OF_BOOLEAN_TRUE
// This flag allows to customize the exact value of boolean true.

// AS_SIZEOF_BOOL
// On some target platforms the sizeof(bool) is 4, but on most it is 1.

// STDCALL
// This is used to declare a function to use the stdcall calling convention.

// AS_USE_NAMESPACE
// Adds the AngelScript namespace on the declarations.

// AS_NO_MEMORY_H
// Some compilers don't come with the memory.h header file.



//
// How to identify different compilers
//-----------------------------------------

// MS Visual C++
//  _MSC_VER   is defined
//  __MWERKS__ is not defined

// Metrowerks
//  _MSC_VER   is defined
//  __MWERKS__ is defined

// GNU C based compilers
//  __GNUC__   is defined



//
// CPU differences
//---------------------------------------

// AS_ALIGN
// Some CPUs require that data words are aligned in some way. This macro
// should be defined if the words should be aligned to boundaries of the same
// size as the word, i.e.
//  1 byte  on 1 byte boundaries
//  2 bytes on 2 byte boundaries
//  4 bytes on 4 byte boundaries
//  8 bytes on 4 byte boundaries (no it's not a typo)

// AS_USE_DOUBLE_AS_FLOAT
// If there is no 64 bit floating point type, then this constant can be defined
// to treat double like normal floats.

// AS_X86
// Use assembler code for the x86 CPU family

// AS_SH4
// Use assembler code for the SH4 CPU family

// AS_MIPS
// Use assembler code for the MIPS CPU family

// AS_PPC
// Use assembler code for the 32bit PowerPC CPU family

// AS_PPC_64
// Use assembler code for the 64bit PowerPC CPU family

// AS_64BIT_PTR
// Define this to make the engine store all pointers in 64bit words.

// AS_BIG_ENDIAN
// Define this for CPUs that use big endian memory layout, e.g. PPC



// 
// Target systems
//--------------------------------
// This group shows a few of the flags used to identify different target systems. 
// Sometimes there are differences on different target systems, while both CPU and 
// compiler is the same for both, when this is so these flags are used to produce the
// right code.

// AS_WIN     - Microsoft Windows
// AS_LINUX   - Linux
// AS_MAC     - Apple Macintosh
// AS_XBOX    - Microsoft XBox 
// AS_XBOX360 - Microsoft XBox 360
// AS_PSP     - Sony Playstation Portable
// AS_PS2     - Sony Playstation 2
// AS_PS3     - Sony Playstation 3
// AS_DC      - Sega Dreamcast




//
// Calling conventions
//-----------------------------------------

// GNU_STYLE_VIRTUAL_METHOD
// This constant should be defined if method pointers store index for virtual
// functions in the same location as the function pointer. In such cases the method
// is identified as virtual if the least significant bit is set.

// MULTI_BASE_OFFSET(x)
// This macro is used to retrieve the offset added to the object pointer in order to
// implicitly cast the object to the base object. x is the method pointer received by
// the register function.

// HAVE_VIRTUAL_BASE_OFFSET
// Define this constant if the compiler stores the virtual base offset in the method
// pointers. If it is not stored in the pointers then AngelScript have no way of
// identifying a method as coming from a class with virtual inheritance.

// VIRTUAL_BASE_OFFSET(x)
// This macro is used to retrieve the offset added to the object pointer in order to
// find the virtual base object. x is the method pointer received by the register
// function;

// COMPLEX_MASK
// This constant shows what attributes determines if an object is returned in memory
// or in the registers as normal structures

// THISCALL_RETURN_SIMPLE_IN_MEMORY
// CDECL_RETURN_SIMPLE_IN_MEMORY
// STDCALL_RETURN_SIMPLE_IN_MEMORY
// When these constants are defined then the corresponding calling convention always
// return classes/structs in memory regardless of size or complexity.

// CALLEE_POPS_HIDDEN_RETURN_POINTER
// This constant should be defined if the callee pops the hidden return pointer,
// used when returning an object in memory.

// THISCALL_PASS_OBJECT_POINTER_ON_THE_STACK
// With this constant defined AngelScript will pass the object pointer on the stack

// THISCALL_CALLEE_POPS_ARGUMENTS
// If the callee pops arguments for class methods then define this constant

// COMPLEX_OBJS_PASSED_BY_REF
// Some compilers always pass certain objects by reference. GNUC for example does
// this if the the class has a defined destructor.



// TODO: 
// First detect hardware and OS, set AngelScript defined variables for these
// Then detect compiler and configure the library based on detected target system (using the AngelScript defined variables)




//
// Detect compiler
//------------------------------------------------

#define VALUE_OF_BOOLEAN_TRUE  1

// Microsoft Visual C++
#if defined(_MSC_VER) && !defined(__MWERKS__)
	#define MULTI_BASE_OFFSET(x) (*((asDWORD*)(&x)+1))
	#define HAVE_VIRTUAL_BASE_OFFSET
	#define VIRTUAL_BASE_OFFSET(x) (*((asDWORD*)(&x)+3))
	#define THISCALL_RETURN_SIMPLE_IN_MEMORY
	#define THISCALL_PASS_OBJECT_POINTER_IN_ECX
	//#define vsnprintf(a, b, c, d) _vsnprintf(a, b, c, d)
	#define THISCALL_CALLEE_POPS_ARGUMENTS
	#define COMPLEX_MASK (asOBJ_CLASS_CONSTRUCTOR | asOBJ_CLASS_DESTRUCTOR | asOBJ_CLASS_ASSIGNMENT)
	#define STDCALL __stdcall
	#define AS_SIZEOF_BOOL 1

	#define ASM_INTEL  // Intel style for inline assembly on microsoft compilers

	#if _XBOX_VER >= 200
		// 360 is PPC, big endian and 64 bit.
		#define AS_XBOX360
		#define AS_PPC_64
	#else
		// Support native calling conventions on x86, but not 64bit yet
		#if defined(_XBOX) || (defined(_M_IX86) && !defined(__LP64__))
			#define AS_X86
		#endif
	#endif

	#if _MSC_VER <= 1200 // MSVC++ 6
		#define I64(x) x##l
	#else
		#define I64(x) x##ll
	#endif

	#define UNREACHABLE_RETURN
#endif

// Metrowerks CodeWarrior (experimental, let me know if something isn't working)
#if defined(__MWERKS__)
	#define MULTI_BASE_OFFSET(x) (*((asDWORD*)(&x)+1))
	#define HAVE_VIRTUAL_BASE_OFFSET
	#define VIRTUAL_BASE_OFFSET(x) (*((asDWORD*)(&x)+3))
	#define THISCALL_RETURN_SIMPLE_IN_MEMORY
	#define THISCALL_PASS_OBJECT_POINTER_IN_ECX
	#define vsnprintf(a, b, c, d) _vsnprintf(a, b, c, d)
	#define THISCALL_CALLEE_POPS_ARGUMENTS
	#define COMPLEX_MASK (asOBJ_CLASS_CONSTRUCTOR | asOBJ_CLASS_DESTRUCTOR | asOBJ_CLASS_ASSIGNMENT)
	#define AS_SIZEOF_BOOL 1
	#define STDCALL __stdcall

	// Support native calling conventions on x86, but not 64bit yet
	#if defined(_M_IX86) && !defined(__LP64__)
		#define AS_X86
		#define ASM_INTEL  // Intel style for inline assembly
	#endif

	#if _MSC_VER <= 1200 // MSVC++ 6
		#define I64(x) x##l
	#else
		#define I64(x) x##ll
	#endif

	#define UNREACHABLE_RETURN
#endif

// SN Systems ProDG (also experimental, let me know if something isn't working)
#if defined(__SNC__) || defined(SNSYS)
	#define GNU_STYLE_VIRTUAL_METHOD
	#define MULTI_BASE_OFFSET(x) (*((asDWORD*)(&x)+1))
	#define CALLEE_POPS_HIDDEN_RETURN_POINTER
	#define COMPLEX_OBJS_PASSED_BY_REF
	#define ASM_AT_N_T  // AT&T style inline assembly
	#define COMPLEX_MASK (asOBJ_CLASS_DESTRUCTOR)
	#define AS_SIZEOF_BOOL 1

	// SN doesnt seem to like STDCALL.
	// Maybe it can work with some fiddling, but I can't imagine linking to 
	// any STDCALL functions with a console anyway...
	#define STDCALL

	// Linux specific
	#ifdef __linux__
		#define THISCALL_RETURN_SIMPLE_IN_MEMORY
		#define CDECL_RETURN_SIMPLE_IN_MEMORY
		#define STDCALL_RETURN_SIMPLE_IN_MEMORY
	#endif

	// Support native calling conventions on x86, but not 64bit yet
	#if defined(i386) && !defined(__LP64__)
		#define AS_X86
	#endif

	#define I64(x) x##ll

	#define UNREACHABLE_RETURN
#endif

// GNU C
#if defined(__GNUC__) && !defined(__SNC__)
	#define GNU_STYLE_VIRTUAL_METHOD
	#define MULTI_BASE_OFFSET(x) (*((asDWORD*)(&x)+1))
	#define CALLEE_POPS_HIDDEN_RETURN_POINTER
	#define COMPLEX_OBJS_PASSED_BY_REF
	#define COMPLEX_MASK (asOBJ_CLASS_DESTRUCTOR)
	#define AS_NO_MEMORY_H
	#define AS_SIZEOF_BOOL 1
	#define STDCALL __attribute__((stdcall))
	#define ASM_AT_N_T

	// MacOSX
	#ifdef __APPLE__
		// The sizeof bool is different depending on the target CPU
		#define AS_MAC
		#undef AS_SIZEOF_BOOL
		#if defined(__ppc__)
			#define AS_SIZEOF_BOOL 4
			// STDCALL is not available on PPC
			#undef STDCALL
			#define STDCALL
		#else
			#define AS_SIZEOF_BOOL 1
		#endif
		#if defined(i386) && !defined(__LP64__)
			// Support native calling conventions on Mac OS X + Intel 32bit CPU
			#define AS_X86
		#elif (defined(__ppc__) || defined(__PPC__)) && !defined(__LP64__)
			// Support native calling conventions on Mac OS X + PPC 32bit CPU
			#define AS_PPC
			#define THISCALL_RETURN_SIMPLE_IN_MEMORY
			#define CDECL_RETURN_SIMPLE_IN_MEMORY
			#define STDCALL_RETURN_SIMPLE_IN_MEMORY
		#elif (defined(__ppc__) || defined(__PPC__)) && defined(__LP64__)
			#define AS_PPC_64
		#else
			// No support for native calling conventions yet
			#define AS_MAX_PORTABILITY
		#endif

	// Windows and Linux
	#elif defined(WIN32) || defined(__linux__)
		#define THISCALL_RETURN_SIMPLE_IN_MEMORY
		#define CDECL_RETURN_SIMPLE_IN_MEMORY
		#define STDCALL_RETURN_SIMPLE_IN_MEMORY
		#if defined(i386) && !defined(__LP64__)
			// Support native calling conventions on Intel 32bit CPU
			#define AS_X86
		#else
			// No support for native calling conventions yet
			#define AS_MAX_PORTABILITY
		#endif
		
	// PSP and PS2
	#elif defined(__PSP__) || defined(__psp__) || defined(_EE_) || defined(_PSP) || defined(_PS2)
		// Support native calling conventions on MIPS architecture
		#if (defined(_MIPS_ARCH) || defined(_mips) || defined(__MIPSEL__)) && !defined(__LP64__)
			#define AS_MIPS
		#else
			#define AS_MAX_PORTABILITY
		#endif
		
	// PS3
	#elif (defined(__PPC__) || defined(__ppc__)) && defined(__PPU__)
		// Support native calling conventions on PS3
		#define AS_PS3
		#define AS_PPC_64
		#define THISCALL_RETURN_SIMPLE_IN_MEMORY
		#define CDECL_RETURN_SIMPLE_IN_MEMORY
		#define STDCALL_RETURN_SIMPLE_IN_MEMORY
		// PS3 doesn't have STDCALL
		#undef STDCALL
		#define STDCALL

	// Dreamcast
	#elif __SH4_SINGLE_ONLY__
		// Support native calling conventions on Dreamcast
		#define AS_DC
		#define AS_SH4
	#endif

	#define I64(x) x##ll

	#define UNREACHABLE_RETURN
#endif


//
// Detect target hardware
//------------------------------------------------

// X86, Intel, AMD, etc, i.e. most PCs
#if defined(__i386__) || defined(_M_IX86)
	// Nothing special here
#endif

// MIPS architecture (generally PS2 and PSP consoles, potentially supports N64 as well)
#if defined(_MIPS_ARCH) || defined(_mips) || defined(__MIPSEL__) || defined(__PSP__) || defined(__psp__) || defined(_EE_) || defined(_PSP) || defined(_PS2)
	#define AS_ALIGN				// align datastructures
	#define AS_USE_DOUBLE_AS_FLOAT	// use 32bit floats instead of doubles
#endif

// PowerPC, e.g. Mac, GameCube, XBox 360 and PS3
#if defined(__PPC__) || defined(__ppc__)
	#define AS_BIG_ENDIAN

	// Gamecube
	#if defined(_GC)
		#define AS_ALIGN
		#define AS_USE_DOUBLE_AS_FLOAT
	#endif
	// XBox 360
	#if (_XBOX_VER >= 200 )
		#define AS_ALIGN
		#define AS_USE_DOUBLE_AS_FLOAT
	#endif
	// PS3
	#if defined(__PPU__)
		#define AS_ALIGN
	#endif
#endif

// Dreamcast console
#ifdef __SH4_SINGLE_ONLY__
	#define AS_ALIGN				// align datastructures
	#define AS_USE_DOUBLE_AS_FLOAT	// use 32bit floats instead of doubles
#endif

// Is the target a 64bit system?
#if defined(__LP64__) || defined(__amd64__)
	#ifndef AS_64BIT_PTR
		#define AS_64BIT_PTR
	#endif
#endif

// If there are no current support for native calling
// conventions, then compile with AS_MAX_PORTABILITY
#if (!defined(AS_X86) && !defined(AS_SH4) && !defined(AS_MIPS) && !defined(AS_PPC) && !defined(AS_PPC_64))
	#ifndef AS_MAX_PORTABILITY
		#define AS_MAX_PORTABILITY
	#endif
#endif





//
// Internal defines (do not change these)
//----------------------------------------------------------------

#ifdef AS_ALIGN
	#define	ALIGN(b) (((b)+3)&(~3))
#else
	#define	ALIGN(b) (b)
#endif

#define	ARG_W(b)    ((asWORD*)&b)
#define	ARG_DW(b)   ((asDWORD*)&b)
#define	ARG_QW(b)   ((asQWORD*)&b)
#define	BCARG_W(b)  ((asWORD*)&(b)[1])
#define	BCARG_DW(b) ((asDWORD*)&(b)[1])
#define	BCARG_QW(b) ((asQWORD*)&(b)[1])

#ifdef AS_64BIT_PTR
	#define PTR_SIZE     2
	#define asPTRWORD    asQWORD
#else
	#define PTR_SIZE     1
	#define asPTRWORD    asDWORD
#endif
#define ARG_PTR(b)   ((asPTRWORD*)&b)
#define BCARG_PTR(b) ((asPTRWORD*)&(b)[1])

// This macro is used to avoid warnings about unused variables.
// Usually where the variables are only used in debug mode.
#define UNUSED_VAR(x) (x)=(x)


#include "angelscript.h"
#include "as_memory.h"

#ifdef AS_USE_NAMESPACE
using namespace AngelScript;
#endif

#endif


