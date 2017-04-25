#include "rosic_asSineOscillatorStereo.h"
using namespace rosic;

asSineOscillatorStereo::asSineOscillatorStereo() : refCount(1)
{

}

asSineOscillatorStereo::~asSineOscillatorStereo()
{
  if( refCount != 0 )
    DEBUG_BREAK;
}

/*
static void SineOscillatorStereo(asSineOscillatorStereo *thisPointer) 
{
  new(thisPointer) asSineOscillatorStereo();
}
*/


void asSineOscillatorStereo::addRef()
{
  refCount++;
}

void asSineOscillatorStereo::release()
{
  if( --refCount == 0 )  
    delete this;
}

// This is where we register the type
void asSineOscillatorStereo::registerObjectType(asIScriptEngine *engine)
{
  int r;
  r = engine->RegisterObjectType("SineOscillatorStereo", sizeof(asSineOscillatorStereo), asOBJ_CLASS); rassert( r >= 0 );
  //r = engine->RegisterObjectType("SineOscillatorStereo", sizeof(asSineOscillatorStereo), asOBJ_CLASS_C); rassert( r >= 0 ); //asOBJ_CLASS_C: class has constructor

  r = engine->RegisterObjectMethod("SineOscillatorStereo", "double getSample()", asMETHOD(asSineOscillatorStereo,getSample), asCALL_THISCALL); rassert( r >= 0 );

  r = engine->RegisterObjectBehaviour("SineOscillatorStereo", asBEHAVE_ADDREF,  "void f()", asMETHOD(asSineOscillatorStereo,addRef),  asCALL_THISCALL); rassert( r >= 0 );
  r = engine->RegisterObjectBehaviour("SineOscillatorStereo", asBEHAVE_RELEASE, "void f()", asMETHOD(asSineOscillatorStereo,release), asCALL_THISCALL); rassert( r >= 0 );


  //r = engine->RegisterObjectBehaviour("SineOscillatorStereo", asBEHAVE_CONSTRUCT,  "void f()",    asFUNCTION(ConstructSineOscillatorStereo), asCALL_CDECL_OBJLAST);
  /*
  r = engine->RegisterObjectBehaviour("SineOscillatorStereo", asBEHAVE_CONSTRUCT,  "void f(const File& in)",      asFUNCTION(ConstructSineOscillatorStereoByFile), asCALL_CDECL_OBJLAST); assert( r >= 0 );
  r = engine->RegisterObjectBehaviour("SineOscillatorStereo", asBEHAVE_ADDREF,     "void f()",                    asMETHOD(asSineOscillatorStereo,addRef), asCALL_THISCALL); assert( r >= 0 );
  r = engine->RegisterObjectBehaviour("SineOscillatorStereo", asBEHAVE_RELEASE,    "void f()",                    asMETHOD(asSineOscillatorStereo,release), asCALL_THISCALL); assert( r >= 0 );
  r = engine->RegisterObjectMethod("SineOscillatorStereo", "@ getDocumentElement()",                    asMETHODPR(asSineOscillatorStereo,getDocumentElement,(),asXmlElement*), asCALL_THISCALL); assert( r >= 0 );
  r = engine->RegisterObjectMethod("SineOscillatorStereo", "XmlElement@ getDocumentElement(const bool)",          asMETHODPR(asSineOscillatorStereo,getDocumentElement,(const bool),asXmlElement*), asCALL_THISCALL); assert( r >= 0 );
  r = engine->RegisterObjectMethod("SineOscillatorStereo", "String@ getLastParseError()",                         asMETHOD(asSineOscillatorStereo,getLastParseError), asCALL_THISCALL); assert( r >= 0 );
  */
}

void asSineOscillatorStereo::setSampleRate(double newSampleRate)
{

}

void asSineOscillatorStereo::setFrequency(double newFrequency)
{

}

void asSineOscillatorStereo::setAmplitude(double newAmplitude)
{

}

void asSineOscillatorStereo::trigger()
{

}

double asSineOscillatorStereo::getSample()
{
  double left, right;
  sineOscillatorStereo.getSampleFrameStereo(&left, &right);
  return left;
}

void asSineOscillatorStereo::getSampleFrameStereo(double *outL, double *outR)
{

  sineOscillatorStereo.getSampleFrameStereo(outL, outR);
}