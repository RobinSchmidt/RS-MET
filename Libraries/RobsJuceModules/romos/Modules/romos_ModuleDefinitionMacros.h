#ifndef romos_ModuleDefinitionMacros_h
#define romos_ModuleDefinitionMacros_h

/** This file contains preprocessor macros that facilitate the definition a new module class by 
generating much of the boilerplate code that otherwise would have to written for each module.

todo: Maybe instead of using the pre-processor, write the code generation functions directly in
C++. This has the advantage that we can actually inspect the generated code and set debug 
breakpoints. When new atomic modules are added, we just need to run the code-generator once. We can
create a special project for that. Maybe it can be made part of the Liberty test project. */





// This macro should be used in each Module subclass declaration in order to make the constructor, 
// copy-constructor, assignment operator and destructor protected and to declare the ModuleFactor 
// as friend in order to enforce creation/deletion through the ModuleFactory class - maybe make 
// them even private later:
#define ENFORCE_FACTORY_USAGE(ClassName)                      \
  public:                                                     \
    ClassName() {}                                            \
  protected:                                                  \
    ClassName(const ClassName&) {}                            \
    virtual ~ClassName() {}                                   \
    ClassName& operator=(const ClassName&) { return *this; }  \
    friend class ModuleFactory;                               \

/* // old - doesn't work with new type registry:
#define ENFORCE_FACTORY_USAGE(ClassName)                               \
  protected:                                                           \
    ClassName() {}                                                     \
    ClassName(const ClassName &other) {}                               \
    virtual ~ClassName() {}                                            \
    ClassName& operator=(const ClassName &other) { return *this; }     \
    friend class ModuleFactory;                                        \
*/


// This macro can be used in a Module subclass to declare all 4 static processing functions and 
// the (then necessary) override for assignProcessingFunctions():
#define DECLARE_PROCESSING_FUNCTIONS                                                      \
  protected:                                                                              \
    static INLINE void processMonoFrame(Module *module, int voiceIndex);                  \
    static INLINE void processPolyFrame(Module *module, int voiceIndex);                  \
    static INLINE void processMonoBlock(Module *module, int voiceIndex,  int blockSize);  \
    static INLINE void processPolyBlock(Module *module, int voiceIndex,  int blockSize);  \
    virtual void assignProcessingFunctions();                                             \


// This macro can be used in place of ENFORCE_FACTORY_USAGE, DECLARE_PROCESSING_FUNCTIONS and the 
// declaration of the override of initialize for Module subclasses where all 3 are needed (which is
// the common case):
#define CREATE_COMMON_DECLARATIONS(ClassName)   \
  ENFORCE_FACTORY_USAGE(ClassName);             \
  DECLARE_PROCESSING_FUNCTIONS;                 \
  virtual void initialize();                    \






// creates the assignFunctionPointers() member function for the class with the given name:
#define CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)       \
  void ClassName ::assignProcessingFunctions()           \
  {                                                      \
    if( isPolyphonic() )                                 \
    {                                                    \
      processFrame = & ClassName ::processPolyFrame;     \
      processBlock = & ClassName ::processPolyBlock;     \
    }                                                    \
    else                                                 \
    {                                                    \
      processFrame = & ClassName ::processMonoFrame;     \
      processBlock = & ClassName ::processMonoBlock;     \
    }                                                    \
  }                                                      \


// given a function process(Module *module, double *out, int voiceIndex), these macros create the 
// mono/poly, frame/block processing functions respectively for modules without input pins:
#define CREATE_MONO_FRAME_FUNCTION_0(ClassName)                                                \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                             \
  {                                                                                            \
    ClassName::process(module,                                                                 \
                       module->audioOutputs,                                                   \
                       voiceIndex);                                                            \
  }                                                                                            \

#define CREATE_POLY_FRAME_FUNCTION_0(ClassName)                                                \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                             \
  {                                                                                            \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();            \
    double *outPointer = module->audioOutputs;                                                 \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)      \
    {                                                                                          \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                    \
      ClassName::process(module,                                                               \
                         outPointer + voiceIndex * outVoiceStride,                             \
                         voiceIndex);                                                          \
    }                                                                                          \
  }                                                                                            \

#define CREATE_MONO_BLOCK_FUNCTION_0(ClassName)                                                \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)              \
  {                                                                                            \
    int outFrameStride = module->outFrameStride;                                               \
    double *outPointer = module->audioOutputs;                                                 \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                              \
    {                                                                                          \
      ClassName::process(module, outPointer, voiceIndex);                                      \
      outPointer += outFrameStride;                                                            \
    }                                                                                          \
  }                                                                                            \

#define CREATE_POLY_BLOCK_FUNCTION_0(ClassName)                                                \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)              \
  {                                                                                            \
    int outFrameStride       = module->outFrameStride;                                         \
    int outVoiceStride       = outFrameStride * processingStatus.getBufferSize();              \
    double *outPointerVoice0 = module->audioOutputs;                                           \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)      \
    {                                                                                          \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                 \
      double *outPointer = outPointerVoice0 + voiceIndex * outVoiceStride;                     \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                            \
      {                                                                                        \
        ClassName::process(module, outPointer, voiceIndex);                                    \
        outPointer += outFrameStride;                                                          \
      }                                                                                        \
    }                                                                                          \
  }                                                                                            \

#define CREATE_DERIVED_FUNCTIONS_0(ClassName)                                                  \
  CREATE_MONO_FRAME_FUNCTION_0(ClassName);                                                     \
  CREATE_POLY_FRAME_FUNCTION_0(ClassName);                                                     \
  CREATE_MONO_BLOCK_FUNCTION_0(ClassName);                                                     \
  CREATE_POLY_BLOCK_FUNCTION_0(ClassName);                                                     \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_0(ClassName)                                    \
  CREATE_DERIVED_FUNCTIONS_0(ClassName)                                                        \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                   \

#define CREATE_COMMON_DECLARATIONS_0(ClassName)                                                \
  ENFORCE_FACTORY_USAGE(ClassName);                                                            \
  DECLARE_PROCESSING_FUNCTIONS;                                                                \
  virtual void initialize();                                                                   \
  static INLINE void process(Module *module, double *out, int voiceIndex);                     \


// given a function process(Module *module, double *in, double *out, int voiceIndex), these macros 
// create the mono/poly, frame/block processing functions respectively for modules with one input 
// pin:
#define CREATE_MONO_FRAME_FUNCTION_1(ClassName)                                                \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                             \
  {                                                                                            \
    ClassName::process(module,                                                                 \
                       module->inputPins[0].outputPointer,                                     \
                       module->audioOutputs,                                                   \
                       voiceIndex);                                                            \
  }                                                                                            \

#define CREATE_POLY_FRAME_FUNCTION_1(ClassName)                                                \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                             \
  {                                                                                            \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();            \
    int inVoiceStride  = module->inputPins[0].outputVoiceStride;                               \
    double *inPointer  = module->inputPins[0].outputPointer;                                   \
    double *outPointer = module->audioOutputs;                                                 \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)      \
    {                                                                                          \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                    \
      ClassName::process(module,                                                               \
                         inPointer  + voiceIndex * inVoiceStride,                              \
                         outPointer + voiceIndex * outVoiceStride,                             \
                         voiceIndex);                                                          \
    }                                                                                          \
  }                                                                                            \

#define CREATE_MONO_BLOCK_FUNCTION_1(ClassName)                                                \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)              \
  {                                                                                            \
    int outFrameStride = module->outFrameStride;                                               \
    int inFrameStride  = module->inputPins[0].outputFrameSize;                                 \
    double *inPointer  = module->inputPins[0].outputPointer;                                   \
    double *outPointer = module->audioOutputs;                                                 \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                              \
    {                                                                                          \
      ClassName::process(module, inPointer, outPointer, voiceIndex);                           \
      inPointer  += inFrameStride;                                                             \
      outPointer += outFrameStride;                                                            \
    }                                                                                          \
  }                                                                                            \

#define CREATE_POLY_BLOCK_FUNCTION_1(ClassName)                                                \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)              \
  {                                                                                            \
    int outFrameStride       = module->outFrameStride;                                         \
    int outVoiceStride       = outFrameStride * processingStatus.getBufferSize();              \
    int inFrameStride        = module->inputPins[0].outputFrameSize;                           \
    int inVoiceStride        = module->inputPins[0].outputVoiceStride;                         \
    double *inPointerVoice0  = module->inputPins[0].outputPointer;                             \
    double *outPointerVoice0 = module->audioOutputs;                                           \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)      \
    {                                                                                          \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                 \
      double *inPointer  = inPointerVoice0  + voiceIndex * inVoiceStride;                      \
      double *outPointer = outPointerVoice0 + voiceIndex * outVoiceStride;                     \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                            \
      {                                                                                        \
        ClassName::process(module, inPointer, outPointer, voiceIndex);                         \
        inPointer  += inFrameStride;                                                           \
        outPointer += outFrameStride;                                                          \
      }                                                                                        \
    }                                                                                          \
  }                                                                                            \

#define CREATE_DERIVED_FUNCTIONS_1(ClassName)                                                  \
  CREATE_MONO_FRAME_FUNCTION_1(ClassName);                                                     \
  CREATE_POLY_FRAME_FUNCTION_1(ClassName);                                                     \
  CREATE_MONO_BLOCK_FUNCTION_1(ClassName);                                                     \
  CREATE_POLY_BLOCK_FUNCTION_1(ClassName);                                                     \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(ClassName)                                    \
  CREATE_DERIVED_FUNCTIONS_1(ClassName)                                                        \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                   \

#define CREATE_COMMON_DECLARATIONS_1(ClassName)                                                \
  ENFORCE_FACTORY_USAGE(ClassName);                                                            \
  DECLARE_PROCESSING_FUNCTIONS;                                                                \
  virtual void initialize();                                                                   \
  static INLINE void process(Module *module, double *in, double *out, int voiceIndex);         \



// given a function process(Module *module, double *in1, double *in2, double *out, int voiceIndex), 
// these macros create the mono/poly, frame/block processing functions respectively for modules 
// with two input pins:
#define CREATE_MONO_FRAME_FUNCTION_2(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_2(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_2(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module, inPointer0, inPointer1, outPointer, voiceIndex);                     \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_2(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module, inPointer0, inPointer1, outPointer, voiceIndex);                   \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_2(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_2(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_2(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_2(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_2(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_2(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_2(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_2(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *out, int voiceIndex);  \


// given a function process(Module *module, double *in1, double *in2, double *in3, double *out, int voiceIndex), these macros create the 
// mono/poly, frame/block processing functions respectively for modules with 3 input pins:
#define CREATE_MONO_FRAME_FUNCTION_3(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->inputPins[2].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_3(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    int inVoiceStride2 = module->inputPins[2].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         inPointer2 + voiceIndex * inVoiceStride2,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_3(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    int inFrameStride2 = module->inputPins[2].outputFrameSize;                                        \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module,                                                                      \
                         inPointer0,                                                                  \
                         inPointer1,                                                                  \
                         inPointer2,                                                                  \
                         outPointer,                                                                  \
                         voiceIndex);                                                                 \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      inPointer2 += inFrameStride2;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_3(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    int inFrameStride2          = module->inputPins[2].outputFrameSize;                               \
    int inVoiceStride2          = module->inputPins[2].outputVoiceStride;                             \
    double *inPointerVoice0Pin2 = module->inputPins[2].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *inPointer2 = inPointerVoice0Pin2 + voiceIndex * inVoiceStride2;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module,                                                                    \
                           inPointer0,                                                                \
                           inPointer1,                                                                \
                           inPointer2,                                                                \
                           outPointer,                                                                \
                           voiceIndex);                                                               \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        inPointer2 += inFrameStride2;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_3(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_3(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_3(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_3(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_3(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_3(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_3(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *in3, double *out,      \
                             int voiceIndex);                                                         \


// ...and so on - 4 input pins:
#define CREATE_MONO_FRAME_FUNCTION_4(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->inputPins[2].outputPointer,                                            \
                       module->inputPins[3].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_4(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    int inVoiceStride2 = module->inputPins[2].outputVoiceStride;                                      \
    int inVoiceStride3 = module->inputPins[3].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         inPointer2 + voiceIndex * inVoiceStride2,                                    \
                         inPointer3 + voiceIndex * inVoiceStride3,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_4(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    int inFrameStride2 = module->inputPins[2].outputFrameSize;                                        \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    int inFrameStride3 = module->inputPins[3].outputFrameSize;                                        \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module,                                                                      \
                         inPointer0,                                                                  \
                         inPointer1,                                                                  \
                         inPointer2,                                                                  \
                         inPointer3,                                                                  \
                         outPointer,                                                                  \
                         voiceIndex);                                                                 \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      inPointer2 += inFrameStride2;                                                                   \
      inPointer3 += inFrameStride3;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_4(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    int inFrameStride2          = module->inputPins[2].outputFrameSize;                               \
    int inVoiceStride2          = module->inputPins[2].outputVoiceStride;                             \
    double *inPointerVoice0Pin2 = module->inputPins[2].outputPointer;                                 \
    int inFrameStride3          = module->inputPins[3].outputFrameSize;                               \
    int inVoiceStride3          = module->inputPins[3].outputVoiceStride;                             \
    double *inPointerVoice0Pin3 = module->inputPins[3].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *inPointer2 = inPointerVoice0Pin2 + voiceIndex * inVoiceStride2;                         \
      double *inPointer3 = inPointerVoice0Pin3 + voiceIndex * inVoiceStride3;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module,                                                                    \
                           inPointer0,                                                                \
                           inPointer1,                                                                \
                           inPointer2,                                                                \
                           inPointer3,                                                                \
                           outPointer,                                                                \
                           voiceIndex);                                                               \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        inPointer2 += inFrameStride2;                                                                 \
        inPointer3 += inFrameStride3;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_4(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_4(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_4(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_4(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_4(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_4(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_4(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_4(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *in3, double *in4,      \
                             double *out, int voiceIndex);                                            \



// 5 input pins:
#define CREATE_MONO_FRAME_FUNCTION_5(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->inputPins[2].outputPointer,                                            \
                       module->inputPins[3].outputPointer,                                            \
                       module->inputPins[4].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_5(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    int inVoiceStride2 = module->inputPins[2].outputVoiceStride;                                      \
    int inVoiceStride3 = module->inputPins[3].outputVoiceStride;                                      \
    int inVoiceStride4 = module->inputPins[4].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         inPointer2 + voiceIndex * inVoiceStride2,                                    \
                         inPointer3 + voiceIndex * inVoiceStride3,                                    \
                         inPointer4 + voiceIndex * inVoiceStride4,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_5(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    int inFrameStride2 = module->inputPins[2].outputFrameSize;                                        \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    int inFrameStride3 = module->inputPins[3].outputFrameSize;                                        \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    int inFrameStride4 = module->inputPins[4].outputFrameSize;                                        \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module,                                                                      \
                         inPointer0,                                                                  \
                         inPointer1,                                                                  \
                         inPointer2,                                                                  \
                         inPointer3,                                                                  \
                         inPointer4,                                                                  \
                         outPointer,                                                                  \
                         voiceIndex);                                                                 \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      inPointer2 += inFrameStride2;                                                                   \
      inPointer3 += inFrameStride3;                                                                   \
      inPointer4 += inFrameStride4;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_5(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    int inFrameStride2          = module->inputPins[2].outputFrameSize;                               \
    int inVoiceStride2          = module->inputPins[2].outputVoiceStride;                             \
    double *inPointerVoice0Pin2 = module->inputPins[2].outputPointer;                                 \
    int inFrameStride3          = module->inputPins[3].outputFrameSize;                               \
    int inVoiceStride3          = module->inputPins[3].outputVoiceStride;                             \
    double *inPointerVoice0Pin3 = module->inputPins[3].outputPointer;                                 \
    int inFrameStride4          = module->inputPins[4].outputFrameSize;                               \
    int inVoiceStride4          = module->inputPins[4].outputVoiceStride;                             \
    double *inPointerVoice0Pin4 = module->inputPins[4].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *inPointer2 = inPointerVoice0Pin2 + voiceIndex * inVoiceStride2;                         \
      double *inPointer3 = inPointerVoice0Pin3 + voiceIndex * inVoiceStride3;                         \
      double *inPointer4 = inPointerVoice0Pin4 + voiceIndex * inVoiceStride4;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module,                                                                    \
                           inPointer0,                                                                \
                           inPointer1,                                                                \
                           inPointer2,                                                                \
                           inPointer3,                                                                \
                           inPointer4,                                                                \
                           outPointer,                                                                \
                           voiceIndex);                                                               \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        inPointer2 += inFrameStride2;                                                                 \
        inPointer3 += inFrameStride3;                                                                 \
        inPointer4 += inFrameStride4;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_5(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_5(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_5(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_5(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_5(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_5(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_5(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_5(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *in3, double *in4,      \
                             double *in5, double *out, int voiceIndex);                               \


// 6 input pins:
#define CREATE_MONO_FRAME_FUNCTION_6(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->inputPins[2].outputPointer,                                            \
                       module->inputPins[3].outputPointer,                                            \
                       module->inputPins[4].outputPointer,                                            \
                       module->inputPins[5].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_6(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    int inVoiceStride2 = module->inputPins[2].outputVoiceStride;                                      \
    int inVoiceStride3 = module->inputPins[3].outputVoiceStride;                                      \
    int inVoiceStride4 = module->inputPins[4].outputVoiceStride;                                      \
    int inVoiceStride5 = module->inputPins[5].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    double *inPointer5 = module->inputPins[5].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         inPointer2 + voiceIndex * inVoiceStride2,                                    \
                         inPointer3 + voiceIndex * inVoiceStride3,                                    \
                         inPointer4 + voiceIndex * inVoiceStride4,                                    \
                         inPointer5 + voiceIndex * inVoiceStride5,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_6(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    int inFrameStride2 = module->inputPins[2].outputFrameSize;                                        \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    int inFrameStride3 = module->inputPins[3].outputFrameSize;                                        \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    int inFrameStride4 = module->inputPins[4].outputFrameSize;                                        \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    int inFrameStride5 = module->inputPins[5].outputFrameSize;                                        \
    double *inPointer5 = module->inputPins[5].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module,                                                                      \
                         inPointer0,                                                                  \
                         inPointer1,                                                                  \
                         inPointer2,                                                                  \
                         inPointer3,                                                                  \
                         inPointer4,                                                                  \
                         inPointer5,                                                                  \
                         outPointer,                                                                  \
                         voiceIndex);                                                                 \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      inPointer2 += inFrameStride2;                                                                   \
      inPointer3 += inFrameStride3;                                                                   \
      inPointer4 += inFrameStride4;                                                                   \
      inPointer5 += inFrameStride5;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_6(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    int inFrameStride2          = module->inputPins[2].outputFrameSize;                               \
    int inVoiceStride2          = module->inputPins[2].outputVoiceStride;                             \
    double *inPointerVoice0Pin2 = module->inputPins[2].outputPointer;                                 \
    int inFrameStride3          = module->inputPins[3].outputFrameSize;                               \
    int inVoiceStride3          = module->inputPins[3].outputVoiceStride;                             \
    double *inPointerVoice0Pin3 = module->inputPins[3].outputPointer;                                 \
    int inFrameStride4          = module->inputPins[4].outputFrameSize;                               \
    int inVoiceStride4          = module->inputPins[4].outputVoiceStride;                             \
    double *inPointerVoice0Pin4 = module->inputPins[4].outputPointer;                                 \
    int inFrameStride5          = module->inputPins[5].outputFrameSize;                               \
    int inVoiceStride5          = module->inputPins[5].outputVoiceStride;                             \
    double *inPointerVoice0Pin5 = module->inputPins[5].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *inPointer2 = inPointerVoice0Pin2 + voiceIndex * inVoiceStride2;                         \
      double *inPointer3 = inPointerVoice0Pin3 + voiceIndex * inVoiceStride3;                         \
      double *inPointer4 = inPointerVoice0Pin4 + voiceIndex * inVoiceStride4;                         \
      double *inPointer5 = inPointerVoice0Pin5 + voiceIndex * inVoiceStride5;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module,                                                                    \
                           inPointer0,                                                                \
                           inPointer1,                                                                \
                           inPointer2,                                                                \
                           inPointer3,                                                                \
                           inPointer4,                                                                \
                           inPointer5,                                                                \
                           outPointer,                                                                \
                           voiceIndex);                                                               \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        inPointer2 += inFrameStride2;                                                                 \
        inPointer3 += inFrameStride3;                                                                 \
        inPointer4 += inFrameStride4;                                                                 \
        inPointer5 += inFrameStride5;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_6(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_6(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_6(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_6(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_6(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_6(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_6(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_6(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *in3, double *in4,      \
                             double *in5, double *in6, double *out, int voiceIndex);                  \


// 7 input pins:
#define CREATE_MONO_FRAME_FUNCTION_7(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->inputPins[2].outputPointer,                                            \
                       module->inputPins[3].outputPointer,                                            \
                       module->inputPins[4].outputPointer,                                            \
                       module->inputPins[5].outputPointer,                                            \
                       module->inputPins[6].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_7(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    int inVoiceStride2 = module->inputPins[2].outputVoiceStride;                                      \
    int inVoiceStride3 = module->inputPins[3].outputVoiceStride;                                      \
    int inVoiceStride4 = module->inputPins[4].outputVoiceStride;                                      \
    int inVoiceStride5 = module->inputPins[5].outputVoiceStride;                                      \
    int inVoiceStride6 = module->inputPins[6].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    double *inPointer5 = module->inputPins[5].outputPointer;                                          \
    double *inPointer6 = module->inputPins[6].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         inPointer2 + voiceIndex * inVoiceStride2,                                    \
                         inPointer3 + voiceIndex * inVoiceStride3,                                    \
                         inPointer4 + voiceIndex * inVoiceStride4,                                    \
                         inPointer5 + voiceIndex * inVoiceStride5,                                    \
                         inPointer6 + voiceIndex * inVoiceStride6,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_7(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    int inFrameStride2 = module->inputPins[2].outputFrameSize;                                        \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    int inFrameStride3 = module->inputPins[3].outputFrameSize;                                        \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    int inFrameStride4 = module->inputPins[4].outputFrameSize;                                        \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    int inFrameStride5 = module->inputPins[5].outputFrameSize;                                        \
    double *inPointer5 = module->inputPins[5].outputPointer;                                          \
    int inFrameStride6 = module->inputPins[6].outputFrameSize;                                        \
    double *inPointer6 = module->inputPins[6].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module,                                                                      \
                         inPointer0,                                                                  \
                         inPointer1,                                                                  \
                         inPointer2,                                                                  \
                         inPointer3,                                                                  \
                         inPointer4,                                                                  \
                         inPointer5,                                                                  \
                         inPointer6,                                                                  \
                         outPointer,                                                                  \
                         voiceIndex);                                                                 \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      inPointer2 += inFrameStride2;                                                                   \
      inPointer3 += inFrameStride3;                                                                   \
      inPointer4 += inFrameStride4;                                                                   \
      inPointer5 += inFrameStride5;                                                                   \
      inPointer6 += inFrameStride6;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_7(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    int inFrameStride2          = module->inputPins[2].outputFrameSize;                               \
    int inVoiceStride2          = module->inputPins[2].outputVoiceStride;                             \
    double *inPointerVoice0Pin2 = module->inputPins[2].outputPointer;                                 \
    int inFrameStride3          = module->inputPins[3].outputFrameSize;                               \
    int inVoiceStride3          = module->inputPins[3].outputVoiceStride;                             \
    double *inPointerVoice0Pin3 = module->inputPins[3].outputPointer;                                 \
    int inFrameStride4          = module->inputPins[4].outputFrameSize;                               \
    int inVoiceStride4          = module->inputPins[4].outputVoiceStride;                             \
    double *inPointerVoice0Pin4 = module->inputPins[4].outputPointer;                                 \
    int inFrameStride5          = module->inputPins[5].outputFrameSize;                               \
    int inVoiceStride5          = module->inputPins[5].outputVoiceStride;                             \
    double *inPointerVoice0Pin5 = module->inputPins[5].outputPointer;                                 \
    int inFrameStride6          = module->inputPins[6].outputFrameSize;                               \
    int inVoiceStride6          = module->inputPins[6].outputVoiceStride;                             \
    double *inPointerVoice0Pin6 = module->inputPins[6].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *inPointer2 = inPointerVoice0Pin2 + voiceIndex * inVoiceStride2;                         \
      double *inPointer3 = inPointerVoice0Pin3 + voiceIndex * inVoiceStride3;                         \
      double *inPointer4 = inPointerVoice0Pin4 + voiceIndex * inVoiceStride4;                         \
      double *inPointer5 = inPointerVoice0Pin5 + voiceIndex * inVoiceStride5;                         \
      double *inPointer6 = inPointerVoice0Pin6 + voiceIndex * inVoiceStride6;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module,                                                                    \
                           inPointer0,                                                                \
                           inPointer1,                                                                \
                           inPointer2,                                                                \
                           inPointer3,                                                                \
                           inPointer4,                                                                \
                           inPointer5,                                                                \
                           inPointer6,                                                                \
                           outPointer,                                                                \
                           voiceIndex);                                                               \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        inPointer2 += inFrameStride2;                                                                 \
        inPointer3 += inFrameStride3;                                                                 \
        inPointer4 += inFrameStride4;                                                                 \
        inPointer5 += inFrameStride5;                                                                 \
        inPointer6 += inFrameStride6;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_7(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_7(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_7(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_7(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_7(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_7(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_7(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_7(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *in3, double *in4,      \
                             double *in5, double *in6, double *in7, double *out, int voiceIndex);     \


// 8 input pins:
#define CREATE_MONO_FRAME_FUNCTION_8(ClassName)                                                       \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    ClassName::process(module,                                                                        \
                       module->inputPins[0].outputPointer,                                            \
                       module->inputPins[1].outputPointer,                                            \
                       module->inputPins[2].outputPointer,                                            \
                       module->inputPins[3].outputPointer,                                            \
                       module->inputPins[4].outputPointer,                                            \
                       module->inputPins[5].outputPointer,                                            \
                       module->inputPins[6].outputPointer,                                            \
                       module->inputPins[7].outputPointer,                                            \
                       module->audioOutputs,                                                          \
                       voiceIndex);                                                                   \
  }                                                                                                   \

#define CREATE_POLY_FRAME_FUNCTION_8(ClassName)                                                       \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                    \
  {                                                                                                   \
    int outVoiceStride = module->outFrameStride * processingStatus.getBufferSize();                   \
    int inVoiceStride0 = module->inputPins[0].outputVoiceStride;                                      \
    int inVoiceStride1 = module->inputPins[1].outputVoiceStride;                                      \
    int inVoiceStride2 = module->inputPins[2].outputVoiceStride;                                      \
    int inVoiceStride3 = module->inputPins[3].outputVoiceStride;                                      \
    int inVoiceStride4 = module->inputPins[4].outputVoiceStride;                                      \
    int inVoiceStride5 = module->inputPins[5].outputVoiceStride;                                      \
    int inVoiceStride6 = module->inputPins[6].outputVoiceStride;                                      \
    int inVoiceStride7 = module->inputPins[7].outputVoiceStride;                                      \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    double *inPointer5 = module->inputPins[5].outputPointer;                                          \
    double *inPointer6 = module->inputPins[6].outputPointer;                                          \
    double *inPointer7 = module->inputPins[7].outputPointer;                                          \
    double *outPointer = module->audioOutputs;                                                        \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex  = voiceAllocator.getPlayingVoiceIndices()[playIndex];                           \
      ClassName::process(module,                                                                      \
                         inPointer0 + voiceIndex * inVoiceStride0,                                    \
                         inPointer1 + voiceIndex * inVoiceStride1,                                    \
                         inPointer2 + voiceIndex * inVoiceStride2,                                    \
                         inPointer3 + voiceIndex * inVoiceStride3,                                    \
                         inPointer4 + voiceIndex * inVoiceStride4,                                    \
                         inPointer5 + voiceIndex * inVoiceStride5,                                    \
                         inPointer6 + voiceIndex * inVoiceStride6,                                    \
                         inPointer7 + voiceIndex * inVoiceStride7,                                    \
                         outPointer + voiceIndex * outVoiceStride,                                    \
                         voiceIndex);                                                                 \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_MONO_BLOCK_FUNCTION_8(ClassName)                                                       \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride = module->outFrameStride;                                                      \
    double *outPointer = module->audioOutputs;                                                        \
    int inFrameStride0 = module->inputPins[0].outputFrameSize;                                        \
    double *inPointer0 = module->inputPins[0].outputPointer;                                          \
    int inFrameStride1 = module->inputPins[1].outputFrameSize;                                        \
    double *inPointer1 = module->inputPins[1].outputPointer;                                          \
    int inFrameStride2 = module->inputPins[2].outputFrameSize;                                        \
    double *inPointer2 = module->inputPins[2].outputPointer;                                          \
    int inFrameStride3 = module->inputPins[3].outputFrameSize;                                        \
    double *inPointer3 = module->inputPins[3].outputPointer;                                          \
    int inFrameStride4 = module->inputPins[4].outputFrameSize;                                        \
    double *inPointer4 = module->inputPins[4].outputPointer;                                          \
    int inFrameStride5 = module->inputPins[5].outputFrameSize;                                        \
    double *inPointer5 = module->inputPins[5].outputPointer;                                          \
    int inFrameStride6 = module->inputPins[6].outputFrameSize;                                        \
    double *inPointer6 = module->inputPins[6].outputPointer;                                          \
    int inFrameStride7 = module->inputPins[7].outputFrameSize;                                        \
    double *inPointer7 = module->inputPins[7].outputPointer;                                          \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                     \
    {                                                                                                 \
      ClassName::process(module,                                                                      \
                         inPointer0,                                                                  \
                         inPointer1,                                                                  \
                         inPointer2,                                                                  \
                         inPointer3,                                                                  \
                         inPointer4,                                                                  \
                         inPointer5,                                                                  \
                         inPointer6,                                                                  \
                         inPointer7,                                                                  \
                         outPointer,                                                                  \
                         voiceIndex);                                                                 \
      inPointer0 += inFrameStride0;                                                                   \
      inPointer1 += inFrameStride1;                                                                   \
      inPointer2 += inFrameStride2;                                                                   \
      inPointer3 += inFrameStride3;                                                                   \
      inPointer4 += inFrameStride4;                                                                   \
      inPointer5 += inFrameStride5;                                                                   \
      inPointer6 += inFrameStride6;                                                                   \
      inPointer7 += inFrameStride7;                                                                   \
      outPointer += outFrameStride;                                                                   \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_POLY_BLOCK_FUNCTION_8(ClassName)                                                       \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                     \
  {                                                                                                   \
    int outFrameStride          = module->outFrameStride;                                             \
    int outVoiceStride          = outFrameStride * processingStatus.getBufferSize();                  \
    int inFrameStride0          = module->inputPins[0].outputFrameSize;                               \
    int inVoiceStride0          = module->inputPins[0].outputVoiceStride;                             \
    double *inPointerVoice0Pin0 = module->inputPins[0].outputPointer;                                 \
    int inFrameStride1          = module->inputPins[1].outputFrameSize;                               \
    int inVoiceStride1          = module->inputPins[1].outputVoiceStride;                             \
    double *inPointerVoice0Pin1 = module->inputPins[1].outputPointer;                                 \
    int inFrameStride2          = module->inputPins[2].outputFrameSize;                               \
    int inVoiceStride2          = module->inputPins[2].outputVoiceStride;                             \
    double *inPointerVoice0Pin2 = module->inputPins[2].outputPointer;                                 \
    int inFrameStride3          = module->inputPins[3].outputFrameSize;                               \
    int inVoiceStride3          = module->inputPins[3].outputVoiceStride;                             \
    double *inPointerVoice0Pin3 = module->inputPins[3].outputPointer;                                 \
    int inFrameStride4          = module->inputPins[4].outputFrameSize;                               \
    int inVoiceStride4          = module->inputPins[4].outputVoiceStride;                             \
    double *inPointerVoice0Pin4 = module->inputPins[4].outputPointer;                                 \
    int inFrameStride5          = module->inputPins[5].outputFrameSize;                               \
    int inVoiceStride5          = module->inputPins[5].outputVoiceStride;                             \
    double *inPointerVoice0Pin5 = module->inputPins[5].outputPointer;                                 \
    int inFrameStride6          = module->inputPins[6].outputFrameSize;                               \
    int inVoiceStride6          = module->inputPins[6].outputVoiceStride;                             \
    double *inPointerVoice0Pin6 = module->inputPins[6].outputPointer;                                 \
    int inFrameStride7          = module->inputPins[7].outputFrameSize;                               \
    int inVoiceStride7          = module->inputPins[7].outputVoiceStride;                             \
    double *inPointerVoice0Pin7 = module->inputPins[7].outputPointer;                                 \
    double *outPointerVoice0    = module->audioOutputs;                                               \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)             \
    {                                                                                                 \
      int voiceIndex     = voiceAllocator.getPlayingVoiceIndices()[playIndex];                        \
      double *inPointer0 = inPointerVoice0Pin0 + voiceIndex * inVoiceStride0;                         \
      double *inPointer1 = inPointerVoice0Pin1 + voiceIndex * inVoiceStride1;                         \
      double *inPointer2 = inPointerVoice0Pin2 + voiceIndex * inVoiceStride2;                         \
      double *inPointer3 = inPointerVoice0Pin3 + voiceIndex * inVoiceStride3;                         \
      double *inPointer4 = inPointerVoice0Pin4 + voiceIndex * inVoiceStride4;                         \
      double *inPointer5 = inPointerVoice0Pin5 + voiceIndex * inVoiceStride5;                         \
      double *inPointer6 = inPointerVoice0Pin6 + voiceIndex * inVoiceStride6;                         \
      double *inPointer7 = inPointerVoice0Pin7 + voiceIndex * inVoiceStride7;                         \
      double *outPointer = outPointerVoice0    + voiceIndex * outVoiceStride;                         \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                   \
      {                                                                                               \
        ClassName::process(module,                                                                    \
                           inPointer0,                                                                \
                           inPointer1,                                                                \
                           inPointer2,                                                                \
                           inPointer3,                                                                \
                           inPointer4,                                                                \
                           inPointer5,                                                                \
                           inPointer6,                                                                \
                           inPointer7,                                                                \
                           outPointer,                                                                \
                           voiceIndex);                                                               \
        inPointer0 += inFrameStride0;                                                                 \
        inPointer1 += inFrameStride1;                                                                 \
        inPointer2 += inFrameStride2;                                                                 \
        inPointer3 += inFrameStride3;                                                                 \
        inPointer4 += inFrameStride4;                                                                 \
        inPointer5 += inFrameStride5;                                                                 \
        inPointer6 += inFrameStride6;                                                                 \
        inPointer7 += inFrameStride7;                                                                 \
        outPointer += outFrameStride;                                                                 \
      }                                                                                               \
    }                                                                                                 \
  }                                                                                                   \

#define CREATE_DERIVED_FUNCTIONS_8(ClassName)                                                         \
  CREATE_MONO_FRAME_FUNCTION_8(ClassName);                                                            \
  CREATE_POLY_FRAME_FUNCTION_8(ClassName);                                                            \
  CREATE_MONO_BLOCK_FUNCTION_8(ClassName);                                                            \
  CREATE_POLY_BLOCK_FUNCTION_8(ClassName);                                                            \

#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_8(ClassName)                                           \
  CREATE_DERIVED_FUNCTIONS_8(ClassName)                                                               \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                                                          \

#define CREATE_COMMON_DECLARATIONS_8(ClassName)                                                       \
  ENFORCE_FACTORY_USAGE(ClassName);                                                                   \
  DECLARE_PROCESSING_FUNCTIONS;                                                                       \
  virtual void initialize();                                                                          \
  static INLINE void process(Module *module, double *in1, double *in2, double *in3, double *in4,      \
                             double *in5, double *in6, double *in7, double *in8, double *out,         \
                             int voiceIndex);                                                         \
















// given the ClassName::process function, this macro creates the corresponding monophonic per-frame 
// processing function:
#define CREATE_MONO_FRAME_FUNCTION_N(ClassName)                                                                           \
  void ClassName::processMonoFrame(Module *module, int voiceIndex)                                                        \
  {                                                                                                                       \
    for(unsigned int pinIndex = 0; pinIndex < ((Module*) module)->getNumInputsInlined(); pinIndex++)                      \
      WorkArea::tmpInFrame[pinIndex] = *(((Module*) module)->inputPins[pinIndex].outputPointer);                          \
    ClassName::process(module, WorkArea::tmpInFrame, ((Module*) module)->audioOutputs, voiceIndex);                       \
  }                                                                                                                       \

// given the ClassName::process function, this macro creates the corresponding polyphonic per-frame processing function:
#define CREATE_POLY_FRAME_FUNCTION_N(ClassName)                                                                           \
  void ClassName::processPolyFrame(Module *module, int voiceIndex)                                                        \
  {                                                                                                                       \
    int outVoiceStride = ((Module*) module)->outFrameStride * processingStatus.getBufferSize();                           \
    double *outPointer = ((Module*) module)->audioOutputs;                                                                \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)                                 \
    {                                                                                                                     \
      int voiceIndex = voiceAllocator.getPlayingVoiceIndices()[playIndex];                                                \
      for(unsigned int pinIndex = 0; pinIndex < ((Module*) module)->getNumInputsInlined(); pinIndex++)                    \
      {                                                                                                                   \
        WorkArea::tmpInFrame[pinIndex] = *(((Module*) module)->inputPins[pinIndex].outputPointer                          \
                                           + voiceIndex * ((Module*) module)->inputPins[pinIndex].outputVoiceStride);     \
      }                                                                                                                   \
      ClassName::process(module, WorkArea::tmpInFrame, outPointer + voiceIndex * outVoiceStride, voiceIndex);             \
    }                                                                                                                     \
  }                                                                                                                       \

// given the ClassName::process function, this macro creates the corresponding monophonic per-blaoc processing function:
#define CREATE_MONO_BLOCK_FUNCTION_N(ClassName)                                                                           \
  void ClassName::processMonoBlock(Module *module, int voiceIndex, int blockSize)                                         \
  {                                                                                                                       \
    int outFrameStride = ((Module*) module)->outFrameStride;                                                              \
    double *outPointer = ((Module*) module)->audioOutputs;                                                                \
    for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                                         \
    {                                                                                                                     \
      for(unsigned int pinIndex = 0; pinIndex < ((Module *) module)->getNumInputsInlined(); pinIndex++)                   \
      {                                                                                                                   \
        WorkArea::tmpInFrame[pinIndex] = * (((Module *) module)->inputPins[pinIndex].outputPointer                        \
                                            + frameIndex * ((Module *) module)->inputPins[pinIndex].outputFrameSize);     \
      }                                                                                                                   \
      ClassName::process(module, WorkArea::tmpInFrame, outPointer + outFrameStride * frameIndex, voiceIndex);             \
    }                                                                                                                     \
  }                                                                                                                       \

// given the ClassName::process function, this macro creates the corresponding polyphonic per-block processing function:
#define CREATE_POLY_BLOCK_FUNCTION_N(ClassName)                                                                               \
  void ClassName::processPolyBlock(Module *module, int voiceIndex, int blockSize)                                             \
  {                                                                                                                           \
    unsigned int pinIndex;                                                                                                    \
    int outFrameStride = ((Module*) module)->outFrameStride;                                                                  \
    int outVoiceStride = outFrameStride * processingStatus.getBufferSize();                                                   \
    double *outVoiceFramePointer;                                                                                             \
    for(int playIndex = 0; playIndex < voiceAllocator.getNumPlayingVoices(); playIndex++)                                     \
    {                                                                                                                         \
      int voiceIndex       = voiceAllocator.getPlayingVoiceIndices()[playIndex];                                              \
      outVoiceFramePointer = ((Module*) module)->audioOutputs + voiceIndex * outVoiceStride;                                  \
      for(pinIndex = 0; pinIndex < ((Module *) module)->getNumInputsInlined(); pinIndex++)                                    \
      {                                                                                                                       \
        WorkArea::inVoiceFramePointer[pinIndex] = ((Module *) module)->inputPins[pinIndex].outputPointer                      \
                                                  + voiceIndex * ((Module *) module)->inputPins[pinIndex].outputVoiceStride   \
                                                  - ((Module *) module)->inputPins[pinIndex].outputFrameSize;                 \
          /* points 1 position before the actual start, to allow increment before dereferencing  */                           \
      }                                                                                                                       \
      for(int frameIndex = 0; frameIndex < blockSize; frameIndex++)                                                           \
      {                                                                                                                       \
        for(pinIndex = 0; pinIndex < ((Module *) module)->getNumInputsInlined(); pinIndex++)                                  \
        {                                                                                                                     \
          WorkArea::inVoiceFramePointer[pinIndex] += ((Module *) module)->inputPins[pinIndex].outputFrameSize;                \
          WorkArea::tmpInFrame[pinIndex]           = *(WorkArea::inVoiceFramePointer[pinIndex]);                              \
        }                                                                                                                     \
        ClassName::process(module, WorkArea::tmpInFrame, outVoiceFramePointer, voiceIndex);                                   \
        outVoiceFramePointer   += outFrameStride;                                                                             \
      }                                                                                                                       \
    }                                                                                                                         \
  }                                                                                                                           \

// given a monophonic per-frame processing function, this macro creates the corresponding other 3 processing functions:
#define CREATE_DERIVED_FUNCTIONS_N(ClassName)    \
  CREATE_MONO_FRAME_FUNCTION_N(ClassName);       \
  CREATE_POLY_FRAME_FUNCTION_N(ClassName);       \
  CREATE_MONO_BLOCK_FUNCTION_N(ClassName);       \
  CREATE_POLY_BLOCK_FUNCTION_N(ClassName);       \

// combines CREATE_DERIVED_FUNCTIONS and CREATE_ASSIGN_FUNCTION_POINTERS into a single macro:
#define CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_N(ClassName) \
  CREATE_DERIVED_FUNCTIONS_N(ClassName)                     \
  CREATE_ASSIGN_FUNCTION_POINTERS(ClassName)                \

#define CREATE_COMMON_DECLARATIONS_N(ClassName)                                            \
  ENFORCE_FACTORY_USAGE(ClassName);                                                        \
  DECLARE_PROCESSING_FUNCTIONS;                                                            \
  virtual void initialize();                                                               \
  static INLINE void process(Module *module, double  *ins, double *outs, int voiceIndex);  \

#endif 
