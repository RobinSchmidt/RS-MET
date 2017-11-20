////////////////////////////////////////////////////////////////////////////////////
// /*! \file elastiqueMCAPI.h */ 
//
//  Copyright (c) 2000-2011  
//      zplane.development GmbH & Co. KG
//
//  CONFIDENTIALITY:
//
//      This file is the property of zplane.development.
//      It contains information that is regarded as privilege
//      and confidential by zplane.development.
//      It may not be publicly disclosed or communicated in any way without 
//      prior written authorization by zplane.development.
//      It cannot be copied, used, or modified without obtaining
//      an authorization from zplane.development.
//      If such an authorization is provided, any modified version or
//      copy of the software has to contain this header.
//
//  WARRANTIES: 
//      This software is provided as << is >>, zplane.development 
//      makes no warranty express or implied with respect to this software, 
//      its fitness for a particular purpose or its merchantability. 
//      In no case, shall zplane.development be liable for any 
//      incidental or consequential damages, including but not limited 
//      to lost profits.
//
//      zplane.development shall be under no obligation or liability in respect of 
//      any infringement of statutory monopoly or intellectual property 
//      rights of third parties by the use of elements of such software 
//      and User shall in any case be entirely responsible for the use 
//      to which he puts such elements.	
//
////////////////////////////////////////////////////////////////////////////////////
//  CVS INFORMATION
//
//  $RCSfile: elastiqueMCAPI.h,v $
//  $Author: flea $
//  $Date: 2008-07-21 14:58:43 $
//
//
//
////////////////////////////////////////////////////////////////////////////////////


#if !defined(__libELASTIQUEMCAPI_HEADER_INCLUDED__)
#define __libELASTIQUEMCAPI_HEADER_INCLUDED__



/*!
CLASS
    

    This class provides the interface for zplane's élastique class. 

USAGE
    
*/
class CElastiqueMCIf
{
public:
    
    virtual ~CElastiqueMCIf() {};

    enum _elastique_proc_mode
    {
        _TRANSIENT_MODE,    //!< mode optimized for audio rich of transients
        _TONAL_MODE,        //!< DEPRECATED: mode optimized for audio rich of tonal components
        _E_PRO_MODE
    };

    enum _ElastiqueStereoInputMode_
    {
        _PLAIN_STEREO   = 0,    //!< normal LR stereo mode
        _MS_STEREO              //!< MS stereo mode M must be in channel 0 and S in channel 1
    };
    

    enum Version_t
    {
        kMajor,
        kMinor,
        kPatch,
        kBuild,

        kNumVersionInts
    };

    static const int  GetVersion (const Version_t eVersionIdx);
    static const char* GetBuildDate ();

    /*!
     * creates an instance of zplane's élastique class
     *
     * @param cCElastique : returns a pointer to the class instance
     * @param iOutputBufferSize : desired number of frames at the output 
     * @param iNumOfChannels : number of channels (1..2)
     * @param fSampleRate : input samplerate
     * @param eMode : either transient or tonal optimized mode
     *
     * @return static int  : returns some error code otherwise NULL 
     */
    static int      CreateInstance (CElastiqueMCIf*& cCElastique, int iOutputBufferSize, int iNumOfChannels, float fSampleRate, _elastique_proc_mode eMode);
    
    /*!
     * destroys an instance of the zplane's élastique class
     *
     * @param cCElastique : pointer to the instance to be destroyed
     *
     * @return static int  : returns some error code otherwise NULL
     */
    static int      DestroyInstance(CElastiqueMCIf*  cCElastique);
    

    /*!
     * does the actual processing if the number of frames provided is as retrieved by CElastiqueMCIf::GetFramesNeeded()
     * this function always returns the number of frames as specified when calling CElastiqueMCIf::CreateInstance(..)
     *
     * @param ppInSampleData : double pointer to the input buffer of samples [channels][samples]
     * @param iNumOfInFrames :  the number of input frames
     * @param ppOutSampleData : double pointer to the output buffer of samples [channels][samples]
     *
     * @return virtual int  : returns the error code, otherwise 0
     */
    virtual int     ProcessData(float** ppInSampleData, int iNumOfInFrames, float** ppOutSampleData) = 0;


    /*!
     * returns the number of frames needed for the next processing step according to the required buffersize, this function should be always called directly before CElastiqueMCIf::ProcessData(..)
     * this function changes the internal output buffer size as specified when calling CElastiqueMCIf::CreateInstance(..)
     * use this function if want to change the output size.
     *
     * @param iOutBufferSize : required output buffer size
     *
     * @return virtual int  : returns the number of frames required
     */
    virtual int     GetFramesNeeded(int iOutBufferSize) = 0;


    /*!
     * returns the number of frames needed for the next processing step, this function should be always called directly before CElastiqueMCIf::ProcessData(..)
     *
     * @param none
     *
     * @return virtual int  : returns the number of frames required
     */
    virtual int     GetFramesNeeded() = 0;

    /*!
    * returns the maximum number of frames needed for a given minimum factor
    *
    * @param float fMinStretchFactor : the minimum stretch factor to be used
    * @param float fMinPitchFactor   : the minimum pitch factor to be used
    *
    * @return virtual int : returns number of frames
    */
    virtual int     GetMaxFramesNeeded(float fMinStretchFactor = 0.1, float fMinPitchFactor = 1.0) = 0;

    /*!
    * sets the internal stretch & pitch factor. A value between 0.1 and 10.0 for stretching is valid. The stretch factor is quantized. 
    * The product of the stretch factor and the pitch factor must be between 0.1 and 10.0.
    *
    * @param fStretchFactor : stretch factor 0.1 - 10.0 (1.0 is default and does nothing) 
    * @param fPitchFactor: pitch factor 0.25 - 4.0 (1.0 is default and does nothing) 
    * @param bUsePitchSync: synchronizes timestretch and pitchshifting (default is set to _FALSE in order to preserve old behavior, see documentation for V2.1), this is ignored in SOLO modes
    *
    * @return int  : returns the error code, otherwise 0
    */
    virtual int     SetStretchQPitchFactor(float& fStretchFactor, float fPitchFactor, bool bUsePitchSync= false) = 0;

    /*!
    * sets the internal stretch & pitch factor. A value between 0.1 and 10.0 for stretching is valid. The pitch factor is quantized
    * The product of the stretch factor and the pitch factor must be between 0.1 and 10.0.
    *
    * @param fStretchFactor : stretch factor 0.1 - 10.0 (1.0 is default and does nothing) 
    * @param fPitchFactor: pitch factor 0.25 - 4.0 (1.0 is default and does nothing) 
    * @param bUsePitchSync: synchronizes timestretch and pitchshifting (default is set to _FALSE in order to preserve old behavior, see documentation for V2.1), this is ignored in SOLO modes
    *
    * @return int  : returns the error code, otherwise 0
    */
    virtual int     SetStretchPitchQFactor(float fStretchFactor, float& fPitchFactor, bool bUsePitchSync= false) = 0;

    /*!
     * resets the internal state of the élastique algorithm
     *
     * @param none
     *
     * @return int  : returns some error code otherwise NULL
     */
    virtual void    Reset() = 0;
    

    /*!
     * gets the last frames in the internal buffer, ProcessData must have returned -1 before.
     *
     * @param ppfOutSampleData: double pointer to the output buffer of samples [channels][samples]
     *
     * @return int  : returns some error code otherwise NULL
     */
    virtual int     FlushBuffer(float** ppfOutSampleData) = 0;

    
    /*!
     * returns the number of frames already buffered in the output buffer. 
     * The number of frames is divided by the stretch factor, so that one is able to calculate the current position within the source file.
     *
     * @param none
     *
     * @return virtual int  : returns the number of frames buffered 
     */
    virtual int     GetFramesBuffered() = 0;

  
    /*!
     * sets the stereo input mode
     *
     * @param eStereoInputMode : sets the mode according to _ElastiqueStereoInputMode_
     * 
     * @return static int  : returns some error code otherwise NULL 
     */
    virtual int SetStereoInputMode(_ElastiqueStereoInputMode_ eStereoInputMode) = 0;

	/*!
     * sets a cutoff frequency in order to save processing time
     *
     * @param fFreq : cutoff freq in Hz
     *
     * @return virtual int  : returns some error code otherwise NULL 
     */
	virtual int SetCutOffFreq(float fFreq) = 0;

	/*!
     * gets the current input time position
     *
     * @return virtual int  : returns the time position
     */
    virtual int GetTimePos() = 0;


};


#endif // #if !defined(__libELASTIQUEMCAPI_HEADER_INCLUDED__)



