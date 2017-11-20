////////////////////////////////////////////////////////////////////////////////////
//     /*! \file aufTAKT_If.h: \brief interface of the CaufTAKT_If class. */
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
//  $RCSfile: aufTAKT_If.h,v $
//  $Author: flea $
//  $Date: 2009-09-08 15:16:45 $
//
//
//
////////////////////////////////////////////////////////////////////////////////////

#if !defined(__libaufTAKT_If_HEADER_INCLUDED__)
#define __libaufTAKT_If_HEADER_INCLUDED__

#ifndef _BEAT_STRUCTS
#define _BEAT_STRUCTS
struct stBeatInfo
{
    int        lPos;               //!< position of the beat mark in samples
    float      fBPM;               //!< current estimate of the tempo at the beat mark
    float      fProbability;       //!< probability of the current beat mark (experimental)
};

struct stBeatInfoEntry
{
    int        lPos;               //!< position of the onset mark in samples
    float      fTransitionEnergy;  //!< transition energy for the onset marks (only used for onset marks)
    float      fBPM;               //!< current estimate of the tempo at the beat mark (only used for beat marks)
    float      fProbability;       //!< probability of the current beat mark (only used for beat marks)
};

struct stBeatGrid
{
    int     iLastIdx;
    
    float   fPos,
            fLastInter,
            fNextBeat,                
            fNoteLen,
            fInvNoteLen,
            fBPM,
            fMeanProbability,
            fProbability;
};
#endif

/*!
    This class provides the interface for the [aufTAKT] beat tracking SDK
*/
class CaufTAKT_If
{
public:

    enum eaufTAKTVersion_t
    { 
        kaufTAKT2,
        kaufTAKT3 
    };

    /*!
    * only V3: allows the user to enforce (the beat marks to) a constant tempo
    */
    enum eTempoModes_t
    {
        kDynamicTempo,      //!< default mode: [aufTAKT] allows the tempo to vary slightly over time
        kStraightTempo,     //!< new V3 mode: [aufTAKT] forces the resulting tempo to be constant over time

        kNumOfTempoModes
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
     * creates an instance of the beat tracking class
     *
     * @param pCaufTAKT_If : pointer to the instance to be created
     * @param iSampleRate : sample rate of input signal
     * @param iNumOfChannels : number of audio channels
     * @param eVersion :  selects either aufTAKT version 2 or 3 (kaufTAKT2 or kaufTAKT3)
     * @param fDownBeatBlockSizeInSec: defines the blocksize for the downbeat detection and initial bpm estimation in aufTAKT3 (default = 30 sec), has no influence on aufTAKT2 - use at your own risk!
     *
     * @return static int  : returns some error code otherwise NULL 
     */
    static int      CreateInstance (CaufTAKT_If*& pCaufTAKT_If, int iSampleRate, int iNumOfChannels, eaufTAKTVersion_t eVersion = kaufTAKT3, float fDownBeatBlockSizeInSec = 30.0);
    
    /*!
     * destroys an instance of the synthesis class
     *
     * @param pCaufTAKT_If : pointer to the instance to be destroyed
     *
     * @return static int  : returns some error code otherwise NULL
     */
    static int      DestroyInstance (CaufTAKT_If*&  pCaufTAKT_If);
    
    
    /*!
     * resets the internal state of the synthesis
     *
     * @return int  : returns error code otherwise NULL
     */
    virtual int     Reset () = 0;

    /*!
     * does the onset extraction from the audio signal in blocks
     *
     * @param **ppfInputBuffer : pointer to array of input buffers (ppfInputBuffer[0] points to buffer of first channel with length iNumFrames, ppfInputBuffer[1] points to buffer of second channel with the same length
     * @param iNumOfFrames : number of frames in channel buffers
     *
     * @return virtual int  : returns some error code otherwise NULL
     */
    virtual int     PreAnalysis (float **ppfInputBuffer, int iNumOfFrames) = 0; 

    /*!
     * tells the library that the pre-processing has finished and calculated initial BPM estimate
     *
     * @param bDoInitialBPMEstimate : if the initial BPM estimate shall not be calculated, set this value to 0
     *
     * @return virtual int  : returns some error code otherwise NULL
     */
    virtual int     FinishPreAnalysis (bool bDoInitialBPMEstimate = true) = 0; 

    /*!
     * allows user to check if the initial BPM estimate can be requested before the call of FinishPreAnalysis (optional)
     *
     * @return virtual bool   : true if ::CalculateInitialBPMEstimate is ready to be called
     */
    virtual bool    IsBPMEstimateReady ()   = 0;          

    /*!
     * calculates the initial BPM estimate even if the user has not called FinishPreAnalysis (optional), a subsequent call of FinishPreAnalysis may then be called with parameter bDoInitialBPMEstimate == false
     *
     * @return virtual float  : initial BPM Estimate
     */
    virtual float   CalculateInitialBPMEstimate ()  = 0;          

    /*!
     * after the pre-processing, the results can be saved to a buffer/file to avoid the preprocessing step the next time
     *
     * @param *&pstOnsetInfo : returns pointer to structure ::BeatInfoEntry with onset data
     * @param *piNumOfBeatInfoEntries : number of entries in pstOnsetInfo (complete size of memory in bytes can be calculated with (*piNumOfBeatInfoEntries * sizeof(BeatInfoEntry)))
     * @param *&pstBeatTrackState : returns pointer to internal state information
     * @param *piNumOfGridPoints : number of state data to store in grid points (complete size of memory in bytes can be calculated with (*piNumOfGridPoints * sizeof(stBeatGrid)))
     * @param *pfBeatTrackLock : another internal state variable
     * @param *piDownBeatLocation: most probable location of the downbeat (only valid for aufTAKT3)
     * @param *piBarLength: most probable length of a bar (only valid for aufTAKT3)
     *
     * @return virtual int  : returns some error code otherwise NULL
     */
    virtual int     GetPreAnalysisResult (stBeatInfoEntry *&pstOnsetInfo, int *piNumOfBeatInfoEntries, stBeatGrid *&pstBeatTrackState, int *piNumOfGridPoints, float *pfBeatTrackLock, int *piDownBeatLocation = 0, int *piBarLength = 0) = 0; 

    /*!
     * to avoid the preprocessing step, the once extracted data can be imported before the processing (see function GetPreAnalysisResult)
     *
     * @param *pstOnsetInfo : pointer to structure ::BeatInfoEntry with onset data
     * @param iNumOfBeatInfoEntries : number of entries in pstOnsetInfo (complete size of memory in bytes can be calculated with (iNumOfBeatInfoEntries * sizeof(BeatInfoEntry)))
     * @param *pstBeatTrackState : pointer to internal state information
     * @param iNumOfGridPoints : number of state data in grid points (complete size of memory in bytes can be calculated with (iNumOfGridPoints * sizeof(stBeatGrid)))
     * @param fBeatTrackLock : another internal state variable
     * @param iDownBeatLocation: most probable location of the downbeat (only valid for aufTAKT3)
     * @param iBarLength: most probable length of a bar (only valid for aufTAKT3)
     *
     * @return virtual int  : returns some error code otherwise NULL
     */
    virtual int     SetPreAnalysisResult (stBeatInfoEntry *pstOnsetInfo, int iNumOfBeatInfoEntries, stBeatGrid *pstBeatTrackState, int iNumOfGridPoints, float fBeatTrackLock, int iDownBeatLocation = 0, int iBarLength = 0) = 0; 

    /*!
     * returns the initial BPM estimate after pre-processing
     *
     * @return virtual float  : initial BPM estimate
     */
    virtual float   GetInitialBPMEstimate () = 0;          

    /*!
     * set the initial BPM estimate before processing
     *
     * @param fInitialBPM : BPM value
     *
     * @return virtual int  : returns some error code otherwise NULL
     */
    virtual int     SetInitialBPMEstimate (float fInitialBPM) = 0;          

    /*!
     * does the actual processing for beat mark extraction
     *
     * @param fAdaptSpeed: optional adaption parameter (range 1..50, 1 means fast adaptation (loose tempo), 50 means very slow adaption (straight tempo), default 0 means use internal value)
     * @param iDownBeatPos: optional downbeat position if internal guess is wrong, if iDownBeatPos == 0 internal guess is taken
     * @param bStartFromHere : optionally start beat tracking from downbeat position in both directions
     * @param bUseBPMAdaption : enables adaption of some internal params to current bpm estimate
     * 
     * @return int  : returns some error code otherwise NULL
     */
    virtual int     Process (float fAdaptSpeed  = 0, int iDownBeatPos = 0, bool bStartFromHere = false, bool bUseBPMAdaption = false) = 0;          

    /*!
     * returns the number of extracted beats
     *
     * @param eTempoMode : assume dynamic (default) or constant tempo
     *
     * @return virtual int  : number of beats
     */
    virtual int     GetNumBeatMarks (eTempoModes_t eTempoMode = kDynamicTempo) = 0;

    /*!
     * calculates the most probable overall tempo guess (only meaningful for files with relatively constant tempo). Must not be called before process
     *
     * @param eTempoMode : assume dynamic (default) or constant tempo
     *
     * @return virtual int  : estimated overall tempo
     */
    virtual float   GetOverallTempo (eTempoModes_t eTempoMode = kDynamicTempo) = 0;

    /*!
    * returns how many beats are within one bar (will always be 8 for auftakt2)
    *
    * @param eTempoMode : assume dynamic (default) or constant tempo
    * @return virtual int  : estimated number of beats per bar
    */
    virtual int     GetTimeSignatureIdx(eTempoModes_t eTempoMode = kDynamicTempo) = 0;

    /*!
     * estimates the first beat which is a downbeat
     *
     * @param iBeatDivisor : the number of beats a measure is supposed to be divided in (note: set to 0 to use internal estimate
     * @param eTempoMode : assume dynamic (default) or constant tempo
     *
     * @return virtual int  : index of first downbeat estimate
     */
    virtual int     GetFirstDownBeatEstimate (int iBeatDivisor = 0, eTempoModes_t eTempoMode = kDynamicTempo) = 0;

    /*!
     * returns information about one beat
     *
     * @param *pBeat : structure for beat information
     * @param iIdx : index of beat
     * @param eTempoMode : assume dynamic (default) or constant tempo
    *
     * @return int  : returns some error code otherwise NULL
     */
    virtual int     GetBeatMark (stBeatInfo *pBeat, int iIdx, eTempoModes_t eTempoMode = kDynamicTempo) = 0;

	 /*!
     * Enable special mode optimized for beat based music. This mode favours the detection of 4/4 resp 8/8 beats. This must en/disabled before preprocessing. 
	 * In case one decides to change this settings the whole preprocessing process must be repeated! By default this is disabled.
     *
     * @param bEnable : set this to true in order to enable the beat optimized preprocessing.
     *
     */
    virtual void    SetBeatMusicOpt(bool bEnable) = 0; 

};


#endif // __libaufTAKT_If_HEADER_INCLUDED__
