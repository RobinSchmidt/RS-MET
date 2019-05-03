// TuningMap.h: Interface of the class CTuningMap.
//
// (C)opyright in 2003 by Mark Henning, Germany
//
// Read TuningMap.cpp for more informations about this class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TUNINGMAP_H__15693DC4_FB37_11D6_A827_F4C607C10000__INCLUDED_)
#define AFX_TUNINGMAP_H__15693DC4_FB37_11D6_A827_F4C607C10000__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <fstream>
#include <math.h>
using namespace std;

class CTuningMap  
{
public:
	CTuningMap();
	virtual ~CTuningMap();

	void	Reset();

	bool	WriteToFile(const char * szFilepath, bool bSaveBaseFreq = false);
	bool	ReadFromFile(const char * szFilepath);

	double	GetBaseFreq() const; // BaseFreq in Hz
	double	GetNoteFreq(int nNoteIndex);       // Absolute tune in Hz
	double	GetRelativeTune(int nNoteIndex);   // Relative tune in cents

	bool	SetBaseFreq(double dblBaseFreq); // BaseFreq in Hz
	bool	SetRelativeTune(int nNoteIndex, double dblTune); // Relative tune in cents

	// Call this directly, when one of the above functions failed:
	const char *	GetLastError() const {return m_szErrorString;}

private:
	double	m_dblTunes[128];	// Unit: Cents
	double	m_dblBaseFreq;		// Unit: Hz

	char	m_szErrorString[400];

	enum eSection {
		SEC_None = 0,
		SEC_Unknown = 1,
		SEC_Tuning = 2,
		SEC_ExactTuning = 3
	};
};

#endif // !defined(AFX_TUNINGMAP_H__15693DC4_FB37_11D6_A827_F4C607C10000__INCLUDED_)
