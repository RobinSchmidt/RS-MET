// TuningMap.cpp: Implementation of the class CTuningMap.
//
// (C)opyright in 2003 by Mark Henning, Germany
//
// Contact email: info@anamark.de or mh@homolog.de
//
// You may use this code for free. If you find an error or make some
// interesting changes, please let me know.
//
// This class deals with the reading/writing of AnaMark / VAZ 1.5 Plus-
// compatible tuning files. Be carefull with changes of the functions
// WriteToFile/ReadFromFile because this may lead to incompatibilities!
//
// I think, the source-code is rather self-explaining.
//
// The specifications of the AnaMark / VAZ 1.5 Plus-compatible
// tuning file format can be found at http://www.anamark.de
//
// IMPORTANT:
// Please note, that
//     VAZ 1.5 Plus tuning files
// is NOT the same as
//     AnaMark / VAZ 1.5 Plus-compatible tuning files
//
// This code deals with the latter one, not with the first one.
//
// Have fun!
//
//////////////////////////////////////////////////////////////////////



#include "TuningMap.h"





//////////////////////////////////////////////////////////////////////
// Konstruktion/Destruktion
//////////////////////////////////////////////////////////////////////





CTuningMap::CTuningMap()
{
	// Provide a standard tuning
	Reset();
}



CTuningMap::~CTuningMap()
{

}





//////////////////////////////////////////////////////////////////////
// Public Functions
//////////////////////////////////////////////////////////////////////





void CTuningMap::Reset()
{
	// This function _must_ never produce an error, so we don't need
	// to return a bool value...
	strcpy(m_szErrorString, "");

	m_dblBaseFreq = 8.1757989156437073336; // Means A = 440Hz

	for ( int i = 0 ; i < 128 ; ++i )
		m_dblTunes[i] = 100 * i;
}



bool CTuningMap::WriteToFile(const char * szFilepath, bool bSaveBaseFreq /*= false*/)
{
	strcpy(m_szErrorString, "");

  ofstream	ofs(szFilepath, ios::trunc | ios::out);
	int			i;

	ofs << ";" << endl;
	ofs << "; AnaMark / VAZ 1.5 Plus-compatible tuning map file" << endl;
	ofs << ";" << endl;
	ofs << ";" << endl;
	ofs << "; 1. VAZ-section with quantized tunings" << endl;
	ofs << ";" << endl;

	ofs << "[Tuning]" << endl;
	for ( i = 0 ; i < 128 ; ++i )
		ofs << "note " << i << "=" << (long)(m_dblTunes[i]) << endl;

	ofs << ";" << endl;
	ofs << "; 2. AnaMark-specific section with exact tunings" << endl;
	ofs << ";" << endl;
	ofs << "[Exact Tuning]" << endl;
	ofs.precision(10);
	if ( bSaveBaseFreq )
		ofs << "basefreq = " << m_dblBaseFreq << endl;
	for ( i = 0 ; i < 128 ; ++i )
		ofs << "note " << i << "=" << m_dblTunes[i] << endl;

	ofs.close();

	// Always returns true, because currently there's no error processing
	return true;
}



bool CTuningMap::ReadFromFile(const char * szFilepath)
{
	char		szLine[512]; // Max. line length is 255 characters, but to be on the save side...

	eSection	secCurr = SEC_None;

	long		lTunes[128];		// For temporary use, unit: Cents

	bool		bTuningFound = false;
	bool		bExactTuningFound = false;

	long		lET_LastNoteFound = -1;

	long		lLineCount = 0;

	// Initialize data
	// Important, because notes not listed in the tuning file
	// should always have standard tuning.
	Reset();
	for ( int i = 0 ; i < 128 ; ++i )
		lTunes[i] = (long)m_dblTunes[i];

	// Now open the file
  ifstream	ifs(szFilepath, ios::in | ios::binary);
	//ifstream	ifs(szFilepath, ios::in | ios::nocreate | ios::binary);

	if ( !ifs )
	{
		sprintf(m_szErrorString, "Error opening the file  '%s'", szFilepath);
		return false;
	}

	while ( ifs )
	{
		// Increase Line counter to make it easier detecting errors, if
		// a more detailed output is wanted.
		++lLineCount;

		// Read line, until '\n', '\r' or '\0' is reached
		// Thus it is able to read WIN/DOS-style as well as UNIX-style files
		// By the way: Skip empty lines or multiple line-end-characters
		// Is not case sensitive, so all chars are converted to lower ones
		int	nCurrPos = 0;
		do
		{
			while ( (ifs) && (nCurrPos < 510) )
			{
				char	ch = '\0';
				ifs.read(&ch, 1);
				if ( (ch == '\0') || (ch == '\r') || (ch == '\n') )
					break;
				szLine[nCurrPos++] = (char) tolower(ch);
			}
		} while ( (ifs) && (nCurrPos == 0) );
		if ( nCurrPos >= 510 )
		{
			sprintf(m_szErrorString, "Line too long (line %d)", (int) lLineCount);
			return false; // Line too long
		}
		szLine[nCurrPos] = '\0';

		// Skip leading spaces/tabs
		const char	* szCurr = szLine;
		while ( (*szCurr != '\0') && ((*szCurr == ' ') || (*szCurr == '\t')) )
			++szCurr;

		// Skip empty lines
		if ( szCurr[0] == '\0' )
			continue;

		// Skip comment lines
		if ( szCurr[0] == ';' )
			continue;

		// Skip trailing spaces/tabs
		char	* szLast = (char *) &szCurr[strlen(szCurr)-1];
		while ( (szLast > szCurr) && ((*szLast == ' ') || (*szLast == '\t')) )
			--szLast;
		*(szLast+1) = '\0';

		// Check for new section
		if ( szCurr[0] == '[' )
		{
			if ( szCurr[strlen(szCurr)-1] != ']' )
			{
				sprintf(m_szErrorString, "Syntax error: Section-tag must be the only string in the line! (line %d)", (int) lLineCount);
				return false; // error in section-tag! Must be the only one string in the line!
			}
			// Known section found?
			secCurr = SEC_Unknown;
			if ( strcmp(&szCurr[1], "tuning]") == 0 )
			{
				secCurr = SEC_Tuning;
				bTuningFound = true;
			}
			if ( strcmp(&szCurr[1], "exact tuning]") == 0 )
			{
				secCurr = SEC_ExactTuning;
				bExactTuningFound = true;
			}
			// Now process next line
			continue;
		}

		// Skip all lines which are in none or in an unknown section
		if ( (secCurr == SEC_None) || (secCurr == SEC_Unknown) )
			continue;

		// Separate parameter name and value
		const char	* szParam = szCurr;
		const char	* szValue = strchr(szCurr, '=');
		if ( szValue == NULL )
		{
			sprintf(m_szErrorString, "Syntax error: '=' missing! (line %d)", (int) lLineCount);
			return false; // definitely an error: '=' missing!
		}
		++szValue; // Set the pointer to the first char behind the '='
		// Now skip trailing spaces/tabs of parameter name:
		szLast = (char *) &szValue[-2];
		while ( (szLast > szParam) && ((*szLast == ' ') || (*szLast == '\t')) )
			--szLast;
		*(szLast+1) = '\0';
		// Now skip leading spaces/tabs of value:
		while ( (*szValue != '\0') && ((*szValue == ' ') || (*szValue == '\t')) )
			++szValue;

		// Now process the different sections:
		switch ( secCurr )
		{
		case SEC_Tuning:
			{
				// Check for note-tag
				if ( memcmp(szParam, "note", 4) != 0 )
					continue; // note-tag not found, ignore line in case that it's
							  // an option of a later version
							  // If you want, you can return here false and process it as an error

				// Get MIDI-Note number
				long	lNoteIndex = atol(&szParam[4]);
				// Check for correct range [0;127] and ignore it, if it's out
				// of range.
				if ( (lNoteIndex < 0) || (lNoteIndex > 127) )
					continue;
				lTunes[lNoteIndex] = atol(szValue);
			}
			break;

		case SEC_ExactTuning:
			{
				// Check for note-tag
				if ( memcmp(szParam, "note", 4) == 0 )
				{
					// note-tag found
					// Get MIDI-Note number
					long	lNoteIndex = atol(&szParam[4]);
					// Check for correct range [0;127] and ignore it, if it's out
					// of range.
					if ( (lNoteIndex < 0) || (lNoteIndex > 127) )
						continue;
					m_dblTunes[lNoteIndex] = atof(szValue);
					if ( lET_LastNoteFound < lNoteIndex )
						lET_LastNoteFound = lNoteIndex;
					// O.K. -> Process next line
					continue;
				}
				// Check for basefreq parameter
				if ( strcmp(szParam, "basefreq") == 0 )
				{
					// basefreq found
					m_dblBaseFreq = atof(szValue);
					// O.K. -> Process next line
					continue;
				}
				// No known parameter found, so skip it
				continue;
			}
			break;

		default: // This part of the code should never be reached!
			assert(false);
			sprintf(m_szErrorString, "Compilation error! Section not coded! (line %d)", (int) lLineCount);
			return false;
		}
	}

	if ( (!bTuningFound) && (!bExactTuningFound) )
	{
		sprintf(m_szErrorString, "No tuning data found!");
		return false; // No tuning data found at all should be worth an error...
	}

	if ( !bExactTuningFound )
	{
		// There are no exact tuning values, so map the quantized
		// values to the exact ones:
		for ( int i = 0 ; i < 128 ; ++i )
			m_dblTunes[i] = (double) lTunes[i];
	}
	else
	{
		// [Exact Tuning] section found, so ignore the values found
		// in the [Tuning] section and do the "auto expand":
		if ( (lET_LastNoteFound >= 0) && (lET_LastNoteFound < 127) )
		{
			// Now loop the given data (auto expand):
			int	H = (int)lET_LastNoteFound;	// Highest MIDI note number
			double	P = m_dblTunes[H];		// Period length
			for ( int i = H ; i < 128 ; ++i )
				m_dblTunes[i] = m_dblTunes[i-H] + P;
		}
	}

	return true; // Everything nice!
}



double CTuningMap::GetBaseFreq() const
{
	return m_dblBaseFreq;
}



double CTuningMap::GetNoteFreq(int nNoteIndex)
{
	return GetBaseFreq() * pow(2, GetRelativeTune(nNoteIndex) / 1200.);
}



double CTuningMap::GetRelativeTune(int nNoteIndex)
{
	strcpy(m_szErrorString, "");
	// First make sure, that the note index is in the valid range
	// If not, return a "standard value"
	if ( (nNoteIndex >= 0) && (nNoteIndex <= 127) )
		return m_dblTunes[nNoteIndex];
	else
		return 100 * nNoteIndex;
}



bool CTuningMap::SetBaseFreq(double dblBaseFreq)
{
	strcpy(m_szErrorString, "");
	// First make sure, that the base frequency is in the valid range
	// If not, return false;
	if ( dblBaseFreq > 0 )
	{
		m_dblBaseFreq = dblBaseFreq;
		return true;
	}
	else
	{
		sprintf(m_szErrorString, "Base frequency out of range: %f", dblBaseFreq);
		return false;
	}
}



bool CTuningMap::SetRelativeTune(int nNoteIndex, double dblTune)
{
	strcpy(m_szErrorString, "");
	// First make sure, that the note index is in the valid range
	// If not, return false;
	if ( (nNoteIndex >= 0) && (nNoteIndex <= 127) )
	{
		m_dblTunes[nNoteIndex] = dblTune;
		return true;
	}
	else
	{
		sprintf(m_szErrorString, "Note index out of range: %d", nNoteIndex);
		return false;
	}
}
