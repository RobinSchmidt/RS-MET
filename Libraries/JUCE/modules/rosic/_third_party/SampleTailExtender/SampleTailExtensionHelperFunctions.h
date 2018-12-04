//=================================================================================
//  SampleTailExtension
//
//  Written by Adam Stark
//  Copyright © 2017 Adam Stark. All rights reserved.
//=================================================================================

#ifndef SampleTailExtensionHelperFunctions_h
#define SampleTailExtensionHelperFunctions_h

#include <iostream>

//=================================================================================
std::string getFileNameFromPath (std::string filePath)
{
    return filePath.substr (filePath.find_last_of("/\\") + 1);
}

//=================================================================================
int getPitchClassFromFileName (std::string fileName)
{
    // find the . from the extension
    std::size_t dotIndex = fileName.find_last_of (".");
    
    // remove the extension
    if (dotIndex != std::string::npos)
        fileName = fileName.substr (0, dotIndex);
    
    if (fileName.find ("C#") != std::string::npos)
    {
        return 1;
    }
    else if (fileName.find ("D#") != std::string::npos)
    {
        return 3;
    }
    else if (fileName.find ("F#") != std::string::npos)
    {
        return 6;
    }
    else if (fileName.find ("G#") != std::string::npos)
    {
        return 8;
    }
    else if (fileName.find ("A#") != std::string::npos)
    {
        return 10;
    }
    else if (fileName.find ("C") != std::string::npos)
    {
        return 0;
    }
    else if (fileName.find ("D") != std::string::npos)
    {
        return 2;
    }
    else if (fileName.find ("E") != std::string::npos)
    {
        return 4;
    }
    else if (fileName.find ("F") != std::string::npos)
    {
        return 5;
    }
    else if (fileName.find ("G") != std::string::npos)
    {
        return 7;
    }
    else if (fileName.find ("A") != std::string::npos)
    {
        return 9;
    }
    else if (fileName.find ("B") != std::string::npos)
    {
        return 11;
    }
    
    // oh dear, the file name doesn't contain the note name so
    // we can't work out the pitch class...
    assert (false);
    
    return -1;
}

#endif /* SampleTailExtensionHelperFunctions_h */
