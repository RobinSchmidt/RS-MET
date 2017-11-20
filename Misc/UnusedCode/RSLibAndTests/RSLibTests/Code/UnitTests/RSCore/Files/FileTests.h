#ifndef RSCORE_FILETESTS_H
#define RSCORE_FILETESTS_H

#include "../Text/StringTests.h"

bool testFile(std::string &reportString);

bool testFileTextReadWrite(std::string &reportString);  // tests, if we can write a string into a file and retrieve it
                                                        // again, the string must not contain non-printable characters
bool testFileStreamReadWrite(std::string &reportString);
bool testFileWaveReadWrite(std::string &reportString);


//bool testDirectory



#endif
