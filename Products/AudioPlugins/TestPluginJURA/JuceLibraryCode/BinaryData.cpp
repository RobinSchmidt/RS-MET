/* ==================================== JUCER_BINARY_RESOURCE ====================================

   This is an auto-generated file: Any edits you make may be overwritten!

*/

namespace BinaryData
{

//================== DebugNotes.txt ==================
static const unsigned char temp_binary_data_0[] =
"For an easy workflow for debugging, build the PluginHost (in the \r\n"
"AudioApplications folder) and set up the visual studio project such that it\r\n"
"runs the just built Plugin Host. For this to work, set on the Property settings\r\n"
"under \"Debugging\" the \"Command\" field to the .exe that results from the build.\r\n"
"On my machine, the path is:\r\n"
"\r\n"
"E:\\Programming\\C++\\RS-MET\\Products\\AudioApplications\\PluginHost\\Builds\\VisualStudio2015\\x64\\Debug\\Plugin Host.exe\r\n"
"\r\n"
"\r\n"
"to build a release version with the linux makefile, open a terminal in the \r\n"
"folder containing the makefile and call:\r\n"
"\r\n"
"make CONFIG=Release\r\n"
"\r\n"
"just calling:\r\n"
"\r\n"
"make\r\n"
"\r\n"
"will build a debug version\r\n"
"\r\n"
"the following developer packages have to be installed (via the command \r\n"
"line terminal):\r\n"
"sudo apt-get install libasound2-dev \r\n"
"sudo apt-get install libfreetype6-dev\r\n"
"sudo apt-get install libx11-dev\r\n"
"sudo apt-get install libxrandr-dev \r\n"
"sudo apt-get install libxinerama-dev\r\n"
"sudo apt-get install libxcursor-dev\r\n";

const char* DebugNotes_txt = (const char*) temp_binary_data_0;

//================== ToDo.txt ==================
static const unsigned char temp_binary_data_1[] =
"";

const char* ToDo_txt = (const char*) temp_binary_data_1;


const char* getNamedResource (const char*, int&) throw();
const char* getNamedResource (const char* resourceNameUTF8, int& numBytes) throw()
{
    unsigned int hash = 0;
    if (resourceNameUTF8 != 0)
        while (*resourceNameUTF8 != 0)
            hash = 31 * hash + (unsigned int) *resourceNameUTF8++;

    switch (hash)
    {
        case 0x4358baff:  numBytes = 966; return DebugNotes_txt;
        case 0x80091737:  numBytes = 0; return ToDo_txt;
        default: break;
    }

    numBytes = 0;
    return 0;
}

const char* namedResourceList[] =
{
    "DebugNotes_txt",
    "ToDo_txt"
};

}
