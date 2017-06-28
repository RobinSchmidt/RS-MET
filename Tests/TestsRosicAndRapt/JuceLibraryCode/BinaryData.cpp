/* ==================================== JUCER_BINARY_RESOURCE ====================================

   This is an auto-generated file: Any edits you make may be overwritten!

*/

namespace BinaryData
{

//================== ReadMe.txt ==================
static const unsigned char temp_binary_data_0[] =
"This is the place for experimental prototype code and code that is currently under construction\r\n"
"and not yet production ready. When it reaches a reasonable level of maturity, it may be dragged \r\n"
"over into the library. Thereby it must possibly be adapted (templatized, optimized, etc.). In \r\n"
"some cases, the prototype can be kept here for reference, testing and education, too.";

const char* ReadMe_txt = (const char*) temp_binary_data_0;


const char* getNamedResource (const char*, int&) throw();
const char* getNamedResource (const char* resourceNameUTF8, int& numBytes) throw()
{
    unsigned int hash = 0;
    if (resourceNameUTF8 != 0)
        while (*resourceNameUTF8 != 0)
            hash = 31 * hash + (unsigned int) *resourceNameUTF8++;

    switch (hash)
    {
        case 0x4de88c9f:  numBytes = 376; return ReadMe_txt;
        default: break;
    }

    numBytes = 0;
    return 0;
}

const char* namedResourceList[] =
{
    "ReadMe_txt"
};

}
