namespace rosic {
namespace Sampler {

void rsRemoveLineComments(std::string& str, char commentStart)
{
  // Removes everything between the given character that starts a comment and the next linebreak 
  // except when the commentStart character is part of a string which is detected by looking for 
  // a ". The commentStart character itself is also removed but the linebreak is kept.

  // Our states:                     read-index is in...
  static const int inText    = 1; // ...regular text
  static const int inString  = 2; // ...a string
  static const int inComment = 3; // ...a comment
  static const int inRhs     = 4; // ...a right hand side of an assignment

  int state = inText;             // current state
  size_t ri = 0;                  // read index
  size_t wi = 0;                  // write index

  while(ri < str.size()) {

    // Check if we need to change our state:
    char c = str[ri];
    if(c == '"' && state != inComment) {       // Detect start and end of strings
      if(state != inString)
        state = inString;
      else
        state = inText;   }
    else if(c == '=' && state == inText)       // Detect start of assignment
      state = inRhs;
    else if(c == ' ' && state == inRhs)        // Detect end of assignment
      state = inText;
    else if(c == commentStart) {               // Detect start of comments
      if(state != inString && state != inRhs)
        state = inComment; }
    else if(c == '\n') {                       // Detect end of comments
      state = inText; }

    // Advance to next input character, possibly copying current input char into output:
    if(state == inComment)
      ri++;
    else {
      str[wi] = str[ri];
      ri++;
      wi++; }
  }

  str.resize(wi);

  // Notes:
  // The '/' character marks the start of a line comment and the comment extends for the rest of
  // the line, i.e. until the next newline '\n' is encountered. But: If the '/' occurs inside the
  // path of a filename, it does not start a comment. We solve this by introducing a state
  // inRhs in which the '/' is ignored (i.e. doesn't put us into comment-state)
  // 

  // Bugs:
  // -using ' ' to detect the end of an assignment does not work: when we have a comment slash
  //  followed by a space, we'll jump out of inComment state - but shouldn't

  // Maybe this is overkill - at least for the wavefile names. They are actually not given in 
  // quotes. However, maybe later we will need to handle strings anyway - for example, for the
  // formula opcode. ...that also means we do not yet have test coverage for sfzs  containing 
  // strings ...what if a newline occurs within a string? We don't handle this case and maybe
  // we don't need to. Do we?

  // Wait: this will not work when the filenames containing the slash are not wrapped in quotes!
  // Maybe we should also the '=' as symbol for entering the "inString" state (which should then be 
  // renamed) and the ' ' for leaving that state. Or maybe introduce the 4th state "inRhs" that 
  // behaves similar to "inString"...but actually, "inString" won't be needed then anymore.
  // Maybe take 3 chars as input: the symbol '/' that starts a comment, the symbol '=' that starts 
  // a section within which '/' is ignored, the symbol ' ' that ends the ignore section. Make a 
  // unit test containing filenames using subdirectories...ok done...added the inRhs stuff - it now
  // needs tests
}
// this needs more unit tests - it doesn't work correctly yet

void rsReplaceCharacter(std::string& str, char oldChar, char newChar)
{
  for(size_t i = 0; i < str.size(); i++) {
    if(str[i] == oldChar)
      str[i] = newChar; }
}

void rsRemoveRepeats(std::string& s, char c)
{
  // Removes repeated occurrences of character c, i.e. replaces each sequence of multiple cs with a 
  // single c.

  if(s.size() < 2)
    return;
  size_t i = 0, j = 1;                  // i: read index, j: write index
  for(i = 1; i < s.size(); i++) {
    if(!(s[i] == c && s[i-1] == c)) {
      s[j] = s[i]; 
      j++; }}
  s.resize(j);
}
// move into rosic, write unit test

int firstDigit(const std::string& str)
{
  for(int i = 0; i < (int) str.size(); i++)
    if(isDigit(str[i]))
      return i;
  return -1;
}

int lastDigitInSeq(const std::string& str, int startIndex)  // remove the "InSeq" from name
{
  RAPT::rsAssert(isDigit(str[startIndex]));
  int i = startIndex;
  while(i < ((int)str.size())-1) {
    if(!isDigit(str[i+1]))
      return i;
    i++; }
  return i;
}

int lastNonDigit(const std::string& str, int startIndex)
{
  for(int i = startIndex; i < (int)str.size(); i++)
  {
    if(isDigit(str[i]))
      return i-1;
  }
  return (int)str.size() - 1;
}
// needs unit test, be sure to test also edge cases (empty string, no digits present, 1st char is 
// digit etc.)
// maybe a function firstDigit would be more useful? but no - we need something that works with 
// strings that contain no digits, too

int suffixStart(const std::string& str)
{
  for(int i = (int)str.length()-1; i >= 0; --i)
    if(isDigit(str[i]) || str[i] == '.')   // numbers end in a digit or a dot
      return i+1;
  return 0;
}
// maybe rename, needs tests

int parseNaturalNumber(const std::string& str, int startIndex, int endIndex)
{
  int num    = 0;
  int scaler = 1;
  for(int i = endIndex; i >= startIndex; i--)
  {
    num += scaler * charToInt(str[i]);
    scaler *= 10;
  }
  return num;
}
// maybe rename to rsStringToUnsigned, needs unit tests

float rsStringToFloat(const std::string& str) 
{ 
  //rsAssert(rsIsFloatNumber(str)); 
  // Maybe we should do this to avoid more severe crashes in stof when a non-number string is
  // passed. stof may throw std::bad_alloc exceptions in such cases. Or maybe we should just
  // try and catch them?

  //return std::stof(str);
  // Old. It's not a good idea to use stof raw. It doesn't seem to be error tolerant at all. It
  // may crash the plugin, if the input is not a valid float ..hmm..well...it just raises an assert 
  // but still - oh - no - it actually also triggers a std::bad_alloc exception. Maybe use a custom
  // parser....maybe include the rsBigInt/rsBigFloat classes from RSLib into RAPT. They contain 
  // such a parser. But maybe for that, we'll also move the RSLib::rsString code into RAPT where it
  // more or less will parallel/obviate rosic::rsString. Or Maybe move the rsBigInt/Float code into
  // rosic and on not (yet) templatize it on the integer type. If really needed, this can be done 
  // later.
  // ...hmm...not sure what's best

  // Maybe not the best solution yet, but good enough for now:
  try
  {
    return std::stof(str);
  }
  catch(std::invalid_argument)
  {
    RAPT::rsError("Failed to parse string in rsStringToFloat");
    return RS_NAN(float);  // Or maybe we should return 0? Or maybe the caller should check against
    // NaN, and turn it into 0, if it wants to? API-wise, it seems most sensible to indeed return 
    // NaN from a general stringToFloat function in case of an invalid argument. But then it would 
    // necessarily be inconsistent with a corresponding strToint function because int doesn't have 
    // a concept of NaN
  }
}

// Not declared in header file because we currently only use it internally:
bool hasImplicitFirstGroup(const std::string& code)
{
  // An implicit first group definition (starting at the begin of the string) occurs, if the first
  // <region> header appears before the first <group> header.

  size_t firstGroupStart  = code.find("<group>");
  if(firstGroupStart == string::npos)
    return true;
    // No <group> header was found at all. In this case, there is just one group and that group is
    // indeed implicit.

  size_t firstRegionStart = code.find("<region>");
  if(firstRegionStart == string::npos)
    return false;
    // There is a <group> header but no <region> header. In this case, we have an explicit but 
    // empty (first) group.

  return firstGroupStart > firstRegionStart;
  // When both <group> and <region> headers actually do exist, we need to compare their positions. 
  // If the first region header comes befeore the first group header, the first group is implicit.
}

void findSfzGroup(const std::string& code, int groupIndex, int* startIndex, int* endIndex)
{
  *startIndex = -1;
  *endIndex   = -1;
  if(code.empty())
    return;

  std::string pattern = "<group>";       // Search pattern
  size_t L = pattern.length();
  int foundIndex = -1;
  size_t start = 0;

  // Handle special case when the first group header is missing and therefore implicit:
  if(hasImplicitFirstGroup(code)) {
    *startIndex = 0;
    foundIndex  = 0; }

  // Find the start:
  while(foundIndex < groupIndex) {
    start = code.find(pattern, start);
    if(start != string::npos) {          // npos is returned when no match was found
      foundIndex++;
      *startIndex = (int) start;
      start += L;
      if(foundIndex == groupIndex)
        break; }
    else {
      *startIndex = -1;
      return;  }}           // No group with given index could be found in the given sfz code

  // Find the end. The end is defined to be either the last character in the string or the 
  // character immediately before the subsequent <group> header:
  start = code.find(pattern, start);
  if(start != string::npos)
    *endIndex = (int) start - 1;
  else
    *endIndex = (int) code.length() - 1;
}

void findSfzRegion(const std::string& code, int regionIndex, int searchStart, int searchEnd,
  int* startIndex, int* endIndex)
{
  *startIndex = -1;
  *endIndex   = -1;
  if(code.size() < searchStart)  // maybe <=?
    return;

  std::string pattern = "<region>";       // Search pattern
  size_t L = pattern.length();
  int foundIndex = -1;
  size_t start = (size_t) searchStart;

  // Find the start:
  while(foundIndex < regionIndex && (int) start <= searchEnd) {
    start = code.find(pattern, start);
    if(start != string::npos) {          // npos is returned when no match was found
      foundIndex++;
      *startIndex = (int) start;
      start += L;
      if(foundIndex == regionIndex)
        break; }
    else {
      *startIndex = -1;
      return;  }}  

  // Handle case when we left the loop due to the 2nd loop condition. In this case the caller is
  // searching for a region with an index higher than present in the code:
  if(foundIndex < regionIndex) {
    *startIndex = -1;
    return; }
    // ToDo: document, why we don't have/need a similar handling in findSfzGroup. I think, it's 
    // because there, the 2nd loop condition doesn't exist

  // Find the end. The end is defined to be either the searchEnd or the character immediately 
  // before the subsequent <region> header:
  start = code.find(pattern, start);
  if(start != string::npos)
    *endIndex = (int) start - 1;
  else
    *endIndex = searchEnd;
}
// needs tests
// Can we get rid of the code duplication? Maybe like this:
// -Pass the search pattern as argument
// -Let searchstart/End be optional parameters that defualt to 0 and code.length()-1
// -I think the conditionals:
//    if(code.size() < searchStart) 
//    while(foundIndex < regionIndex && (int) start <= searchEnd)
//  in findSfzRegion are compatible with those in findSfzGroup. 
// -The if(hasImplicitFirstGroup(code)) in findSfzGroup could be problematic, though






}
}