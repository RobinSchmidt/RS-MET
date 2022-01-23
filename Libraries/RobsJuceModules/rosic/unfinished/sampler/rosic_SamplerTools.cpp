namespace rosic {
namespace Sampler {

void rsRemoveLineComments(std::string& str, char commentStart)
{
  // Removes everything between the given character that starts a comment and the next linebreak 
  // except when the commentStart character is part of a string which is detected by looking for 
  // a ". The commentStart character itself is also removed but the linebreak is kept.

                                  // read-index is in...
  static const int inText    = 1; // ...regular text
  static const int inString  = 2; // ...a string
  static const int inComment = 3; // ...a comment
  int state = inText;             // current state
  size_t ri = 0;                  // read index
  size_t wi = 0;                  // write index
  //size_t length = str.size();

  while(ri < str.size()) {

    // Check if we need to change our state:
    char c = str[ri];
    if(c == '"') {                 // Detect start and end of strings
      if(state != inString)
        state = inString;
      else
        state = inText;   }
    else if(c == commentStart) {  // Detect start of comments
      if(state != inString)
        state = inComment; }
    else if(c == '\n') {          // Detect end of comments
      state = inText; }

    // Advance to next input character, possibly copying current into output:
    if(state == inComment)
      ri++;
    else {
      str[wi] = str[ri];
      ri++;
      wi++; }
  }

  str.resize(wi);

  // Maybe this is overkill - at least for the wavefile names. They are actually not given in 
  // quotes. However, maybe later we will need to handle strings anyway - for example, for the
  // function or formula opcode. ...that also means we do not yet have test coverage for sfzs 
  // containing strings
}

void rsReplaceCharacter(std::string& str, char oldChar, char newChar)
{
  for(size_t i = 0; i < str.size(); i++) {
    if(str[i] == oldChar)
      str[i] = newChar; }
}

void rsRemoveRepeats(std::string& s, char c)
{
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

}
}