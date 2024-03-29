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
    if(c == '"') {                             // Detect start and end of strings
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

}
}