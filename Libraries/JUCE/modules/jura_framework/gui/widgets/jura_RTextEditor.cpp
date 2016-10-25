//#include "rojue_RTextEditor.h"

#define SHOULD_WRAP(x, wrapwidth)       (((x) - 0.0001f) >= (wrapwidth))

// a word or space that can't be broken down any further:
struct TextAtom
{
  String atomText;
  int    width;
  uint16 numChars;

  bool  isWhitespace() const throw() { return CharacterFunctions::isWhitespace (atomText[0]); }
  bool  isNewLine() const throw() { return atomText[0] == '\r' || atomText[0] == '\n'; }
  const String getText() const throw() { return atomText; }
  const String getTrimmedText() const throw() { return atomText.substring (0, numChars); }
};

//=================================================================================================
// class UniformTextSection - a run of text with a single font and colour:

class UniformTextSection
{
public:

  UniformTextSection(const String& text, BitmapFont const* font_, const Colour& colour_) throw() 
    : font (font_), colour (colour_)
  {
    initialiseAtoms(text);
  }

  UniformTextSection (const UniformTextSection& other) throw() : font (other.font), 
    colour (other.colour)
  {
    for(int i = 0; i < other.atoms.size(); ++i)
      atoms.add (new TextAtom (*(const TextAtom*)other.atoms.getUnchecked(i)));
  }

  ~UniformTextSection() throw()
  {
    // (no need to delete the atoms, as they're explicitly deleted by the caller)
  }

  void clear() throw()
  {
    for(int i = atoms.size(); --i >= 0;)
    {
      TextAtom* const atom = getAtom(i);
      delete atom;
    }
    atoms.clear();
  }

  int getNumAtoms() const throw()
  {
    return atoms.size();
  }

  TextAtom* getAtom (const int index) const throw()
  {
    return (TextAtom*)atoms.getUnchecked (index);
  }

  void append(const UniformTextSection& other) throw()
  {
    if(other.atoms.size() > 0)
    {
      TextAtom* const lastAtom = (TextAtom*)atoms.getLast();
      int i = 0;
      if(lastAtom != 0)
      {
        if(!CharacterFunctions::isWhitespace(lastAtom->atomText.getLastCharacter()))
        {
          TextAtom* const first = other.getAtom(0);
          if(!CharacterFunctions::isWhitespace(first->atomText[0]))
          {
            lastAtom->atomText += first->atomText;
            lastAtom->numChars = (uint16)(lastAtom->numChars + first->numChars);
            lastAtom->width = 
              font->getTextPixelWidth(lastAtom->getText(), font->getDefaultKerning());
            delete first;
            ++i;
          }
        }
      }
      while(i < other.atoms.size())
      {
        atoms.add (other.getAtom(i));
        ++i;
      }
    }
  }

  UniformTextSection* split(const int indexToBreakAt) throw()
  {
    UniformTextSection* const section2 = new UniformTextSection(String::empty, font, colour);
    int index = 0;
    for(int i=0; i<atoms.size(); ++i)
    {
      TextAtom* const atom = getAtom(i);
      const int nextIndex = index + atom->numChars;
      if(index == indexToBreakAt)
      {
        int j;
        for(j = i; j < atoms.size(); ++j)
          section2->atoms.add(getAtom(j));
        for(j = atoms.size(); --j >= i;)
          atoms.remove(j);
        break;
      }
      else if(indexToBreakAt >= index && indexToBreakAt < nextIndex)
      {
        TextAtom* const secondAtom = new TextAtom();
        secondAtom->atomText = atom->atomText.substring(indexToBreakAt - index);
        secondAtom->width = 
          font->getTextPixelWidth(secondAtom->getText(), font->getDefaultKerning());
        secondAtom->numChars = (uint16)secondAtom->atomText.length();
        section2->atoms.add (secondAtom);
        atom->atomText = atom->atomText.substring(0, indexToBreakAt - index);
        atom->width    = font->getTextPixelWidth(atom->getText(), font->getDefaultKerning());
        atom->numChars = (uint16)(indexToBreakAt - index);
        int j;
        for(j=i+1; j<atoms.size(); ++j)
          section2->atoms.add(getAtom (j));
        for(j=atoms.size(); --j > i;)
          atoms.remove(j);
        break;
      }
      index = nextIndex;
    }
    return section2;
  }

  const String getAllText() const throw()
  {
    String s;
    s.preallocateBytes(getTotalLength());
    for(int i=0; i<atoms.size(); ++i)
      s += getAtom(i)->getText();
    return s;
  }
  // this supposed to be more efficient function does not work anmore sicne switching to juce 4:
  //const String getAllText() const throw()
  //{
  //  String s;
  //  //s.preallocateStorage(getTotalLength());
  //  s.preallocateBytes(getTotalLength());
  //  //tchar* endOfString = (tchar*)&(s[0]);
  //  juce_wchar endOfString = (juce_wchar)&(s[0]);
  //  for(int i=0; i<atoms.size(); ++i)
  //  {
  //    const TextAtom* const atom = getAtom(i);
  //    memcpy(endOfString, &(atom->atomText[0]), atom->numChars * sizeof (tchar));
  //    endOfString += atom->numChars;
  //  }
  //  *endOfString = 0;
  //  jassert((endOfString - (tchar*)&(s[0])) <= getTotalLength());
  //  return s;
  //}

  const String getTextSubstring(const int startCharacter, const int endCharacter) const throw()
  {
    int index    = 0;
    int totalLen = 0;
    int i;
    for(i=0; i<atoms.size(); ++i)
    {
      const TextAtom* const atom = getAtom (i);
      const int nextIndex = index + atom->numChars;
      if(startCharacter < nextIndex)
      {
        if(endCharacter <= index)
          break;
        const int start = jmax(0, startCharacter - index);
        const int end = jmin(endCharacter - index, (int)atom->numChars);
        jassert(end >= start);
        totalLen += end - start;
      }
      index = nextIndex;
    }
    String s;
    s.preallocateBytes(totalLen + 1);
    //juce_wchar* psz = (juce_wchar*)(const juce_wchar*)s;
    //juce_wchar* psz = s.getCharPointer();
    //String::CharPointerType psz = s.getCharPointer();
    index = 0;
    for(i=0; i<atoms.size(); ++i)
    {
      const TextAtom* const atom = getAtom(i);
      const int nextIndex = index + atom->numChars;
      if(startCharacter < nextIndex)
      {
        if(endCharacter <= index)
          break;
        const int start = jmax(0, startCharacter - index);
        const int len   = jmin(endCharacter - index, (int)atom->numChars) - start;
        s += atom->getText().substring(start, start+len-1);
        // probably very inefficient and also untested
        //memcpy(psz, ((const tchar*)atom->atomText) + start, len * sizeof (tchar));
        //psz += len;
        //*psz = 0;
      }
      index = nextIndex;
    }
    return s;
  }
  //const String getTextSubstring(const int startCharacter, const int endCharacter) const throw()
  //{
  //  int index    = 0;
  //  int totalLen = 0;
  //  int i;
  //  for(i=0; i<atoms.size(); ++i)
  //  {
  //    const TextAtom* const atom = getAtom (i);
  //    const int nextIndex = index + atom->numChars;
  //    if(startCharacter < nextIndex)
  //    {
  //      if(endCharacter <= index)
  //        break;
  //      const int start = jmax(0, startCharacter - index);
  //      const int end = jmin(endCharacter - index, (int)atom->numChars);
  //      jassert(end >= start);
  //      totalLen += end - start;
  //    }
  //    index = nextIndex;
  //  }
  //  String s;
  //  s.preallocateStorage (totalLen + 1);
  //  tchar* psz = (tchar*)(const tchar*)s;
  //  index = 0;
  //  for(i=0; i<atoms.size(); ++i)
  //  {
  //    const TextAtom* const atom = getAtom(i);
  //    const int nextIndex = index + atom->numChars;
  //    if(startCharacter < nextIndex)
  //    {
  //      if(endCharacter <= index)
  //        break;
  //      const int start = jmax(0, startCharacter - index);
  //      const int len   = jmin(endCharacter - index, (int)atom->numChars) - start;
  //      memcpy(psz, ((const tchar*)atom->atomText) + start, len * sizeof (tchar));
  //      psz += len;
  //      *psz = 0;
  //    }
  //    index = nextIndex;
  //  }
  //  return s;
  //}

  int getTotalLength() const throw()
  {
    int c = 0;
    for(int i = atoms.size(); --i >= 0;)
      c += getAtom(i)->numChars;
    return c;
  }

  void setFont(const BitmapFont* newFont) throw()
  {
    if(font != newFont)
    {
      font = newFont;
      for(int i = atoms.size(); --i >= 0;)
      {
        TextAtom* const atom = (TextAtom*)atoms.getUnchecked(i);
        atom->width = newFont->getTextPixelWidth(atom->getText(), newFont->getDefaultKerning());
      }
    }
  }

  juce_UseDebuggingNewOperator;

  BitmapFont const* font;
  Colour            colour;

private:

  Array<void*> atoms;

  void initialiseAtoms(const String& textToParse) throw()
  {
    int i = 0;
    const int len = textToParse.length();
    //const tchar* const text = (const tchar*)textToParse;
    //const juce_wchar* const text = (const juce_wchar*)textToParse;
    const String::CharPointerType text = textToParse.getCharPointer();
    while(i < len)
    {
      int start = i;

      // create a whitespace atom unless it starts with non-ws
      if(CharacterFunctions::isWhitespace(text[i]) && text[i] != '\r' && text[i] != '\n')
      {
        while(i < len && CharacterFunctions::isWhitespace(text[i]) && text[i] != '\r' 
          && text[i] != '\n')
          ++i;
      }
      else
      {
        if(text[i] == '\r')
        {
          ++i;
          if((i < len) && (text[i] == '\n'))
          {
            ++start;
            ++i;
          }
        }
        else if(text[i] == '\n')
          ++i;
        else
        {
          while((i < len) && !CharacterFunctions::isWhitespace (text[i]))
            ++i;
        }
      }
      TextAtom* const atom = new TextAtom();
      atom->atomText = String (text + start, i - start);
      atom->width    = font->getTextPixelWidth(atom->getText(), font->getDefaultKerning());
      atom->numChars = (uint16)(i - start);
      atoms.add (atom);
    }
  }

  const UniformTextSection& operator= (const UniformTextSection& other);
};

//=================================================================================================
// class RTextEditorIterator:

class RTextEditorIterator
{

public:

  RTextEditorIterator(const Array<void*>& sections_, const int wordWrapWidth_) throw() 
    : indexInText(0), lineY(0), lineHeight(0), maxDescent(0), atomX(0), atomRight(0), atom(0), 
    currentSection(0), sections(sections_), sectionIndex(0), atomIndex(0), 
    wordWrapWidth(wordWrapWidth_)
  {
    jassert (wordWrapWidth_ > 0);
    if(sections.size() > 0)
      currentSection = (const UniformTextSection*)sections.getUnchecked (sectionIndex);
    if(currentSection != 0)
    {
      lineHeight = currentSection->font->getFontHeight()+textEditorLineSpacing;
      maxDescent = currentSection->font->getFontDescent();
    }
  }

  RTextEditorIterator (const RTextEditorIterator& other) throw() : indexInText(other.indexInText), 
    lineY(other.lineY), lineHeight(other.lineHeight), maxDescent(other.maxDescent), 
    atomX(other.atomX), atomRight(other.atomRight), atom(other.atom),
    currentSection(other.currentSection), sections(other.sections), 
    sectionIndex(other.sectionIndex), atomIndex(other.atomIndex),
    wordWrapWidth(other.wordWrapWidth), tempAtom(other.tempAtom)
  {

  }

  ~RTextEditorIterator() throw()
  {

  }

  bool next() throw()
  {
    if(atom == &tempAtom)
    {
      const int numRemaining = tempAtom.atomText.length() - tempAtom.numChars;

      if(numRemaining > 0)
      {
        tempAtom.atomText = tempAtom.atomText.substring (tempAtom.numChars);
        atomX = 0;
        if(tempAtom.numChars > 0)
          lineY += lineHeight;
        indexInText += tempAtom.numChars;

        // \todo: this is still preliminary (used mostly 'as is' from juce::TextEditor) - we have 
        // to replace that glyphArrangement-stuff with something in BitmapFont....seems like it 
        // is used here to determine the right side of some glyph....
        GlyphArrangement g;
        Font dummyFont(14.f);
        g.addLineOfText(dummyFont, atom->getText(), 0.0f, 0.0f);
        int split;
        for(split = 0; split < g.getNumGlyphs(); ++split)
          if(SHOULD_WRAP (g.getGlyph (split).getRight(), wordWrapWidth))
            break;
        if(split > 0 && split <= numRemaining)
        {
          tempAtom.numChars = (uint16)split;
          tempAtom.width    = (int)ceil(g.getGlyph(split - 1).getRight());
          atomRight         = atomX + tempAtom.width;
          return true;
        }
      }
    }
    bool forceNewLine = false;
    if(sectionIndex >= sections.size())
    {
      moveToEndOfLastAtom();
      return false;
    }
    else if(atomIndex >= currentSection->getNumAtoms() - 1)
    {
      if(atomIndex >= currentSection->getNumAtoms())
      {
        if(++sectionIndex >= sections.size())
        {
          moveToEndOfLastAtom();
          return false;
        }
        atomIndex      = 0;
        currentSection = (const UniformTextSection*)sections.getUnchecked (sectionIndex);
        lineHeight     = jmax((int)lineHeight, currentSection->font->getFontHeight()+textEditorLineSpacing);
        maxDescent     = jmax((int)maxDescent, currentSection->font->getFontDescent());
      }
      else
      {
        const TextAtom* const lastAtom = currentSection->getAtom(atomIndex);
        if(!lastAtom->isWhitespace())
        {
          // handle the case where the last atom in a section is actually part of the same word as
          // the first atom of the next section:
          int right       = atomRight + lastAtom->width;
          int lineHeight2 = lineHeight;
          int maxDescent2 = maxDescent;
          for(int section = sectionIndex+1; section<sections.size(); ++section)
          {
            const UniformTextSection* const s = 
              (const UniformTextSection*)sections.getUnchecked(section);
            if(s->getNumAtoms() == 0)
              break;
            const TextAtom* const nextAtom = s->getAtom(0);
            if(nextAtom->isWhitespace())
              break;
            right += nextAtom->width;
            lineHeight2 = jmax((int)lineHeight2, s->font->getFontHeight()+textEditorLineSpacing);
            maxDescent2 = jmax((int)maxDescent2, s->font->getFontDescent());
            if(SHOULD_WRAP(right, wordWrapWidth))
            {
              lineHeight = lineHeight2;
              maxDescent = maxDescent2;
              forceNewLine = true;
              break;
            }
            if(s->getNumAtoms() > 1)
              break;
          }
        }
      }
    }
    if(atom != 0)
    {
      atomX = atomRight;
      indexInText += atom->numChars;
      if(atom->isNewLine())
      {
        atomX = 0;
        lineY += lineHeight;
      }
    }
    atom = currentSection->getAtom (atomIndex);
    atomRight = atomX + atom->width;
    ++atomIndex;
    if(SHOULD_WRAP(atomRight, wordWrapWidth) || forceNewLine)
    {
      if(atom->isWhitespace())
        atomRight = jmin(atomRight, wordWrapWidth); 
        // leave whitespace at the end of a line, but truncate it to avoid scrolling
      else
        return wrapCurrentAtom();
    }
    return true;
  }

  bool wrapCurrentAtom() throw()
  {
    atomRight = atom->width;
    if(SHOULD_WRAP(atomRight, wordWrapWidth))  // atom too big to fit on a line, so break it up..
    {
      tempAtom          = *atom;
      tempAtom.width    = 0;
      tempAtom.numChars = 0;
      atom              = &tempAtom;
      if(atomX > 0)
      {
        atomX = 0;
        lineY += lineHeight;
      }
      return next();
    }
    atomX = 0;
    lineY += lineHeight;
    return true;
  }

  void draw(Graphics& g, const UniformTextSection*& lastSection) const throw()
  {
    if(!atom->isWhitespace())
    {
      if(lastSection != currentSection)
      {
        lastSection = currentSection;
        g.setColour(currentSection->colour);
        //g.setFont(currentSection->font);
      }
      jassert(atom->getTrimmedText().isNotEmpty());
      drawBitmapFontText(g, atomX, lineY, atom->getTrimmedText(), (currentSection->font), 
        currentSection->colour);
    }
  }

  void drawSelection(Graphics& g, const int selectionStart, const int selectionEnd) const throw()
  {
    const int startX = indexToX(selectionStart);
    const int endX   = indexToX(selectionEnd);
    const int y      = lineY;
    const int nextY  = lineY + lineHeight;
    g.fillRect(startX, y, endX-startX, nextY-y);
  }

  void drawSelectedText(Graphics& g, const int selectionStart, const int selectionEnd, 
    const Colour& selectedTextColour) const throw()
  {
    if(!atom->isWhitespace())
    {
      /*
      Font dummyFont(16);
      GlyphArrangement ga;
      ga.addLineOfText(dummyFont, atom->getTrimmedText (passwordCharacter), atomX,
        (float) roundFloatToInt (lineY + lineHeight - maxDescent));

      if (selectionEnd < indexInText + atom->numChars)
      {
        GlyphArrangement ga2 (ga);
        ga2.removeRangeOfGlyphs (0, selectionEnd - indexInText);
        ga.removeRangeOfGlyphs (selectionEnd - indexInText, -1);

        g.setColour (currentSection->colour);
        ga2.draw (g);
      }

      if (selectionStart > indexInText)
      {
        GlyphArrangement ga2 (ga);
        ga2.removeRangeOfGlyphs (selectionStart - indexInText, -1);
        ga.removeRangeOfGlyphs (0, selectionStart - indexInText);

        g.setColour (currentSection->colour);
        ga2.draw (g);
      }

      g.setColour (selectedTextColour);
      ga.draw (g);
      */

      //drawBitmapFontText(g, atomX, roundFloatToInt(lineY), 
        //atom->getTrimmedText(passwordCharacter), *(currentSection->font), Colours::yellow);
      drawBitmapFontText(g, atomX, lineY, atom->getTrimmedText(), currentSection->font, 
        selectedTextColour);
    }
  }

  int indexToX (const int indexToFind) const throw()
  {
    if(indexToFind <= indexInText)
      return atomX;
    if(indexToFind >= indexInText + atom->numChars)
      return atomRight;
    return currentSection->font->glyphIndexToX(atom->getText(), indexToFind-indexInText, 
      currentSection->font->getDefaultKerning());
  }

  int xToIndex(const int xToFind) const throw()
  {
    if(xToFind <= atomX || atom->isNewLine())
      return indexInText;
    if(xToFind >= atomRight)
      return indexInText + atom->numChars;
    return currentSection->font->xToGlyphIndex(atom->getText(), xToFind, 
      currentSection->font->getDefaultKerning());
  }

  void updateLineHeight() throw()
  {
    int x = atomRight;
    int tempSectionIndex = sectionIndex;
    int tempAtomIndex    = atomIndex;
    const UniformTextSection* currentSection = 
      (const UniformTextSection*)sections.getUnchecked(tempSectionIndex);
    while(!SHOULD_WRAP (x, wordWrapWidth))
    {
      if(tempSectionIndex >= sections.size())
        break;
      bool checkSize = false;
      if(tempAtomIndex >= currentSection->getNumAtoms())
      {
        if(++tempSectionIndex >= sections.size())
          break;
        tempAtomIndex = 0;
        currentSection = (const UniformTextSection*)sections.getUnchecked(tempSectionIndex);
        checkSize = true;
      }
      const TextAtom* const atom = currentSection->getAtom(tempAtomIndex);
      if(atom == 0)
        break;
      x += atom->width;
      if(SHOULD_WRAP (x, wordWrapWidth) || atom->isNewLine())
        break;
      if(checkSize)
      {
        lineHeight = 
          jmax((int)lineHeight, currentSection->font->getFontHeight()+textEditorLineSpacing);
        maxDescent = jmax((int)maxDescent, currentSection->font->getFontDescent());
      }
      ++tempAtomIndex;
    }
  }

  bool getCharPosition(const int index, int& cx, int& cy, int& lineHeight_) throw()
  {
    while(next())
    {
      if(indexInText + atom->numChars >= index)
      {
        updateLineHeight();
        if(indexInText + atom->numChars > index)
        {
          cx = indexToX (index);
          cy = lineY;
          lineHeight_ = lineHeight;
          return true;
        }
      }
    }
    cx = atomX;
    cy = lineY;
    lineHeight_ = lineHeight;
    return false;
  }

  juce_UseDebuggingNewOperator;

  int   indexInText;
  int   lineY, lineHeight, maxDescent;
  int   atomX, atomRight;
  const TextAtom* atom;
  const UniformTextSection* currentSection;

private:

  const Array<void*>& sections;
  int sectionIndex, atomIndex;
  const int wordWrapWidth;
  TextAtom tempAtom;
  const RTextEditorIterator& operator= (const RTextEditorIterator&);
  void moveToEndOfLastAtom() throw()
  {
    if(atom != 0)
    {
      atomX = atomRight;
      if(atom->isNewLine())
      {
        atomX = 0;
        lineY += lineHeight;
      }
    }
  }
};

//=================================================================================================
// classes for various actions:

class RTextEditorInsertAction : public UndoableAction
{

  RTextEditor& owner;
  const String text;
  const int insertIndex, oldCaretPos, newCaretPos;
  BitmapFont const* font; // get rid
  const Colour colour;    // get rid

  RTextEditorInsertAction (const RTextEditorInsertAction&);
  const RTextEditorInsertAction& operator= (const RTextEditorInsertAction&);

public:

  RTextEditorInsertAction(RTextEditor& owner_, const String& text_, const int insertIndex_, 
    BitmapFont const* font_, const Colour& colour_, const int oldCaretPos_, 
    const int newCaretPos_) throw()
    : owner(owner_), text(text_), insertIndex(insertIndex_), oldCaretPos(oldCaretPos_), 
    newCaretPos(newCaretPos_), font(font_), colour (colour_)
  {

  }

  ~RTextEditorInsertAction()
  {

  }

  bool perform()
  {
    owner.insert(text, insertIndex, font, colour, 0, newCaretPos);
    return true;
  }

  bool undo()
  {
    owner.remove(insertIndex, insertIndex + text.length(), 0, oldCaretPos);
    return true;
  }

  int getSizeInUnits() { return text.length() + 16; }

};


class RTextEditorRemoveAction : public UndoableAction
{

  RTextEditor& owner;
  const int startIndex, endIndex, oldCaretPos, newCaretPos;
  Array<void*> removedSections;

  RTextEditorRemoveAction (const RTextEditorRemoveAction&);
  const RTextEditorRemoveAction& operator= (const RTextEditorRemoveAction&);

public:

  RTextEditorRemoveAction(RTextEditor& owner_, const int startIndex_, const int endIndex_, 
    const int oldCaretPos_, const int newCaretPos_, 
    const Array<void*>& removedSections_) throw() 
    : owner(owner_), startIndex(startIndex_), endIndex(endIndex_), oldCaretPos(oldCaretPos_), 
    newCaretPos(newCaretPos_), removedSections(removedSections_)
  {

  }

  ~RTextEditorRemoveAction()
  {
    for(int i = removedSections.size(); --i >= 0;)
    {
      UniformTextSection* const section = (UniformTextSection*)removedSections.getUnchecked (i);
      section->clear();
      delete section;
    }
  }

  bool perform()
  {
    owner.remove(startIndex, endIndex, 0, newCaretPos);
    return true;
  }

  bool undo()
  {
    owner.reinsert(startIndex, removedSections);
    owner.moveCursorTo(oldCaretPos, false);
    return true;
  }

  int getSizeInUnits()
  {
    int n = 0;
    for(int i = removedSections.size(); --i >= 0;)
    {
      UniformTextSection* const section = (UniformTextSection*)removedSections.getUnchecked (i);
      n += section->getTotalLength();
    }
    return n + 16;
  }

};

//=================================================================================================
// class RTextHolderComponent (what is this class good for?):

class RTextHolderComponent : public Component, public Timer
{

  RTextEditor* const owner;

  RTextHolderComponent (const RTextHolderComponent&);
  const RTextHolderComponent& operator= (const RTextHolderComponent&);

public:

  RTextHolderComponent(RTextEditor* const owner_) : owner(owner_)
  {
    setWantsKeyboardFocus(false);
    setInterceptsMouseClicks(false, true);
  }

  ~RTextHolderComponent()
  {

  }

  void paint(Graphics& g)            { owner->drawContent (g); }
  void timerCallback()               { owner->timerCallbackInt(); }
  //const MouseCursor getMouseCursor() { return owner->getMouseCursor(); }
  MouseCursor getMouseCursor() { return owner->getMouseCursor(); }

};

//=================================================================================================
// class RTextEditorViewport :

class RTextEditorViewport : public Viewport
{

  RTextEditor* const owner;
  int lastWordWrapWidth;

  RTextEditorViewport (const RTextEditorViewport&);
  const RTextEditorViewport& operator= (const RTextEditorViewport&);

public:

  RTextEditorViewport(RTextEditor* const owner_) : owner(owner_), lastWordWrapWidth(0)
  {

  }

  ~RTextEditorViewport()
  {

  }

  virtual void visibleAreaChanged(const Rectangle<int>& newVisibleArea)
  {
    const int wordWrapWidth = owner->getWordWrapWidth();
    if(wordWrapWidth != lastWordWrapWidth)
    {
      lastWordWrapWidth = wordWrapWidth;
      owner->updateTextHolderSize();
    }
  }
  //void visibleAreaChanged(int, int, int, int)
  //{
  //  const int wordWrapWidth = owner->getWordWrapWidth();
  //  if(wordWrapWidth != lastWordWrapWidth)
  //  {
  //    lastWordWrapWidth = wordWrapWidth;
  //    owner->updateTextHolderSize();
  //  }
  //}

};

//=================================================================================================
// class RTextEditor:

const int flashSpeedIntervalMs = 380;
const int textChangeMessageId  = 0x10003001;
const int returnKeyMessageId   = 0x10003002;
const int escapeKeyMessageId   = 0x10003003;
const int focusLossMessageId   = 0x10003004;

RTextEditor::RTextEditor(const String& name)
  : RWidget(name),
  borderSize(1, 1, 1, 3),
  readOnly(false),
  multiline(false),
  wordWrap(false),
  returnKeyStartsNewLine(false),
  caretVisible(true),
  popupMenuEnabled(true),
  selectAllTextWhenFocused(false),
  scrollbarVisible(true),
  wasFocused(false),
  caretFlashState(true),
  keepCursorOnScreen(true),
  tabKeyUsed(false),
  menuActive(false),
  cursorX(0),
  cursorY(0),
  cursorHeight(0),
  maxTextLength(0),
  selectionStart(0),
  selectionEnd(0),
  leftIndent(4),
  topIndent(4),
  lastTransactionTime(0),
  currentFont(0),
  totalNumChars(0),
  caretPosition(0),
  dragType(notDragging)
{
  //currentFont = &boldFont10px;
  currentFont = &BitmapFontRoundedBoldA10D0::instance;
  setOpaque(true);
  addAndMakeVisible(viewport = new RTextEditorViewport (this));
  viewport->setViewedComponent(textHolder = new RTextHolderComponent (this));
  viewport->setWantsKeyboardFocus(false);
  viewport->setScrollBarsShown(false, false);
  setMouseCursor(MouseCursor::IBeamCursor);
  setWantsKeyboardFocus(true);
}

RTextEditor::~RTextEditor()
{
  clearInternal(0);
  delete viewport;
}

void RTextEditor::newTransaction() throw()
{
  lastTransactionTime = Time::getApproximateMillisecondCounter();
  undoManager.beginNewTransaction();
}

void RTextEditor::doUndoRedo(const bool isRedo)
{
  if(!isReadOnly())
  {
    if((isRedo) ? undoManager.redo() : undoManager.undo())
    {
      scrollToMakeSureCursorIsVisible();
      repaint();
      textChanged();
    }
  }
}

void RTextEditor::setMultiLine(const bool shouldBeMultiLine, const bool shouldWordWrap)
{
  multiline = shouldBeMultiLine;
  wordWrap  = shouldWordWrap && shouldBeMultiLine;
  setScrollbarsShown(scrollbarVisible);
  viewport->setViewPosition(0, 0);
  resized();
  scrollToMakeSureCursorIsVisible();
}

bool RTextEditor::isMultiLine() const throw()
{
  return multiline;
}

void RTextEditor::setScrollbarsShown(bool enabled) throw()
{
  scrollbarVisible = enabled;
  enabled = enabled && isMultiLine();
  viewport->setScrollBarsShown (enabled, enabled);
}

void RTextEditor::setReadOnly(const bool shouldBeReadOnly)
{
  readOnly = shouldBeReadOnly;
  enablementChanged();
}

bool RTextEditor::isReadOnly() const throw()
{
  return readOnly || !isEnabled();
}

void RTextEditor::setReturnKeyStartsNewLine(const bool shouldStartNewLine)
{
  returnKeyStartsNewLine = shouldStartNewLine;
}

void RTextEditor::setTabKeyUsedAsCharacter(const bool shouldTabKeyBeUsed) throw()
{
  tabKeyUsed = shouldTabKeyBeUsed;
}

void RTextEditor::setPopupMenuEnabled(const bool b) throw()
{
  popupMenuEnabled = b;
}

void RTextEditor::setSelectAllWhenFocused(const bool b) throw()
{
  selectAllTextWhenFocused = b;
}

const BitmapFont* RTextEditor::getFont() const throw()
{
  return currentFont;
}

void RTextEditor::setFont(const BitmapFont* newFont) throw()
{
  currentFont = newFont;
  scrollToMakeSureCursorIsVisible();
}

void RTextEditor::applyFontToAllText(const BitmapFont* newFont)
{
  currentFont = newFont;
  const Colour overallColour(getTextColour());
  for(int i = sections.size(); --i >= 0;)
  {
    UniformTextSection* const uts = (UniformTextSection*)sections.getUnchecked(i);
    uts->setFont(newFont);
    uts->colour = overallColour;
  }
  coalesceSimilarSections();
  updateTextHolderSize();
  scrollToMakeSureCursorIsVisible();
  repaint();
}

void RTextEditor::colourChanged()
{
  setOpaque(getBackgroundColour().isOpaque());
  repaint();
}

void RTextEditor::setCaretVisible(const bool shouldCaretBeVisible) throw()
{
  caretVisible = shouldCaretBeVisible;
  if(shouldCaretBeVisible)
    textHolder->startTimer (flashSpeedIntervalMs);
  setMouseCursor(shouldCaretBeVisible ? MouseCursor::IBeamCursor : MouseCursor::NormalCursor);
}

void RTextEditor::setInputRestrictions(const int maxLen, const String& chars) throw()
{
  maxTextLength     = jmax(0, maxLen);
  allowedCharacters = chars;
}

void RTextEditor::setTextToShowWhenEmpty(const String& text, const Colour& colourToUse) throw()
{
  textToShowWhenEmpty = text;
  colourForTextWhenEmpty = colourToUse;
}

void RTextEditor::setScrollBarThickness(const int newThicknessPixels)
{
  viewport->setScrollBarThickness(newThicknessPixels);
}

//void RTextEditor::setScrollBarButtonVisibility(const bool buttonsVisible)
//{
//  viewport->setScrollBarButtonVisibility(buttonsVisible);
//}

void RTextEditor::clear()
{
  clearInternal (0);
  updateTextHolderSize();
  undoManager.clearUndoHistory();
}

void RTextEditor::setText(const String& newText, const bool sendTextChangeMessage)
{
  const int newLength = newText.length();
  if(newLength != getTotalNumChars() || getText() != newText)
  {
    const int  oldCursorPos   = caretPosition;
    const bool cursorWasAtEnd = oldCursorPos >= getTotalNumChars();
    clearInternal(0);
    insert(newText, 0, currentFont, getTextColour(), 0, caretPosition);

    // if you're adding text with line-feeds to a single-line text editor, it ain't gonna look 
    // right!
    jassert (multiline || !newText.containsAnyOf("\r\n"));
    if(cursorWasAtEnd && !isMultiLine())
      moveCursorTo (getTotalNumChars(), false);
    else
      moveCursorTo (oldCursorPos, false);
    if(sendTextChangeMessage)
      textChanged();
    repaint();
  }
  updateTextHolderSize();
  scrollToMakeSureCursorIsVisible();
  undoManager.clearUndoHistory();
}

void RTextEditor::textChanged() throw()
{
  updateTextHolderSize();
  postCommandMessage (textChangeMessageId);
}

void RTextEditor::returnPressed()
{
  postCommandMessage(returnKeyMessageId);
}

void RTextEditor::escapePressed()
{
  postCommandMessage(escapeKeyMessageId);
}

void RTextEditor::addListener(RTextEditorListener* const newListener) throw()
{
  jassert(newListener != 0);
  if(newListener != 0)
    listeners.add(newListener);
}

void RTextEditor::removeListener(RTextEditorListener* const listenerToRemove) throw()
{
  listeners.removeFirstMatchingValue(listenerToRemove);
}

void RTextEditor::timerCallbackInt()
{
  const bool newState = (!caretFlashState) && !isCurrentlyBlockedByAnotherModalComponent();
  if(caretFlashState != newState)
  {
    caretFlashState = newState;
    if(caretFlashState)
      wasFocused = true;
    if(caretVisible && hasKeyboardFocus (false) && !isReadOnly())
      repaintCaret();
  }
  const unsigned int now = Time::getApproximateMillisecondCounter();
  if(now > lastTransactionTime + 200)
    newTransaction();
}

void RTextEditor::repaintCaret()
{
  if(!getSpecialColour1().isTransparent())
  {
    repaint(borderSize.getLeft()+textHolder->getX()+leftIndent+cursorX-1,  // x
      borderSize.getTop()+textHolder->getY()+topIndent+cursorY-1,            // y
      4, cursorHeight+2);                                                    // w, h
  }
}

void RTextEditor::repaintText(int textStartIndex, int textEndIndex)
{
  if(textStartIndex > textEndIndex && textEndIndex > 0)
    swapVariables (textStartIndex, textEndIndex);
  int x = 0, y = 0, lh = currentFont->getFontHeight()+textEditorLineSpacing;
  const int wordWrapWidth = getWordWrapWidth();
  if(wordWrapWidth > 0)
  {
    RTextEditorIterator i(sections, wordWrapWidth);
    i.getCharPosition(textStartIndex, x, y, lh);
    const int y1 = (int)y;
    int y2;
    if(textEndIndex >= 0)
    {
      i.getCharPosition(textEndIndex, x, y, lh);
      y2 = (int)(y + lh * 2.0f);
    }
    else
      y2 = textHolder->getHeight();
    textHolder->repaint (0, y1, textHolder->getWidth(), y2 - y1);
  }
}

void RTextEditor::moveCaret(int newCaretPos) throw()
{
  if(newCaretPos < 0)
    newCaretPos = 0;
  else if(newCaretPos > getTotalNumChars())
    newCaretPos = getTotalNumChars();
  if(newCaretPos != getCaretPosition())
  {
    repaintCaret();
    caretFlashState = true;
    caretPosition   = newCaretPos;
    textHolder->startTimer(flashSpeedIntervalMs);
    scrollToMakeSureCursorIsVisible();
    repaintCaret();
  }
}

void RTextEditor::setCaretPosition(const int newIndex) throw()
{
  moveCursorTo (newIndex, false);
}

int RTextEditor::getCaretPosition() const throw()
{
  return caretPosition;
}

int RTextEditor::getWordWrapWidth() const throw()
{
  return (wordWrap) ? (viewport->getMaximumVisibleWidth() - leftIndent - leftIndent / 2) : INT_MAX;
}

void RTextEditor::updateTextHolderSize() throw()
{
  const int wordWrapWidth = getWordWrapWidth();
  if(wordWrapWidth > 0)
  {
    int maxWidth = 0;
    RTextEditorIterator i(sections, wordWrapWidth);
    while(i.next())
      maxWidth = jmax(maxWidth, i.atomRight);
    const int w = leftIndent + maxWidth;
    const int h = topIndent  + jmax(i.lineY+i.lineHeight, 
      currentFont->getFontHeight()+textEditorLineSpacing);
    textHolder->setSize(w + 1, h + 1);
  }
}

int RTextEditor::getTextWidth() const throw()
{
  return textHolder->getWidth();
}

int RTextEditor::getTextHeight() const throw()
{
  return textHolder->getHeight();
}

void RTextEditor::setIndents (const int newLeftIndent, const int newTopIndent) throw()
{
  leftIndent = newLeftIndent;
  topIndent = newTopIndent;
}

void RTextEditor::setBorder(const BorderSize<int>& border) throw()
{
  borderSize = border;
  resized();
}

const BorderSize<int> RTextEditor::getBorder() const throw()
{
  return borderSize;
}

void RTextEditor::setScrollToShowCursor(const bool shouldScrollToShowCursor) throw()
{
  keepCursorOnScreen = shouldScrollToShowCursor;
}

void RTextEditor::scrollToMakeSureCursorIsVisible() throw()
{
  cursorHeight = currentFont->getFontHeight(); // (in case the text is empty and the call below doesn't set this value)
  getCharPosition(caretPosition, cursorX, cursorY, cursorHeight);
  if(keepCursorOnScreen)
  {
    int x = viewport->getViewPositionX();
    int y = viewport->getViewPositionY();
    const int relativeCursorX = cursorX - x;
    const int relativeCursorY = cursorY - y;
    if(relativeCursorX < jmax(1, proportionOfWidth (0.05f)))
      x += relativeCursorX - proportionOfWidth (0.2f);
    else if(relativeCursorX > jmax(0, viewport->getMaximumVisibleWidth() - (wordWrap ? 2 : 10)))
      x += relativeCursorX + (isMultiLine() ? proportionOfWidth (0.2f) : 10) 
      - viewport->getMaximumVisibleWidth();
    x = jlimit (0, jmax (0, textHolder->getWidth() + 8 - viewport->getMaximumVisibleWidth()), x);
    if(!isMultiLine())
      y = (getHeight() - textHolder->getHeight() - topIndent) / -2;
    else
    {
      const int curH = cursorHeight;
      if(relativeCursorY < 0)
        y = jmax (0, relativeCursorY + y);
      else if(relativeCursorY > jmax(0, viewport->getMaximumVisibleHeight() - topIndent - curH))
        y += relativeCursorY + 2 + curH + topIndent - viewport->getMaximumVisibleHeight();
    }
    viewport->setViewPosition (x, y);
  }
}

void RTextEditor::moveCursorTo(const int newPosition, const bool isSelecting) throw()
{
  if(isSelecting)
  {
    moveCaret (newPosition);
    const int oldSelStart = selectionStart;
    const int oldSelEnd   = selectionEnd;
    if(dragType == notDragging)
    {
      if(abs(getCaretPosition() - selectionStart) < abs(getCaretPosition() - selectionEnd))
        dragType = draggingSelectionStart;
      else
        dragType = draggingSelectionEnd;
    }
    if(dragType == draggingSelectionStart)
    {
      selectionStart = getCaretPosition();
      if(selectionEnd < selectionStart)
      {
        swapVariables(selectionStart, selectionEnd);
        dragType = draggingSelectionEnd;
      }
    }
    else
    {
      selectionEnd = getCaretPosition();
      if(selectionEnd < selectionStart)
      {
        swapVariables(selectionStart, selectionEnd);
        dragType = draggingSelectionStart;
      }
    }
    jassert(selectionStart <= selectionEnd);
    jassert(oldSelStart <= oldSelEnd);
    repaintText(jmin (oldSelStart, selectionStart), jmax(oldSelEnd, selectionEnd));
  }
  else
  {
    dragType = notDragging;
    if(selectionEnd > selectionStart)
      repaintText(selectionStart, selectionEnd);
    moveCaret(newPosition);
    selectionStart = getCaretPosition();
    selectionEnd   = getCaretPosition();
  }
}

int RTextEditor::getTextIndexAt(const int x, const int y) throw()
{
  return indexAtPosition((x + viewport->getViewPositionX() - leftIndent), 
    (y + viewport->getViewPositionY() - topIndent));
}

void RTextEditor::insertTextAtCursor(String newText)
{
  if(allowedCharacters.isNotEmpty())
    newText = newText.retainCharacters(allowedCharacters);
  if(!isMultiLine())
    newText = newText.replaceCharacters("\r\n", "  ");
  else
    newText = newText.replace("\r\n", "\n");
  const int newCaretPos = selectionStart + newText.length();
  const int insertIndex = selectionStart;
  remove (selectionStart, selectionEnd, &undoManager, newCaretPos);
  if(maxTextLength > 0)
    newText = newText.substring (0, maxTextLength - getTotalNumChars());
  if(newText.isNotEmpty())
    insert(newText, insertIndex, currentFont, getTextColour(), &undoManager, newCaretPos);
  textChanged();
}

void RTextEditor::setHighlightedRegion(int startPos, int numChars) throw()
{
  moveCursorTo(startPos, false);
  moveCursorTo(startPos + numChars, true);
}

void RTextEditor::copy()
{
  const String selection(getTextSubstring (selectionStart, selectionEnd));
  if(selection.isNotEmpty())
    SystemClipboard::copyTextToClipboard(selection);
}

void RTextEditor::paste()
{
  if(!isReadOnly())
  {
    const String clip(SystemClipboard::getTextFromClipboard());
    if(clip.isNotEmpty())
      insertTextAtCursor(clip);
  }
}

void RTextEditor::cut()
{
  if(!isReadOnly())
  {
    moveCaret(selectionEnd);
    insertTextAtCursor(String::empty);  // is this the right behavior for cut?
  }
}

void RTextEditor::drawContent(Graphics& g)
{
  const int wordWrapWidth = getWordWrapWidth();
  if(wordWrapWidth > 0)
  {
    g.setOrigin(leftIndent, topIndent);
    const Rectangle<int> clip(g.getClipBounds());
    Colour selectedTextColour = getTextColour();
    RTextEditorIterator i(sections, wordWrapWidth);

    while(i.lineY + 200.0 < clip.getY() && i.next())
    {
    }

    if(selectionStart < selectionEnd)
    {
      g.setColour(getHandleColour().withMultipliedAlpha(hasKeyboardFocus (true) ? 1.0f : 0.5f));

      //selectedTextColour = findColour (highlightedTextColourId);
      //selectedTextColour = RWidget::highlightColour;
      //g.setColour(Colours::red); // test

      RTextEditorIterator i2(i);
      while(i2.next() && i2.lineY < clip.getBottom())
      {
        i2.updateLineHeight();
        if(i2.lineY + i2.lineHeight >= clip.getY() && selectionEnd >= i2.indexInText
          && selectionStart <= i2.indexInText + i2.atom->numChars)
        {
          i2.drawSelection (g, selectionStart, selectionEnd);
        }
      }
    }

    const UniformTextSection* lastSection = 0;
    while(i.next() && i.lineY < clip.getBottom())
    {
      i.updateLineHeight();
      if(i.lineY + i.lineHeight >= clip.getY())
      {
        if(selectionEnd >= i.indexInText && selectionStart <= i.indexInText + i.atom->numChars)
        {
          i.drawSelectedText (g, selectionStart, selectionEnd, selectedTextColour);
          lastSection = 0;
        }
        else
          i.draw(g, lastSection);
      }
    }
  }
}

void RTextEditor::paint(Graphics& g)
{
  g.fillAll(getBackgroundColour());
}

void RTextEditor::paintOverChildren(Graphics& g)
{
  if(caretFlashState && hasKeyboardFocus(false) && caretVisible && !isReadOnly())
  {
    g.setColour(getHandleColour());
    g.fillRect(borderSize.getLeft() + textHolder->getX() + leftIndent + cursorX,
      borderSize.getTop() + textHolder->getY() + topIndent + cursorY,
      2, cursorHeight);
  }
  if(textToShowWhenEmpty.isNotEmpty() && (!hasKeyboardFocus(false)) && getTotalNumChars() == 0)
  {
    g.setColour(getTextColour());
    //g.setFont(getFont());
    if(isMultiLine())
      drawBitmapFontText(g, 0, 0, textToShowWhenEmpty, currentFont, getTextColour());
    else
      drawBitmapFontText(g, leftIndent, topIndent, textToShowWhenEmpty, currentFont, 
        getTextColour());
  }
  g.setColour(getOutlineColour());
  g.drawRect(0, 0, getWidth(), getHeight(), 2);
}

void RTextEditor::mouseDown(const MouseEvent& e)
{
  beginDragAutoRepeat(100);
  newTransaction();
  if(wasFocused || !selectAllTextWhenFocused)
  {
    if(!(popupMenuEnabled && e.mods.isPopupMenu()))
      moveCursorTo(getTextIndexAt (e.x, e.y), e.mods.isShiftDown());
    else
    {
      PopupMenu m;  // use RPopUpMenu
      addPopupMenuItems(m, &e);
      menuActive       = true;
      const int result = m.show();
      menuActive       = false;
      if(result != 0)
        performPopupMenuAction(result);
    }
  }
}

void RTextEditor::mouseDrag(const MouseEvent& e)
{
  if(wasFocused || !selectAllTextWhenFocused)
  {
    if(!(popupMenuEnabled && e.mods.isPopupMenu()))
      moveCursorTo(getTextIndexAt(e.x, e.y), true);
  }
}

void RTextEditor::mouseUp(const MouseEvent& e)
{
  newTransaction();
  textHolder->startTimer(flashSpeedIntervalMs);
  if(wasFocused || !selectAllTextWhenFocused)
  {
    if(!(popupMenuEnabled && e.mods.isPopupMenu()))
      moveCaret(getTextIndexAt(e.x, e.y));
  }
  wasFocused = true;
}

void RTextEditor::mouseDoubleClick(const MouseEvent& e)
{
  int tokenEnd   = getTextIndexAt(e.x, e.y);
  int tokenStart = tokenEnd;
  if(e.getNumberOfClicks() > 3)
  {
    tokenStart = 0;
    tokenEnd = getTotalNumChars();
  }
  else
  {
    const String t(getText());
    const int totalLength = getTotalNumChars();
    while(tokenEnd < totalLength)
    {
      if(CharacterFunctions::isLetterOrDigit(t[tokenEnd]))
        ++tokenEnd;
      else
        break;
    }
    tokenStart = tokenEnd;
    while(tokenStart > 0)
    {
      if(CharacterFunctions::isLetterOrDigit(t[tokenStart-1]))
        --tokenStart;
      else
        break;
    }
    if(e.getNumberOfClicks() > 2)
    {
      while(tokenEnd < totalLength)
      {
        if(t[tokenEnd] != '\r' && t[tokenEnd] != '\n')
          ++tokenEnd;
        else
          break;
      }
      while(tokenStart > 0)
      {
        if(t[tokenStart-1] != '\r' && t[tokenStart-1] != '\n')
          --tokenStart;
        else
          break;
      }
    }
  }
  moveCursorTo(tokenEnd, false);
  moveCursorTo(tokenStart, true);
}

void RTextEditor::mouseWheelMove(const MouseEvent &e, const MouseWheelDetails &wheel)
{
  if(!viewport->useMouseWheelMoveIfNeeded(e, wheel))
    Component::mouseWheelMove(e, wheel);
}
//void RTextEditor::mouseWheelMove(const MouseEvent& e, float wheelIncrementX, float wheelIncrementY)
//{
//  if(!viewport->useMouseWheelMoveIfNeeded(e, wheelIncrementX, wheelIncrementY))
//    Component::mouseWheelMove(e, wheelIncrementX, wheelIncrementY);
//}

bool RTextEditor::keyPressed(const KeyPress& key)
{
  if(isReadOnly() && key != KeyPress('c', ModifierKeys::commandModifier, 0))
    return false;
  const bool moveInWholeWordSteps = key.getModifiers().isCtrlDown() || key.getModifiers().isAltDown();
  if(key.isKeyCode(KeyPress::leftKey) || key.isKeyCode(KeyPress::upKey))
  {
    newTransaction();
    int newPos;
    if(isMultiLine() && key.isKeyCode(KeyPress::upKey))
      newPos = indexAtPosition (cursorX, cursorY - 1);
    else if(moveInWholeWordSteps)
      newPos = findWordBreakBefore(getCaretPosition());
    else
      newPos = getCaretPosition() - 1;
    moveCursorTo (newPos, key.getModifiers().isShiftDown());
  }
  else if(key.isKeyCode(KeyPress::rightKey) || key.isKeyCode(KeyPress::downKey))
  {
    newTransaction();
    int newPos;
    if(isMultiLine() && key.isKeyCode(KeyPress::downKey))
      newPos = indexAtPosition(cursorX, cursorY + cursorHeight + 1);
    else if(moveInWholeWordSteps)
      newPos = findWordBreakAfter(getCaretPosition());
    else
      newPos = getCaretPosition() + 1;
    moveCursorTo(newPos, key.getModifiers().isShiftDown());
  }
  else if(key.isKeyCode(KeyPress::pageDownKey) && isMultiLine())
  {
    newTransaction();
    moveCursorTo(indexAtPosition(cursorX, cursorY + cursorHeight + viewport->getViewHeight()), 
      key.getModifiers().isShiftDown());
  }
  else if(key.isKeyCode(KeyPress::pageUpKey) && isMultiLine())
  {
    newTransaction();
    moveCursorTo(indexAtPosition(cursorX, cursorY - viewport->getViewHeight()), 
      key.getModifiers().isShiftDown());
  }
  else if(key.isKeyCode (KeyPress::homeKey))
  {
    newTransaction();
    if(isMultiLine() && !moveInWholeWordSteps)
      moveCursorTo(indexAtPosition(0, cursorY), key.getModifiers().isShiftDown());
    else
      moveCursorTo(0, key.getModifiers().isShiftDown());
  }
  else if(key.isKeyCode(KeyPress::endKey))
  {
    newTransaction();
    if(isMultiLine() && !moveInWholeWordSteps)
      moveCursorTo(indexAtPosition(textHolder->getWidth(), cursorY), 
        key.getModifiers().isShiftDown());
    else
      moveCursorTo(getTotalNumChars(), key.getModifiers().isShiftDown());
  }
  else if(key.isKeyCode(KeyPress::backspaceKey))
  {
    if(moveInWholeWordSteps)
      moveCursorTo(findWordBreakBefore(getCaretPosition()), true);
    else
    {
      if(selectionStart == selectionEnd && selectionStart > 0)
        --selectionStart;
    }
    cut();
  }
  else if(key.isKeyCode(KeyPress::deleteKey))
  {
    if(key.getModifiers().isShiftDown())
      copy();
    if(selectionStart == selectionEnd && selectionEnd < getTotalNumChars())
      ++selectionEnd;
    cut();
  }
  else if(key == KeyPress('c', ModifierKeys::commandModifier, 0)
    || key == KeyPress(KeyPress::insertKey, ModifierKeys::ctrlModifier, 0))
  {
    newTransaction();
    copy();
  }
  else if(key == KeyPress('x', ModifierKeys::commandModifier, 0))
  {
    newTransaction();
    copy();
    cut();
  }
  else if(key == KeyPress('v', ModifierKeys::commandModifier, 0)
    || key == KeyPress(KeyPress::insertKey, ModifierKeys::shiftModifier, 0))
  {
    newTransaction();
    paste();
  }
  else if(key == KeyPress('z', ModifierKeys::commandModifier, 0))
  {
    newTransaction();
    doUndoRedo(false);
  }
  else if(key == KeyPress('y', ModifierKeys::commandModifier, 0))
  {
    newTransaction();
    doUndoRedo(true);
  }
  else if(key == KeyPress('a', ModifierKeys::commandModifier, 0))
  {
    newTransaction();
    moveCursorTo(getTotalNumChars(), false);
    moveCursorTo(0, true);
  }
  else if(key == KeyPress::returnKey)
  {
    newTransaction();
    if(returnKeyStartsNewLine)
      insertTextAtCursor("\n");
    else
      returnPressed();
  }
  else if(key.isKeyCode(KeyPress::escapeKey))
  {
    newTransaction();
    moveCursorTo(getCaretPosition(), false);
    escapePressed();
  }
  else if(key.getTextCharacter() >= ' ' || (tabKeyUsed && (key.getTextCharacter() == '\t')))
  {
    insertTextAtCursor (String::charToString (key.getTextCharacter()));
    lastTransactionTime = Time::getApproximateMillisecondCounter();
  }
  else
    return false;

  return true;
}


const int baseMenuItemID = 0x7fff0000;
void RTextEditor::addPopupMenuItems(PopupMenu& m, const MouseEvent*)
{
  const bool writable = !isReadOnly();

  m.addItem (baseMenuItemID + 1, TRANS("cut"), writable);
  m.addItem (baseMenuItemID + 2, TRANS("copy"), selectionStart < selectionEnd);
  m.addItem (baseMenuItemID + 3, TRANS("paste"), writable);
  m.addItem (baseMenuItemID + 4, TRANS("delete"), writable);
  m.addSeparator();
  m.addItem (baseMenuItemID + 5, TRANS("select all"));
  m.addSeparator();
  m.addItem (baseMenuItemID + 6, TRANS("undo"), undoManager.canUndo());
  m.addItem (baseMenuItemID + 7, TRANS("redo"), undoManager.canRedo());
}

void RTextEditor::performPopupMenuAction(const int menuItemID)
{
  switch(menuItemID)
  {
  case baseMenuItemID + 1:
    copy();
    cut();
    break;
  case baseMenuItemID + 2:
    copy();
    break;
  case baseMenuItemID + 3:
    paste();
    break;
  case baseMenuItemID + 4:
    cut();
    break;
  case baseMenuItemID + 5:
    moveCursorTo (getTotalNumChars(), false);
    moveCursorTo (0, true);
    break;
  case baseMenuItemID + 6:
    doUndoRedo (false);
    break;
  case baseMenuItemID + 7:
    doUndoRedo (true);
    break;
  default:
    break;
  }
}

void RTextEditor::focusGained(FocusChangeType)
{
  newTransaction();
  caretFlashState = true;
  if(selectAllTextWhenFocused)
  {
    moveCursorTo (0, false);
    moveCursorTo (getTotalNumChars(), true);
  }
  repaint();
  if(caretVisible)
    textHolder->startTimer(flashSpeedIntervalMs);
}

void RTextEditor::focusLost(FocusChangeType)
{
  newTransaction();
  wasFocused = false;
  textHolder->stopTimer();
  caretFlashState = false;
  postCommandMessage (focusLossMessageId);
  repaint();
}

void RTextEditor::resized()
{
  viewport->setBoundsInset(borderSize);
  viewport->setSingleStepSizes(16, currentFont->getFontHeight()+textEditorLineSpacing);
  updateTextHolderSize();
  if(!isMultiLine())
    scrollToMakeSureCursorIsVisible();
  else
  {
    cursorHeight = currentFont->getFontHeight(); // in case the text is empty and the call below doesn't set this value
    getCharPosition(caretPosition, cursorX, cursorY, cursorHeight);
  }
}

void RTextEditor::handleCommandMessage(const int commandId)
{
  // my new:
  Component::BailOutChecker checker(this);
  for(int i=listeners.size(); --i >= 0;)
  {
    RTextEditorListener* const tl = (RTextEditorListener*)listeners[i];
    if(tl != 0)
    {
      switch(commandId)
      {
      case textChangeMessageId: tl->rTextEditorTextChanged(*this);      break;
      case returnKeyMessageId:  tl->rTextEditorReturnKeyPressed(*this); break;
      case escapeKeyMessageId:  tl->rTextEditorEscapeKeyPressed(*this); break;
      case focusLossMessageId:  tl->rTextEditorFocusLost(*this);        break;
      default:                  jassertfalse;                           break;
      }
    }
  }
}

void RTextEditor::enablementChanged()
{
  setMouseCursor(MouseCursor (isReadOnly() ? 
    MouseCursor::NormalCursor : MouseCursor::IBeamCursor));
  repaint();
}

void RTextEditor::clearInternal(UndoManager* const um) throw()
{
  remove(0, getTotalNumChars(), um, caretPosition);
}

void RTextEditor::insert(const String& text, const int insertIndex, BitmapFont const* font, 
  const Colour& colour, UndoManager* const um, const int caretPositionToMoveTo) throw()
{
  if(text.isNotEmpty())
  {
    if(um != 0)
      um->perform(new RTextEditorInsertAction(*this, text, insertIndex, font, colour, 
        caretPosition, caretPositionToMoveTo));
    else
    {
      repaintText(insertIndex, -1); // must do this before and after changing the data, in case a line gets moved due to word wrap
      int index = 0;
      int nextIndex = 0;
      for(int i = 0; i < sections.size(); ++i)
      {
        nextIndex = index + ((UniformTextSection*)sections.getUnchecked(i))->getTotalLength();
        if(insertIndex == index)
        {
          sections.insert(i, new UniformTextSection(text, font, colour));
          break;
        }
        else if(insertIndex > index && insertIndex < nextIndex)
        {
          splitSection (i, insertIndex - index);
          sections.insert (i + 1, new UniformTextSection(text, font, colour));
          break;
        }
        index = nextIndex;
      }
      if(nextIndex == insertIndex)
        sections.add(new UniformTextSection(text, font, colour));
      coalesceSimilarSections();
      totalNumChars = -1;
      moveCursorTo(caretPositionToMoveTo, false);
      repaintText(insertIndex, -1);
    }
  }
}

void RTextEditor::reinsert(const int insertIndex, const Array<void*>& sectionsToInsert) throw()
{
  int index = 0;
  int nextIndex = 0;
  for(int i=0; i<sections.size(); ++i)
  {
    nextIndex = index + ((UniformTextSection*)sections.getUnchecked(i))->getTotalLength();
    if(insertIndex == index)
    {
      for(int j = sectionsToInsert.size(); --j >= 0;)
        sections.insert(i, new UniformTextSection(
          *(UniformTextSection*)sectionsToInsert.getUnchecked(j)));
      break;
    }
    else if(insertIndex > index && insertIndex < nextIndex)
    {
      splitSection (i, insertIndex - index);
      for(int j = sectionsToInsert.size(); --j >= 0;)
        sections.insert (i + 1, new UniformTextSection(
          *(UniformTextSection*)sectionsToInsert.getUnchecked(j)));
      break;
    }
    index = nextIndex;
  }
  if(nextIndex == insertIndex)
  {
    for(int j = 0; j < sectionsToInsert.size(); ++j)
      sections.add(new UniformTextSection(*(UniformTextSection*)sectionsToInsert.getUnchecked(j)));
  }
  coalesceSimilarSections();
  totalNumChars = -1;
}

void RTextEditor::remove(const int startIndex, int endIndex, UndoManager* const um, 
  const int caretPositionToMoveTo) throw()
{
  if(endIndex > startIndex)
  {
    int index = 0;
    for(int i=0; i<sections.size(); ++i)
    {
      const int nextIndex = index + ((UniformTextSection*)sections[i])->getTotalLength();
      if(startIndex > index && startIndex < nextIndex)
      {
        splitSection(i, startIndex - index);
        --i;
      }
      else if(endIndex > index && endIndex < nextIndex)
      {
        splitSection(i, endIndex - index);
        --i;
      }
      else
      {
        index = nextIndex;
        if(index > endIndex)
          break;
      }
    }
    index = 0;
    if(um != 0)
    {
      Array<void*> removedSections;
      for(int i=0; i<sections.size(); ++i)
      {
        if(endIndex <= startIndex)
          break;
        UniformTextSection* const section = (UniformTextSection*)sections.getUnchecked (i);
        const int nextIndex = index + section->getTotalLength();
        if(startIndex <= index && endIndex >= nextIndex)
          removedSections.add (new UniformTextSection (*section));
        index = nextIndex;
      }
      um->perform(new RTextEditorRemoveAction(*this, startIndex, endIndex, caretPosition, 
        caretPositionToMoveTo, removedSections));
    }
    else
    {
      for(int i=0; i<sections.size(); ++i)
      {
        if(endIndex <= startIndex)
          break;
        UniformTextSection* const section = (UniformTextSection*)sections.getUnchecked(i);
        const int nextIndex = index + section->getTotalLength();
        if(startIndex <= index && endIndex >= nextIndex)
        {
          sections.remove(i);
          endIndex -= (nextIndex - index);
          section->clear();
          delete section;
          --i;
        }
        else
          index = nextIndex;
      }
      coalesceSimilarSections();
      totalNumChars = -1;
      moveCursorTo(caretPositionToMoveTo, false);
      repaintText(startIndex, -1);
    }
  }
}

const String RTextEditor::getText() const throw()
{
  String t;
  for(int i=0; i<sections.size(); ++i)
    t += ((const UniformTextSection*)sections.getUnchecked(i))->getAllText();
  return t;
}

const String RTextEditor::getTextSubstring(const int startCharacter, 
  const int endCharacter) const throw()
{
  String t;
  int index = 0;
  for(int i=0; i<sections.size(); ++i)
  {
    const UniformTextSection* const s = (const UniformTextSection*)sections.getUnchecked(i);
    const int nextIndex = index + s->getTotalLength();
    if(startCharacter < nextIndex)
    {
      if(endCharacter <= index)
        break;
      const int start = jmax(index, startCharacter);
      t += s->getTextSubstring(start - index, endCharacter - index);
    }
    index = nextIndex;
  }
  return t;
}

const String RTextEditor::getHighlightedText() const throw()
{
  return getTextSubstring(getHighlightedRegionStart(), getHighlightedRegionStart() 
    + getHighlightedRegionLength());
}

int RTextEditor::getTotalNumChars() throw()
{
  if(totalNumChars < 0)
  {
    totalNumChars = 0;
    for(int i = sections.size(); --i >= 0;)
      totalNumChars += ((const UniformTextSection*)sections.getUnchecked(i))->getTotalLength();
  }
  return totalNumChars;
}

bool RTextEditor::isEmpty() const throw()
{
  if(totalNumChars != 0)
  {
    for(int i = sections.size(); --i >= 0;)
      if(((const UniformTextSection*)sections.getUnchecked(i))->getTotalLength() > 0)
        return false;
  }
  return true;
}

void RTextEditor::getCharPosition(const int index, int& cx, int& cy, int& lineHeight) const throw()
{
  const int wordWrapWidth = getWordWrapWidth();
  if(wordWrapWidth > 0 && sections.size() > 0)
  {
    RTextEditorIterator i(sections, wordWrapWidth);
    i.getCharPosition (index, cx, cy, lineHeight);
  }
  else
  {
    cx = cy = 0;
    lineHeight = currentFont->getFontHeight()+textEditorLineSpacing;
  }
}

int RTextEditor::indexAtPosition(const int x, const int y) throw()
{
  const int wordWrapWidth = getWordWrapWidth();
  if(wordWrapWidth > 0)
  {
    RTextEditorIterator i(sections, wordWrapWidth);
    while(i.next())
    {
      if(i.lineY + getHeight() > y)
        i.updateLineHeight();
      if(i.lineY + i.lineHeight > y)
      {
        if(i.lineY > y)
          return jmax (0, i.indexInText - 1);
        if(i.atomX >= x)
          return i.indexInText;
        if(x < i.atomRight)
          return i.xToIndex (x);
      }
    }
  }
  return getTotalNumChars();
}

static int getCharacterCategory(const juce_wchar character) throw()
{
  return CharacterFunctions::isLetterOrDigit(character) ? 
    2 : (CharacterFunctions::isWhitespace(character) ? 0 : 1);
}
//static int getCharacterCategory (const tchar character) throw()
//{
//  return CharacterFunctions::isLetterOrDigit(character) ? 2 : (CharacterFunctions::isWhitespace(character) ? 0 : 1);
//}

int RTextEditor::findWordBreakAfter(const int position) const throw()
{
  const String t(getTextSubstring(position, position + 512));
  const int totalLength = t.length();
  int i = 0;
  while(i < totalLength && CharacterFunctions::isWhitespace(t[i]))
    ++i;
  const int type = getCharacterCategory (t[i]);
  while(i < totalLength && type == getCharacterCategory(t[i]))
    ++i;
  while(i < totalLength && CharacterFunctions::isWhitespace(t[i]))
    ++i;
  return position + i;
}

int RTextEditor::findWordBreakBefore(const int position) const throw()
{
  if(position <= 0)
    return 0;
  const int startOfBuffer = jmax(0, position - 512);
  const String t(getTextSubstring (startOfBuffer, position));
  int i = position - startOfBuffer;
  while(i > 0 && CharacterFunctions::isWhitespace(t[i - 1]))
    --i;
  if(i > 0)
  {
    const int type = getCharacterCategory (t[i - 1]);
    while(i > 0 && type == getCharacterCategory (t[i - 1]))
      --i;
  }
  jassert(startOfBuffer + i >= 0);
  return startOfBuffer + i;
}

void RTextEditor::splitSection (const int sectionIndex, const int charToSplitAt) throw()
{
  jassert(sections[sectionIndex] != 0);
  sections.insert(sectionIndex + 1, 
    ((UniformTextSection*)sections.getUnchecked(sectionIndex))->split(charToSplitAt));
}

void RTextEditor::coalesceSimilarSections() throw()
{
  for(int i=0; i<sections.size() - 1; ++i)
  {
    UniformTextSection* const s1 = (UniformTextSection*)(sections.getUnchecked (i));
    UniformTextSection* const s2 = (UniformTextSection*)(sections.getUnchecked (i + 1));
    if(s1->font == s2->font && s1->colour == s2->colour)
    {
      s1->append(*s2);
      sections.remove(i+1);
      delete s2;
      --i;
    }
  }
}
