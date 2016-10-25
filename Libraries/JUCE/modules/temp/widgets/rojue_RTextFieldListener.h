#ifndef rojue_RTextFieldListener_h
#define rojue_RTextFieldListener_h

namespace rojue
{

  class RTextField;

  class RTextFieldListener : public RTextEditorListener
  {

  public:

    virtual void rTextFieldChanged(RTextField *rTextField) = 0;

    virtual void rTextFieldChanged(RTextEditor &editor);

    virtual void rTextEditorTextChanged(RTextEditor &editor) { }

    virtual void rTextEditorReturnKeyPressed(RTextEditor &editor) { rTextFieldChanged(editor); }

    virtual void rTextEditorEscapeKeyPressed(RTextEditor &editor) { rTextFieldChanged(editor); }

    virtual void rTextEditorFocusLost(RTextEditor &editor) { rTextFieldChanged(editor); }

  };

}

#endif
