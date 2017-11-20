#include "rojue_RTextField.h"
using namespace rojue;

void RTextFieldListener::rTextFieldChanged(RTextEditor &editor)

{
  RTextField* field = dynamic_cast<RTextField*> (&editor);
  rTextFieldChanged(field);
}


