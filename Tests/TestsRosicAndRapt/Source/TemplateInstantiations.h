#ifndef TemplateInstantiationsForTests_h
#define TemplateInstantiationsForTests_h

/** This file does some explicit template instantiations for templates from RAPT which we want to 
test but for which we do not (yet) have an instantiation anywhere else (like in e.g. 
rosic/basics/rosic_TemplateInstantiations.cpp"). As soon as some template gets instantiated there, 
we need to delete the corresponding instantiation here to avoid a double-instantiation with the 
associated linker error. */


//template class RAPT::rsModularInteger<int>;

// Oh - wait - this file is not needed - we do these instantiations in the rs_testing module!




#endif