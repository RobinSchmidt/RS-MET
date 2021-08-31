// Experiments and functions for dealing with meshes, such as computing numerical derivatives, 
// which is an importants step in numerical PDE solvers. These functions were previously located in
// MathExperiments.cpp, but because the amount of code grew so large that at some point, i moved 
// them into this separate file. Their declarations, however, are still in MathExperiments.h 
// because  didn't want to create yet another header file and at some point i want to geht of these
// .h files anyway and instead directly include the .cpp files in the test driver files. Less 
// files, less confusion.