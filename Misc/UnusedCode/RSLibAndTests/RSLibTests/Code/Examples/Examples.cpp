#include "RSAudio/Synthesis/ModalExamples.h"
#include "RSGraphics/Rendering/FractalExamples.h"

/*
void createModalExamples()
{
  //createModalFilterExamples();

  //createModalFilterBankExamples();

  //createBass1();
  //createGong1();
  createPluck1();
}

void createFractalExamples()
{
  createDefaultFractal();
  //createMandelbrotSet1();
  //createBuddhaBrot1();
  int dummy = 0;
}
*/

int main(int argc, char** argv)
{
  // modal synthesis sample-sets:
  createBass1();
  //createGong1();
  //createPluck1();

  // sorting-algorithm synthesis
  //createInsertionSortSound();

  // fractals:
  //createDefaultFractal();
  //createMandelbrotSet1();
  //createBuddhaBrot1();


  //createModalExamples();
  //createFractalExamples();


  std::cout << "\n";
  std::cout << "All examples have been rendered.";

  // ckeck for memory leaks

  if( detectMemoryLeaks() )
    std::cout << "\n\n!!! Memory leaks detected !!! \n";

  getchar();
}

