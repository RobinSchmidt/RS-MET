#include "rojue_DrawingThread.h"
using namespace rojue;

DrawingThread* DrawingThread::soleInstance = NULL;

DrawingThread::DrawingThread() : TimeSliceThread(String(T("DrawingThread")))
{
  startThread(1);
}

DrawingThread* DrawingThread::getInstance()
{
  if( soleInstance == NULL )
    soleInstance = new DrawingThread();
  return soleInstance;
}




/*
DrawingThread drawingThread;

DrawingThread::DrawingThread() : TimeSliceThread(String(T("DrawingThread")))
{
  startThread(1);
}
*/








