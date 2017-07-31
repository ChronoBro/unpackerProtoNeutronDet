#include "det.h"

#include <fstream>

#include <iostream>



//*********************************************************
  /**
   * Constructor
   */
det::det(histo * Histo0)
{
  Histo = Histo0;
  ran = new TRandom;
  Detector = new detector(ran,Histo);

  Rus_count=0;
  S2_count=0;
  RusS2_count=0;

 

}
//************************************************
  /**
   * destructor
   */
det::~det()
{
  delete Detector;
  delete ran;
}
//*************************************************************
  /**
   * unpacks a physics event from the data stream
   \param point0 - pointer to location in data stream
  */
bool det::unpack(unsigned short *point,int runno)
{
  Detector->reset();


  //cout << runno << endl;

  unsigned short  words = *point++;


  Run = runno;
  
  bool stat = false;
  stat = Detector->unpack(point,runno);
  if (!stat) return stat;
  
  return stat;


}
//*********************************

void det::analyze(int event)
{
  //This is where you will analyze the events


}



