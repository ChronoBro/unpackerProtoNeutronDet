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
  
  int channels = 3;
  double relDif;
  double left;
  double right;

  double gain = 2085./1740.;
  gain = 1.;
  double thresholdLeft  = 200.;
  double thresholdRight  = 200.;
  double time;
  double offset = 50.;
  //double E[2];
  //double ID[2];
  double t0 = 2600.; //channel of time at 0 position
  double tSlope = 215.; //channels/ns
  double tCal;
  
  for(int i=0;i<channels;i++)
      {
	if(Detector->DataE[i].id == 0)
	  left = (double)Detector->DataE[i].ienergy*gain + offset;	    
	else if (Detector->DataE[i].id == 1)
	  right = (double)Detector->DataE[i].ienergy;
	else if(Detector->DataE[i].id == 32)
	  time = (double)Detector->DataE[i].ienergy;
      }

  tCal = (time-t0)/tSlope;
  
  Histo->EDet1_gain->Fill(right);
  Histo->EDet0_gain->Fill(left);
  Histo->timeDif->Fill(time);
  Histo->gainMatch->Fill(left,right);
  Histo->timeDifCalibrated->Fill(tCal);
  
  if(left > thresholdLeft && right > thresholdRight){
  
    
    relDif = (left-right)/(left+right);    
    // cout << Detector->DataE[0].ienergy << endl;
    // cout << Detector->DataE[1].ienergy << endl;    
    Histo->relDifference->Fill(relDif);
    Histo->relVtime->Fill(relDif,tCal);
  }
  

}



