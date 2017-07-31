#include "blockCal.h"
#include "math.h"


blockCal::blockCal(string file, int noBlocks)
{
  string name;
  name = "cal/"+file +".cal";
  ifstream ifile(name.c_str());
  if (!ifile.is_open()) 
    {
      cout << "could not open file " << name << endl;
    }

  for(int i=0;i<noBlocks;i++)
    {
      Equad[i]=0;
      ifile  >>  Eslope[i] >> Eint[i];
      
    }
  ifile.close();
  ifile.clear();

}

blockCal::blockCal(string name1, string name2)
{

  ifstream ifile1(name1.c_str()); //linear function of the calibration SLOPES with angle for each CsI
  ifstream ifile2(name2.c_str()); //linear function of the calibration INTERcEPT with angle for each CsI
  if (!ifile1.is_open())
    {
      cout << "could not open file " << name1 << endl;
    }
  if (!ifile2.is_open())
    {
      cout << "could not open file " << name2 << endl;
    }


  for(int i=0;i<16;i++)
    {
      ifile1 >>  Eslope_quad[i] >> Ecal_slope[i] >> Ecal_int[i];
      ifile2 >> Eshift_slope[i] >> Eshift_int[i];
    }


  ifile1.close();
  ifile1.clear();
  ifile2.close();
  ifile2.clear();

}

blockCal::~blockCal()
{

}

//could try a linear interpolation, or fit for the russian since it varies tremendously with each ring
//interpolation should work reasonably since calibration change is rather linear for the russian
//previous method of fitting calibration as function of angle and using that doesn't seem to work that well 6/20
//may need to calibrate for each ring if I want a better calibration

void blockCal::getCal(int blocko)
{
  calSlope = Eslope[blocko-1]; //c++ arrays start at 0 and the output of telescope::block(iring) doesn't
  calInt = Eint[blocko-1];
  calQuad = Equad[blocko-1];
  return;

}

float blockCal::getE(float channel, int blocko)
{
  getCal(blocko);
  energy = calQuad*pow(channel,2.)+ calSlope*channel+calInt;
  return energy;

}

float blockCal::getE(float channel, int CsIhit, float theta)
{
  calSlope = Eslope_quad[CsIhit]*theta*theta +Ecal_slope[CsIhit]*theta + Ecal_int[CsIhit]; //ignoring 2nd order polynomial
  calInt = Eshift_slope[CsIhit]*theta + Eshift_int[CsIhit];
  energy = calSlope*channel+calInt;
  return energy;

}
