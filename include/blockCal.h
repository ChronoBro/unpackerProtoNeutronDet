#ifndef blockCal_
#define blockCal_
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;


class blockCal
{
 public:
  blockCal(string file1, int noBlocks);
  blockCal(string name1,string name2);
  ~blockCal();
  void getCal(int blocko);
  float getE(float channel, int blocko);
  float getE(float channel, int CsIhit, float theta);
  double Eslope_quad[32];
  double Eslope[32];
  double Eint[32];
  double Equad[32];
  float Ecal_slope[32];
  float Ecal_int[32];
  float Eshift_slope[32];
  float Eshift_int[32];
  double calSlope;
  double calInt;
  double energy;
  double calQuad;
  
};
#endif
