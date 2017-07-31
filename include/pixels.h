#ifndef pixels_
#define pixels_

#include <fstream>
#include <iostream>
#include <cmath>
#include "TRandom.h"

using namespace std;

struct location
{
  double x,y,z;
  double r,theta,phi;
  double deltatheta;
  double deltaphi;
  double steradian;
};

struct teleP
{
  location Location[64][48];
};

class pixels
{
 public:
  pixels();
  teleP TeleP[2];
  float thetaRus[32];
  float thetaS2[48];
  location getCenter(int itele);
  float getCsiCenter(int itele, int iCsi);
  float phi;
  float r;
  float getAngle(int itele,int ifront, int iback);
  double targetdist[2];
  float deltatheta;
  float deltaphi;
  //float ** xy;
  
};
#endif
