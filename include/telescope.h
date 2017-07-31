#ifndef telescope_
#define telescope_

#include "elist.h"
#include "histo.h"
#include "pixels.h"
#include "pid.h"
#include <string>
#include <sstream>
#include "solution.h"
#include "calibrate.h"
#include "loss.h"
#include "TRandom.h"
#include "TGraph.h"
#include "blockCal.h"
class telescope
{
 public:


  static bool const relativity;
  telescope(TRandom* ran,int id0, histo * Histo0);
  ~telescope();
  elist Front;
  elist Back;
  elist Pie;
  elist Ring;
  elist Csi;
  elist Silicon;
  

  void analyze(int);
  void reset();
  int multiHitCsi();
  int multiHit();
  void load(int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int);
  void load(int,int);
  void Addfake();
  void Addfake2();
  void getMomentum();
  int Block(int);

  int id;
  float fenergy;
  float benergy;
  float penergy;
  float renergy;
  float Sienergy;
  float CsIenergy;
  float flow;
  float blow;
  float xhit;
  float yhit;
  float theta;
  float phi;
  float triton_par1;
  float alpha_quad;
  float triton_quad;
  double bins[82];

  int fhit;
  int bhit;
  int phit;
  int rhit;
  int CsIhit;

  int Np;
  int N6;

  int multPie;
  int multRing;

  int Event;

 
  solution Solution[10];
  int Nsolution;

  int CsIMap[64];
  int target=0;
  int Csicheck;
  float triton_shift;
  double theta_use;
  int block;
 

 private:
  TRandom * ran;
  histo * Histo;
  pixels  Pixels;

  int NestDim;
  void loop(int);
  int NestArray[50];
  int arrayD[50];
  int arrayB[50];
  float deMin;
  int dstripMin;

  int FrontLow[4];
  int FrontHigh[4];
  int BackLow[4];
  int BackHigh[4];
  int fitCal;
  pid  * Pid[32];


  CLoss * Loss[36];

  calibrate * calCsi;
  calibrate * calCsid;
  calibrate * calCsit;

  calibrate * calCsiHe3;
  calibrate * calCsiA;

  calibrate * calCsiLi6;

  calibrate * calCsiBe7;

  blockCal * blockCalCsid[32];

  blockCal * blockCalCsiA[32];

  blockCal * blockCalCsidRus;
  blockCal * blockCalCsiARus;

  float check_slope[32];
  float check_intercept[32];

};

#endif
