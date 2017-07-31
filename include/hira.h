#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TRandom.h"
#include "histo.h"
#include "calibrate.h"
#include "telescope.h"
#include "caen.h"
#include "pixels.h"


using namespace std;


struct dataE
{
  int ienergy;
  int id;
  int itime;
  float energy;
};

struct dataT
{
  int itime;
  int id;
};

struct map
{
  bool front;
  bool A;
  int itele;
};

class hira
{
 public:
  hira(TRandom* ran, histo * Histo0);
  ~hira();
  bool unpack(unsigned short*& point,int runno);
  bool unpackSi_sis(unsigned short*& point);
  bool unpackSi_adc(unsigned short*& point);
  bool unpackCsi(unsigned short*& point,int runno);
  void analyze();
  telescope **Telescope;
  void reset();

  double Blocker_e;

  int Np;
  int N6;

 private:
  TRandom * ran;
  unsigned short xmarker[3];


  unsigned short chanXLM[3][400];
  unsigned short nXLM[3];

  map Map[3][13];
  
  histo * Histo;
  caen ADC;
  caen TDC;

  dataE DataE[56];
  dataT DataT[56];
 


  calibrate * calBack;
  calibrate * calFront;
  calibrate * calLBack;
  calibrate * calLFront;
  calibrate * calCsi;
  calibrate * calBDual;
  calibrate * calFrontT;
  calibrate *calBackT;
  calibrate * calFrontTLG;
  calibrate * calBackTLG;
  calibrate * calrecal;



  int NE;
  int NT;
  int CsIM;


  //high-low correlation
  double fsumN[2][32];
  double fsumx[2][32];
  double fsumxx[2][32];
  double fsumy[2][32];
  double fsumyx[2][32];

  double bsumN[2][32];
  double bsumx[2][32];
  double bsumxx[2][32];
  double bsumy[2][32];
  double bsumyx[2][32];


  bool fred;

};
