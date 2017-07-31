#ifndef det_
#define det_
#include <iostream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "detector.h"
#include "histo.h"
#include "correl.h"
#include "doppler.h"

using namespace std;

/**
 *!\brief detector readout
 *
 * responsible for unpacking the data stream for physics events
 */


class det
{
 public:
  det(histo * Histo);
  ~det();
  bool unpack(unsigned short *point,int runno);
  detector *Detector;
  histo * Histo;
  TRandom * ran;
  void analyze(int event);

  int Event;

  int Run;

  correl Correl;
  void corr_7Li();
 

  int Nalphat;

  int CsImult;

  int N_IAS;
  int N7Li_ta;
  int N7Li;
  float xg,yg,zg;
  float xL,yL,zL;
  float thetag,phig;
  float thetaL,phiL,thetarel;
  float dot;
  float mag;
  int Rus_count;
  int S2_count;
  int RusS2_count;



 private:
  doppler *Doppler;
  int please;

};
#endif
