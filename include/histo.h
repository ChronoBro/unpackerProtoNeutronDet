#ifndef histo_
#define histo_
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TFile.h"
#include "TDirectory.h"

using namespace std;

class histo
{
 protected:

  TFile * file; //!< output root file


  //CsI
  TDirectoryFile *dirCsI; //!< directory for the CsI info
  TDirectory * dirCsIRaw; //!< CsI energies


 public:
  histo();                  //!< constructor
  ~histo(){};
  void write(); //!< write the root spectra to file


  int NSPies;
  int NSRings;
  int NRPies;
  int NRRings;

  int NCsI;

  //S2



  TH1I ** ECsI;

};
#endif
