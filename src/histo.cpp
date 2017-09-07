#include "histo.h"


histo::histo()
{

  ostringstream outstring;
  string name;

  
  NCsI = 2;

  //create root file
  file = new TFile ("sort.root","RECREATE");
  file->cd();




   //CsI 
  dirCsI = new TDirectoryFile("CsI","CsI");
  EDet0_gain = new TH1D("EDet0_gain","EDet0_gain",1024,0,4095);
  EDet1_gain = new TH1D("EDet1_gain","EDet1_gain",1024,0,4095);

  
  dirCsIRaw = dirCsI->mkdir("CsIRaw","raw");
  
 
 




  ECsI = new TH1I*[NCsI];

  for(int icsi = 0;icsi <NCsI;icsi++)
    {
      outstring.str("");
      outstring << "ECsI_" << icsi;
      name = outstring.str();
      dirCsIRaw->cd();
      ECsI[icsi] = new TH1I(name.c_str(),"",1024,0,4095);


    }

  relDifference = new TH1D("relDif","relDif",200,-1,1);
  timeDif = new TH1D("timeDif","timeDif",1024,0,4095);
  relVtime = new TH2D("relVtime","relVtime",200,-1,1,200,-10,10);
  gainMatch = new TH2D("gainMatch","gainMatch",1024,0,4095,1024,0,4095);
  timeDifCalibrated = new TH1D("timeDifCalibrated","timeDifCalibrated",200,-10,10);

}
//*********************************************
void histo::write()
{
  file->Write();

  cout << "histo written" << endl;
  /*
    for (int i=0;i<Ntele;i++)
    {
    delete red[i];
    }
    delete [] red;
  */
}
