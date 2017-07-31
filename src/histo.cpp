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
