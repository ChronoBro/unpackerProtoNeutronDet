#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TGraph.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;


double gauss(float x,float mean,float sig)
{
double const pi=3.14159;
 return exp(-pow((x-mean)/sig,2)/2.)/(sqrt(2.*pi)*sig);
}



double Peaks(double *x , double *par)
{

  int const Npeak = 9;
  double const intensity[9]={1.0,0.1194,0.1711,0.7077,
			     0.0166,0.1310,0.848,0.2310,0.7690};
  double const energy[9]={3.183,5.106,5.144,5.156,
			  5.388,5.443,5.486,5.763,5.805};
 
  double tot = 0.;
  double e = par[0] + par[1]*x[0];
  for (int i=0;i<Npeak;i++)
    {
      if(i ==0)
	tot+= par[3]*intensity[i]*gauss(e,energy[i],par[2]);
      else if(i >0 && i < 4)
	tot+= par[4]*intensity[i]*gauss(e,energy[i],par[2]);
      else if(i>3 && i < 7)
	tot+= par[5]*intensity[i]*gauss(e,energy[i],par[2]);
      else 
	tot+= par[6]*intensity[i]*gauss(e,energy[i],par[2]);
    }
  return tot;
}

int main()
{

  float const conS =.1/2.35;
  ofstream fout("cal/sRingN.cal");
  ofstream fwhm("cal/fwhmsRing.dat");
  TFile f("sort.root");
  TCanvas* canvas[1];
  int Ntele = 1;
  int Nstrip = 48;
  ostringstream outstring;
  TH2I frame("frame","",10,2.8,6.5,10,0,500);
  frame.SetStats(kFALSE);

  double xx[48];
  double yy[48];

  TF1 *func = new TF1("fit",Peaks,2.5,7,7);
  double para[7];
  ifstream file("cal/S2Rings.cal");
  float intercept, slope;
  int i1,i2;
  string name;

  TH1F con("con","",500,0,10);
  for (int itele=0;itele<Ntele;itele++)
    {
      outstring.str("");
      outstring << "Rus"<<itele;
      name = outstring.str();
      canvas[itele] = new TCanvas(name.c_str());
      canvas[itele]->Divide(6,8);
      for (int istrip =0;istrip<Nstrip;istrip++)
        {
      
          canvas[itele]->cd(istrip+1);
          file >> i1 >> i2  >> slope >> intercept;


          outstring.str("");
          outstring << "S2/S2ringsC/S2RingC"<<"_"<<istrip;
          string name = outstring.str();
          cout <<  name << endl;
          TH1I * hist = (TH1I*) f.Get(name.c_str());
	  frame.Draw();


          hist->SetStats(kFALSE);
          hist->GetXaxis()->SetRangeUser(3.,6.5);
          con.GetXaxis()->SetRangeUser(3.,6.5);
	  for (int i=1;i<=400;i++)
	    for (int j=1;j<400;j++)
	    {
              float deltax = hist->GetBinCenter(i)-con.GetBinCenter(j);
	      if (fabs(deltax) > 10.*conS)continue;
              float fact = gauss(deltax,0.,conS);
	      float y = fact*hist->GetBinContent(i)*hist->GetBinWidth(i);
              con.SetBinContent(j,y+con.GetBinContent(j));
	    }



	  for (int i=1;i<=400;i++) 
	    {
	      // hist->SetBinContent(i,con.GetBinContent(i));
	     con.SetBinContent(i,0.);
	    }

          hist->Draw("same");



          func->SetParameter(0,0);
          func->SetParameter(1,1.);
	  //func->SetParameter(2,conS);
          func->SetParameter(2,0.1);
          func->SetParameter(3,1.);
          func->SetParameter(4,1.);
          func->SetParameter(5,1.);
          func->SetParameter(6,1.);

	  func->SetParLimits(3,0,100);
	  func->SetParLimits(4,0,100);
	  func->SetParLimits(5,0,100);
	  func->SetParLimits(6,0,100);


          func->SetLineColor(2);
          //func->Draw("same");



	  hist->Fit(func,"B");
          func->GetParameters(para);
          cout << "chisq=" << func->GetChisquare() << endl;
           if (fabs(para[1]-1.) < .2) 
	     { 

              slope *= para[1];
              intercept = intercept*para[1] + para[0];
	     }
            fout << itele << " " << istrip << " " 
                 << slope << " " << intercept << endl;
            fwhm << itele << " " << istrip << " " 
                  << para[2]*2.35 << endl;
            int ii = itele*32+istrip;
            xx[ii] = (float)ii;
            yy[ii] = para[2]*2.35;
            cout << para[0] << " " << para[1] << " " << para[2] << endl;

        }
    }

  TFile g("root/MixedSRing.root","RECREATE");
  for (int itele=0;itele<Ntele;itele++) canvas[itele]->Write();
  TCanvas fwhmCan("fwhm");
  TH2I frame2("frame2","",10,0,448,10,0,0.12);
  frame2.SetStats(kFALSE);
  frame2.GetYaxis()->SetTitle("FWHM [MeV]");
  frame2.GetXaxis()->SetTitle("front strip");
  frame2.Draw();
  TGraph graph(48,xx,yy);
  graph.Draw("*");
  graph.Write();
  fwhmCan.Write();
  g.Write();
}
