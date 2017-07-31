#include "hira.h"

 

hira::hira(TRandom * ran0, histo * Histo0)
{
  Histo = Histo0;
  ran = ran0;

  xmarker[0]=0x1ff3;
  xmarker[1]=0x1ff4;
  xmarker[2]=0x1ff5;


  //make map of chips
  ifstream ifile("chipmap.dat");

  if (!ifile.is_open())
    {
      cout << "chip map not found" << endl;
      abort();
    }

  string name;
  getline(ifile,name);
  int i1,iXLM,cbf,cbb,cbfl,cbbl;
  for (;;)
    {
      ifile >> i1 >> iXLM >> cbf>> cbb >> cbfl >> cbbl;
      if (ifile.eof()) break;
      if (ifile.bad()) break;

      iXLM -= 1;
      Map[iXLM][cbf].front = true;
      Map[iXLM][cbf].A = true;
      Map[iXLM][cbf].itele = i1;

      Map[iXLM][cbb].front = false;
      Map[iXLM][cbb].A = true;
      Map[iXLM][cbb].itele = i1;

      Map[iXLM][cbfl].front = true;
      Map[iXLM][cbfl].A = false;
      Map[iXLM][cbfl].itele = i1;

      Map[iXLM][cbbl].front = false;
      Map[iXLM][cbbl].A = false;
      Map[iXLM][cbbl].itele = i1;      
    }
  ifile.close();
  ifile.clear();


  //read in calibrations
  int Ntele = 14;
  int Nstrip = 32;
  name = "cal/front.cal";
  calFront = new calibrate(Ntele,Nstrip,name,1);
  name = "cal/back.cal";
  calBack = new calibrate(Ntele,Nstrip,name,1);

  name = "cal/HL_Front.cal";
  calLFront = new calibrate(2,Nstrip,name,1);
  name = "cal/HL_Back.cal";
  calLBack = new calibrate(2,Nstrip,name,1);

  name = "cal/Si_Back_dual.cal";
  calBDual = new calibrate(2,Nstrip,name,1);

  name = "cal/SiRecalib.cal";
  calrecal = new calibrate(2,Nstrip,name,1);


  name = "cal/Proton_New_ang.cal";
  calCsi = new calibrate(1,56,name,1);

  //TIME Calibrations

  name = "cal/timePulser_front.cal";
  calFrontT = new calibrate(Ntele,Nstrip,name,1);
  name = "cal/timePulser_frontLG.cal";
  calFrontTLG = new calibrate(Ntele,Nstrip,name,1);

  name = "cal/timePulser_back.cal";
  calBackT = new calibrate(Ntele,Nstrip,name,1);
  name = "cal/timePulser_backLG.cal";
  calBackTLG = new calibrate(Ntele,Nstrip,name,1);




  Telescope = new telescope*[14];
  for (int i=0;i<14;i++)
    {
      Telescope[i] = new telescope(ran,i,Histo);
      Telescope[i]->load(0, 0, 15,15,
			 16,16,31,31,
			 0, 15,15,0,
			 16,31,31,16);
    }


  //high low correlations zero arrays
  for (int i=0;i<2;i++)
    for (int j=0;j<32;j++)
      {
	fsumN[i][j] = 0.;
	fsumx[i][j] = 0.;
	fsumxx[i][j] = 0.;
        fsumyx[i][j] = 0.;
        fsumy[i][j] = 0.;
	
	bsumN[i][j] = 0.;
	bsumx[i][j] = 0.;
	bsumxx[i][j] = 0.;
        bsumyx[i][j] = 0.;
        bsumy[i][j] = 0.;
      }

}
//*************************************************************
bool hira::unpack(unsigned short *&point,int runno)
{
  bool stat = true;
  //  stat = unpackSi_sis(point);
   stat = unpackSi_adc(point);
    if (!stat) return stat;
    stat = unpackCsi(point,runno); 
}
//*************************************************************
bool hira::unpackCsi(unsigned short *&point,int runno)
{
  NE = 0;
  CsIM =0;

  // cout << "point in Csi " << hex << *point << endl; 
  
  /*   
  unsigned short * pig = point;
  for (int i=0;i<250;i++)
    {
      cout << dec<< *pig << " " << hex<< *pig << endl;
      pig++;
    }
  cout << endl;
  abort();
  //return true;
  */

  for (int iadc = 0;iadc<2;iadc++)
    {

      //check for ffff's
     unsigned short f3 = *point;
     unsigned short f4 = *(point+1);
     if (f3 == 0xffff && f4 == 0xffff)
       {
	 point += 2;
	 continue;
       } 

      ADC.number = 0;
      point = ADC.read(point);  // suck out the info in the qdc
      for (int i=0;i<ADC.number;i++)
	{


	  if (ADC.underflow[i]) continue;
	  if (ADC.overflow[i]) continue;
          int id = ADC.channel[i] + 32*iadc;
          int ienergy = ADC.data[i];
	  if(id < 56)
	    {

             float energy = calCsi->getEnergy(0,id,ienergy+ran->Rndm());
	     DataE[NE].id = id;
             DataE[NE].ienergy = ienergy;
	     DataE[NE].energy = energy;
	     Histo->ECsI[id]->Fill(ienergy);
	     Histo->ECsISum->Fill(id,ienergy);
	     Histo->ECsICSum->Fill(id,energy);


	      int itele = DataE[NE].id/4;
	      int icsi = DataE[NE].id%4; 
	      
              if (energy > 1 && id !=18)
		{
		     CsIM++;
		     Telescope[itele]->Csi.Add(icsi,0,energy,DataE[NE].ienergy,0.);
		}
	      //cout << Telescope[itele]->Csi.Order[0].strip;
	      //cout << " " << Telescope[itele]->Csi.Order[0].energy << endl;
	      
	     NE++;
	    }
	      
	  else if(id ==56 && runno >= 170 && runno <=172) //Blocker CsI
	    {
	      Histo->Blocker_E->Fill(ienergy);
	      Blocker_e = ienergy;
	    }


	}
      
      Histo->CsIMult->Fill(CsIM);

      //check for ffff's
     unsigned short f1 = *point;
     point++;
     unsigned short f2 = *point;
     point++;
     if(f1 != 0xffff && f2 != 0xffff) return false;
    }



  NT = 0;
  for (int itdc = 0;itdc<2;itdc++)
    {

      //check for ffff's
     unsigned short f3 = *point;
     unsigned short f4 = *(point+1);
     if (f3 == 0xffff && f4 == 0xffff) 
       {
	 point += 2;
        continue;
       } 

      TDC.number = 0;
      point = TDC.read(point);  // suck out the info in the qdc
      for (int i=0;i<TDC.number;i++)
	{


	  if (TDC.underflow[i]) continue;
	  if (TDC.overflow[i]) continue;
          int id = TDC.channel[i] + 32*itdc;
          int itime = TDC.data[i];
	  if (id < 56)
	    {
	     DataT[NT].id = id;
             DataT[NT].itime = itime;
	     Histo->TCsI[id]->Fill(itime);
	     NT++;
	    }
	}


      //check for ffff's
     unsigned short f1 = *point;
     point++;
     unsigned short f2 = *point;
     point++;
     if(f1 != 0xffff && f2 != 0xffff) return false;


    }


  /*
  // match up energies to times
  for (int ie=0;ie<NE;ie++)
    {
      DataE[ie].itime = -1;
      for (int it=0;it<NT;it++)
	{
          if (DataE[ie].id == DataT[it].id ) 	      //we have matched
	    {
	      DataE[ie].itime = DataT[it].itime;
	      int itele = DataE[ie].id/4;
	      int icsi = DataE[ie].id%4;
	      if(DataE[ie].energy >1.)// && DataE[ie].itime > 500 && DataE[ie].itime < 1500)
		Telescope[itele]->Csi.Add(icsi,0,DataE[ie].energy,DataE[ie].ienergy,DataE[ie].itime);

	    }
	  else if (DataE[ie].id < DataT[it].id) break; // no match found
	}
    }
    */  






  bool stat = true;
  return stat;
}
//***************************************************************
  //unpacking the XLM with internal ADC's
bool hira::unpackSi_adc(unsigned short *&point)
{
  
  /*
  unsigned short * pig = point;
  for (int i=0;i<800;i++)
    {
      cout << dec<< *pig << " " << hex<< *pig << endl;
      pig++;
    }
  abort();
  */




  unsigned short marker;
  unsigned short * nextMB = point;
  for (int iMB=0;iMB<3;iMB++)
    {
      point = nextMB;
      //    cout << iMB << " " << hex <<*point;
      marker = *point++;
      //if(marker !=xmarker[iMB])point++;
      if (marker != xmarker[iMB]) return false;
      
     int NstripsRead = 0;
     unsigned short chipWords = *point;
     nextMB = point + chipWords + 2;
     //if (chipWords > 400) return false;  // please fix Kyle
     if (chipWords == 0)
        {
         NstripsRead = 0;
         return (bool) 1;
        }
     point += 2;
     NstripsRead = *point;

     if (NstripsRead > 384) return false; // bad buffer
     point += 5;

     //  cout << chipWords << " " << NstripsRead*3+7 <<  " " << NstripsRead <<  endl;
     for (int istrip = 0;istrip < NstripsRead;istrip++)
        {

          unsigned short id = *point++;
          unsigned short chipNum = (id&0x1FE0)>>5;
          unsigned short chanNum = id & 0x1F;
          unsigned short ienergy = *point++;
          unsigned short itime =  *point++;



          unsigned short underOver = 0;
          if (ienergy == 32767) underOver = -1; // overflow
          if (ienergy == 16384) underOver = 1;// underflow

	  ienergy  = 16383 - ienergy;
	  //by chip board
	  bool secondChip = false;
          if (chipNum%2 == 0)
	    {
	      if (iMB ==2) secondChip = true;
	      // else chanNum += 16;
	      else chanNum = 31- 2*chanNum-1;
	      chipNum /= 2;
	    }
	  else
	    {
	      chipNum = chipNum/2 + 1;
	      if(iMB !=2) chanNum = 31 - 2*chanNum;
	    }
	  
	  //  cout << id << " " << chipNum << " " << chanNum << " " << ienergy << " " << itime << " " << iMB << endl;
          bool bfront = Map[iMB][chipNum].front;
          bool bA = Map[iMB][chipNum].A;
          int itele = Map[iMB][chipNum].itele - 1;



	  //if (iMB == 0) cout << itele << " " << chipNum << endl;

	  if(itele <0)	 return false;
          bool bhigh = true; 

	  if (iMB == 2)
	    {
              if (secondChip) bhigh = false;
              else bhigh = true;

              if (!bA) chanNum += 16;

	    }

	  if(bfront) chanNum = 31- chanNum;


	  if (chanNum > 31)
	    {
	      cout << "chanNum too big" << endl;
              return false;
	    }
          if (chipNum > 12)
	    {
	      cout << "chipNum too big " << chipNum << endl;
              return false;
	    }


	  if (bhigh && bfront)
	    {

	      float energy = calFront->getEnergy(itele,chanNum,ienergy+ran->Rndm());
	      //if(itele ==6 || itele ==7) energy = energy*1.05; //Second Run!!!!!!


	      float time = calFrontT->getEnergy(itele,chanNum,itime+ran->Rndm());
	      Histo->EFTSum[itele]->Fill(chanNum,time);
	      //	      Histo->SiFTime->Fill(time);

	      //Recalibrating
	      if(itele ==6 || itele ==7)
		energy = calrecal->getEnergy(itele-6,chanNum,energy);

	      Histo->EfrontR[itele][chanNum]->Fill(ienergy);
	      Histo->TfrontR[itele][chanNum]->Fill(time);
	      Histo->EfrontC[itele][chanNum]->Fill(energy);
	      //Histo->EFSum[itele]->Fill(chanNum,ienergy);
	      // Histo->EFCSum[itele]->Fill(chanNum,energy);

              if(energy >0.75)
		{
		  Telescope[itele]->Front.Add(chanNum,underOver,energy,ienergy,time);
		}
		
	    }
	  if (bhigh && !bfront)
	    {
	      float energy = calBack->getEnergy(itele,chanNum,ienergy+ran->Rndm());
	      float time = calBackT->getEnergy(itele,chanNum,itime+ran->Rndm());
	      Histo->EBTSum[itele]->Fill(chanNum,time);
	      Histo->SiBTime->Fill(time);




	      Histo->EbackR[itele][chanNum]->Fill(ienergy);
	      Histo->TbackR[itele][chanNum]->Fill(time);
	      Histo->EbackC[itele][chanNum]->Fill(energy);
	      //Histo->EBSum[itele]->Fill(chanNum,ienergy);
	      // Histo->EBCSum[itele]->Fill(chanNum,energy);

	      if(energy > 0.75)
		Telescope[itele]->Back.Add(chanNum,underOver,energy,ienergy,time);

	    }

	  if(bhigh && bfront)
	    {
	      for (int i=0;i<Telescope[itele]->Front.Nstore;i++)
		{
                  if (Telescope[itele]->Front.Order[i].strip == chanNum)
		    Telescope[itele]->Front.Order[i].energylow = 0.;

		}
	    }

	  if(bhigh && !bfront)
	    {
	      for (int i=0;i<Telescope[itele]->Back.Nstore;i++)
		{
                  if (Telescope[itele]->Back.Order[i].strip == chanNum)
		    Telescope[itele]->Back.Order[i].energylow = 0.;
		}
	    }


	  if (!bhigh && bfront)
	    {

	      float time = calFrontTLG->getEnergy(itele,chanNum,itime+ran->Rndm());
	      Histo->SiFTime->Fill(time);

	      Histo->EfrontLGR[itele][chanNum]->Fill(ienergy);
	      Histo->TfrontLG[itele][chanNum]->Fill(time);
	      //Histo->EFLSum[itele]->Fill(chanNum,ienergy);

	      float energy = calLFront->getEnergy(itele-6,chanNum,ienergy+ran->Rndm());
	      energy = calFront->getEnergy(itele,chanNum,energy+ran->Rndm());
	      Histo->EfrontLGC[itele][chanNum]->Fill(energy);

	      //Recalibrating
	      if(itele ==6 || itele ==7)
		energy = calrecal->getEnergy(itele-6,chanNum,energy);


	      for (int i=0;i<Telescope[itele]->Front.Nstore;i++)
		{
                  if (Telescope[itele]->Front.Order[i].strip == chanNum)
		    {
		      Telescope[itele]->Front.Order[i].energyRlow = ienergy;
		      Telescope[itele]->Front.Order[i].energylow = energy;
                      double ienergyHigh=Telescope[itele]->Front.Order[i].energyR;
                      if (ienergyHigh < 15000 && ienergy > 50 && ienergy <15000)
			{
			  float Ratio = 0.;
                          fsumN[itele-6][chanNum]++;
                          fsumx[itele-6][chanNum] += (double)ienergy;
                          fsumxx[itele-6][chanNum] += pow((double)ienergy,2);
                          fsumy[itele-6][chanNum] += ienergyHigh;
                          fsumyx[itele-6][chanNum] += ienergyHigh*(double)ienergy;
			}
                      
		      break;
		    }
		}
	      
	    }
	  
	  
	  if (!bhigh && !bfront)
	    {

	      float time = calBackTLG->getEnergy(itele,chanNum,itime+ran->Rndm());
	      Histo->SiBTime->Fill(time);
	      Histo->EbackLGR[itele][chanNum]->Fill(ienergy);
	      Histo->TbackLG[itele][chanNum]->Fill(time);
              float energy = ienergy;
	      energy = calBDual->getEnergy(itele-6,chanNum,ienergy+ran->Rndm());
	      if(ienergy > energy)energy = ienergy;
	      float corenergy = energy;
	      Histo->EbackLGCC[itele][chanNum]->Fill(energy);
	      energy = calLBack->getEnergy(itele-6,chanNum,energy+ran->Rndm());
	      energy = calBack->getEnergy(itele,chanNum,energy+ran->Rndm());
	      Histo->EbackLGC[itele][chanNum]->Fill(energy);
	      //Histo->EBLSum[itele]->Fill(chanNum,ienergy);

	      for (int i=0;i<Telescope[itele]->Back.Nstore;i++)
		{
		  if (Telescope[itele]->Back.Order[i].strip == chanNum)
		    {

                      Telescope[itele]->Back.Order[i].energyRlow = ienergy;
                      Telescope[itele]->Back.Order[i].energylow = energy;
                      double ienergyHigh=Telescope[itele]->Back.Order[i].energyR;
                      if (ienergyHigh < 15000 && ienergy > 50 && ienergy <15000)
			{
                          bsumN[itele-6][chanNum]++;
                          bsumx[itele-6][chanNum] += (double)corenergy;
                          bsumxx[itele-6][chanNum] += pow((double)corenergy,2);
                          bsumy[itele-6][chanNum] += ienergyHigh;
                          bsumyx[itele-6][chanNum] += ienergyHigh*(double)corenergy;
			}
                      
		      break;
		    }
		}
	      
	    }          
	  
	}

    }
  point = nextMB;
  if (*point == 0xe0fe) return false;
  //  cout << "point in si " << hex<< *point << " " << dec << *point << endl;
  return true;
}


//***************************************************************
//** unpacking with XLM and sis modules
bool hira::unpackSi_sis(unsigned short *&point)
{

  unsigned short marker;
  for (int iXLM=0;iXLM<3;iXLM++)
    {
      marker = *point++;
      unsigned short *epoint = point;
      if (marker == xmarker[iXLM]) //XLM1
        {

	  // in princle words and words2 are both part of 
	  // a int32 parameter, but in practicable terms
	  // words2 should always be zero
          unsigned short words = *epoint++;   // 
          unsigned short words2 = *epoint++;
	  if (words2 != 0) 
	    {
              cout << "words2= " << words2 << endl;
              return false;
	    }
          point = epoint + words;
          unsigned short channelCount = *epoint++;  //number of data in XLM
	  //if (channelCount > 0) cout << iXLM << "  " << channelCount << endl;
          epoint += 4;
	  nXLM[iXLM] = channelCount;

	  for (int idata = 0;idata<channelCount;idata++)
	    {
              chanXLM[iXLM][idata] = *epoint++;
              if (chanXLM[iXLM][idata] == 0) return false;
	    }
       
        }
      else cout << "XLM" << iXLM << "  missing" << endl; 
    }
  int itry = 0;
  for(;;)
    {
      if (itry >10) return false;
      marker = *point++;
      if (marker == 0xfadc) break;
      itry++;
    }

  point++;
  for (int iXLM = 0; iXLM<3;iXLM++)
    {
      unsigned short words = *point++;
      unsigned short words2 = *point++;
      
      //cout << iXLM << " " << words << endl;

      if (words2 != 0) 
	{
         cout << "Words in FADC .ne. 0" << endl;
	 return false;
	}

      if (words != nXLM[iXLM])
        {
          cout << " XLM and FDC channels do not match for XLM " << iXLM <<endl;
          return false;
	}
 

      for (int i=0;i<words;i++)
	{
         
          unsigned short ienergy = *point++;
          unsigned short time = *point++;
          unsigned short chan = chanXLM[iXLM][i];

	  //by chip
          unsigned short chipNum =  (chan >> 5) & 0xff;
          unsigned short chanNum = chan & 0x0F;
	  unsigned short underOver = 0;
          //cout << "chipNum = " << chipNum << " chanNum= " << chanNum << endl;


	  //by chip board
	  bool secondChip = false;
          if (chipNum%2 == 0)
	    {
	      if (iXLM ==2) secondChip = true;
	      else chanNum += 16;
	     chipNum /= 2;
	    }
	  else
	    {
	     chipNum = chipNum/2 + 1;
	    }


	  bool bfront = Map[iXLM][chipNum].front;
          bool bA = Map[iXLM][chipNum].A;
          int itele = Map[iXLM][chipNum].itele - 1;
	  if(itele <0)	 return false;
          bool bhigh = true; 

	  if (iXLM == 2)
	    {
              if (secondChip) bhigh = false;
              else bhigh = true;

              if (!bA) chanNum += 16;

	    }

          /*
          cout << i << " e= " << ienergy << " t= " << time << " chipNum= " << 
	    chipNum << " chanNum= " << chanNum << " bfront= " 
	       << bfront << " bA=" << bA << " itele= " << itele << " bhigh= "
	       << bhigh << endl; 
	  */
	  if (bhigh && bfront)
	    {
	      float energy = calFront->getEnergy(itele,chanNum,ienergy+ran->Rndm());
	      Histo->EfrontR[itele][chanNum]->Fill(ienergy);
	      Histo->TfrontR[itele][chanNum]->Fill(time);
	      Histo->EfrontC[itele][chanNum]->Fill(energy);


              Telescope[itele]->Front.Add(chanNum,underOver,energy,ienergy,time);

	    }
	  if (bhigh && !bfront)
	    {
              float energy = calBack->getEnergy(itele,chanNum,ienergy+ran->Rndm());
	      Histo->EbackR[itele][chanNum]->Fill(ienergy);
	      Histo->TbackR[itele][chanNum]->Fill(time);
	      Histo->EbackC[itele][chanNum]->Fill(energy);

              Telescope[itele]->Back.Add(chanNum,underOver,energy,ienergy,time);

	    }

	  if (!bhigh && bfront)
	    {
	      Histo->EfrontLGR[itele][chanNum]->Fill(ienergy);
	      Histo->TfrontLG[itele][chanNum]->Fill(time);


	      for (int i=0;i<Telescope[itele]->Front.Nstore;i++)
		{
                  if (Telescope[itele]->Front.Order[i].strip == chanNum)
		    {
                      Telescope[itele]->Front.Order[i].energyRlow = ienergy;
                      double ienergyHigh=Telescope[itele]->Front.Order[i].energyR;
                      if (ienergyHigh < 15000 && ienergy > 50 && ienergy <15000)
			{
			  float Ratio = 0.;
                          fsumN[itele-6][chanNum]++;
                          fsumx[itele-6][chanNum] += (double)ienergy;
                          fsumxx[itele-6][chanNum] += pow((double)ienergy,2);
                          fsumy[itele-6][chanNum] += ienergyHigh;
                          fsumyx[itele-6][chanNum] += ienergyHigh*(double)ienergy;
			}
                      
		      break;
		    }
		}


	    }          

	  if (!bhigh && !bfront)
	    {
	      Histo->EbackLGR[itele][chanNum]->Fill(ienergy);
	      Histo->TbackLG[itele][chanNum]->Fill(time);

	      for (int i=0;i<Telescope[itele]->Back.Nstore;i++)
		{
                  if (Telescope[itele]->Back.Order[i].strip == chanNum)
		    {
                      Telescope[itele]->Back.Order[i].energyRlow = ienergy;
                      double ienergyHigh=Telescope[itele]->Back.Order[i].energyR;
                      if (ienergyHigh < 15000 && ienergy > 50 && ienergy <15000)
			{
                          bsumN[itele-6][chanNum]++;
                          bsumx[itele-6][chanNum] += (double)ienergy;
                          bsumxx[itele-6][chanNum] += pow((double)ienergy,2);
                          bsumy[itele-6][chanNum] += ienergyHigh;
                          bsumyx[itele-6][chanNum] += ienergyHigh*(double)ienergy;
			}
                      
		      break;
		    }
		}

	    }          

	  if (chanNum > 31)
	    {
	      cout << "chanNum too big" << endl;
              return false;
	    }
          if (chipNum > 11)
	    {
	      cout << "chipNum too big" << endl;
              return false;
	    }
	}

      marker = *point++;
      if (marker != 0xaaaa) cout << "marker = " << marker << endl;
    }
  



  return true;
}
//***************************************************************
hira::~hira()
{

  cout << "start Hira destr" << endl;
  //high-low correlations

  /*   
  for (int i=0;i<2;i++)
    for (int j=0;j<32;j++)
      {
        double delta = fsumN[i][j]*fsumxx[i][j]- pow(fsumx[i][j],2);
        if (delta == 0) continue;
        double slope = fsumN[i][j]*fsumyx[i][j] - fsumx[i][j]*fsumy[i][j];
        slope /= delta;
        double intercept = (fsumy[i][j] - slope*fsumx[i][j])/fsumN[i][j];
	cout << i << " " << j << " " << slope << " " << intercept << endl;

      }

  for (int i=0;i<2;i++)
    for (int j=0;j<32;j++)
      {
        double delta = bsumN[i][j]*bsumxx[i][j]- pow(bsumx[i][j],2);
        if (delta == 0) continue;
        double slope = bsumN[i][j]*bsumyx[i][j] - bsumx[i][j]*bsumy[i][j];
        slope /= delta;
        double intercept = (bsumy[i][j] - slope*bsumx[i][j])/bsumN[i][j];
	cout << i << " " << j << " " << slope << " " << intercept << endl;

      }
  */  



  delete calFront;
  delete calBack;
  delete calCsi;
  delete calLBack;
  delete calLFront;
  delete calBDual;
  delete calFrontT;
  delete calBackT;
  delete calFrontTLG;
  delete calBackTLG;


  for(int i=0;i<14;i++)
    {
      delete Telescope[i];
    }
  delete [] Telescope;

  cout << "stop Hira destr" << endl;
}

//**********************************************
void hira::reset()
{
  fred = false;
  for (int i=0;i<14;i++) Telescope[i]->reset();
}
