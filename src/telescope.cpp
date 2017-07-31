#include "telescope.h"

bool const telescope::relativity=1;

telescope::telescope(TRandom * ran0, int id0, histo * Histo0)
{
  //  cout << "in telescope constructor" << endl;
  id = id0;
  Histo = Histo0;
  ran = ran0;

  fitCal =0; // set this to 1 to use polynomial fit of block CsI calibrations (only set up for Russian atm)




  Csicheck=8;
  triton_shift = 0.;//-1.5;//-2.2;//-0.8;//0.;//0.6; //remember you subtract this from the energy below
  triton_par1 = 1.;//0.968;//1.; //.98;
  alpha_quad=0.0000002;//000001;
  triton_quad=0.00000025;
  // cout << "alpha_quad = " << alpha_quad << endl;


  ostringstream outstring;
  string name;
  string name1;
  string name2;
  int istart = -1;
  int iend = -1;
  int blocksS2=8;
  int blocksRus=8;

  if(id ==0)
    {
      istart = 0;
      iend = 16;
    }
  else if(id ==1)
    {
      istart =16;
      iend = 32;
    }

  for(int i = istart;i<iend;i++)
    {
      outstring.str("");
      outstring << "pid"<<i;
      name = outstring.str();
      Pid[i] = new pid(name); 
    }

  if(id==1)
    {
      for(int i = istart;i<iend;i++)
	{
	  outstring.str("");
	  outstring << "blockCalCsid"<<i;
	  name = outstring.str();
	  blockCalCsid[i] = new blockCal(name,blocksS2); 
	}
      
      for(int i = istart;i<iend;i++)
	{
	  outstring.str("");
	  outstring << "blockCalCsiA"<<i;
	  name = outstring.str();
	  blockCalCsiA[i] = new blockCal(name,blocksS2); 
	}
    }

  if(id==0)
    {
      for(int i = istart;i<iend;i++)
	{
	  outstring.str("");
	  outstring << "blockCalCsid"<<i;
	  name = outstring.str();
	  blockCalCsid[i] = new blockCal(name,blocksRus); 
	}
      
      for(int i = istart;i<iend;i++)
	{
	  outstring.str("");
	  outstring << "blockCalCsiA"<<i;
	  name = outstring.str();
	  blockCalCsiA[i] = new blockCal(name,blocksRus); 
	}
    }





  //edit the loss files to make target dependent (want to load more loss classes then deal with target switching when loss is calculated)

  name="Hydrogen.loss"; //in 9Be
  Loss[0] = new CLoss(name,1.);
  Loss[1] = new CLoss(name,2);
  Loss[2] = new CLoss(name,3.);
  name="Helium.loss"; //in 9 Be
  Loss[3] = new CLoss(name,3.);
  Loss[4] = new CLoss(name,4);
  name = "Lithium.loss";// in 9 Be
  Loss[5] = new CLoss(name,6.);
  Loss[6] = new CLoss(name,7.);
  Loss[7] = new CLoss(name,8.);

  name="Hydrogen_12C.loss"; //in C
  Loss[8] = new CLoss(name,1.);
  Loss[9] = new CLoss(name,2.);
  Loss[10] = new CLoss(name,3.);
  name="Helium_12C.loss"; //in C
  Loss[11] = new CLoss(name,3.);
  Loss[12] = new CLoss(name,4.);
  name = "Lithium_12C.loss";// in C
  Loss[13] = new CLoss(name,6.);
  Loss[14] = new CLoss(name,7.);
  Loss[15] = new CLoss(name,8.);


  name="Hydrogen_27Al.loss"; //in Al
  Loss[16] = new CLoss(name,1.);
  Loss[17] = new CLoss(name,2.);
  Loss[18] = new CLoss(name,3.);
  name="Helium_27Al.loss"; //in Al
  Loss[19] = new CLoss(name,3.);
  Loss[20] = new CLoss(name,4.);
  name = "Lithium_27Al.loss";// in Al
  Loss[21] = new CLoss(name,6.);
  Loss[22] = new CLoss(name,7.);
  Loss[23] = new CLoss(name,8.);

  name = "Hydrogen_Ta.loss";
  Loss[24] = new CLoss(name,1.);
  Loss[25] = new CLoss(name,2);
  Loss[26] = new CLoss(name,3.);

  name = "Helium_Ta.loss";
  Loss[27] = new CLoss(name,3.);
  Loss[28] = new CLoss(name,4.003);

  name = "Hydrogen_Si.loss";
  Loss[29] = new CLoss(name,1.);
  Loss[30] = new CLoss(name,2);
  Loss[31] = new CLoss(name,3.);
  
  name = "Helium_Si.loss";
  Loss[32] = new CLoss(name,3.);
  Loss[33] = new CLoss(name,4.003);




  //recalibrate depending on isotope the CsI (proton effect energy (MeV)-> Energy (MeV))


  name = "cal/7LiTAMU-H.cal"; //proton
  calCsi = new calibrate(1,32,name,1);
  name = "cal/7LiTAMU-d-CsI.cal";//deuteron
  calCsid = new calibrate(1,32,name,1);
  name = "cal/7LiTAMU-H.cal";//triton
  calCsit = new calibrate(1,32,name,1);

  name = "cal/7LiTAMU-He.cal";//3He
  calCsiHe3 = new calibrate(1,32,name,1);
  name = "cal/7LiTAMU-He.cal";//4He
  calCsiA = new calibrate(1,32,name,1);


  name = "cal/Lithium6_New_ang.cal";//6Li
  calCsiLi6 = new calibrate(1,56,name,1);



  
  name1 = "cal/posCorrect_slope_d.cal";
  name2 = "cal/posCorrect_intercept_d.cal";
  blockCalCsidRus = new blockCal(name1, name2);

  name1 = "cal/posCorrect_slope_a.cal";
  name2 = "cal/posCorrect_intercept_a.cal";
  blockCalCsiARus = new blockCal(name1, name2);
  



  //file to check zlines
  
  ifstream check("zline_check.txt");
  if(!check.is_open()){cout << "Unable to open zline check file" << endl;}
  for(int i=0;i<32;i++)
    {
      check >> check_slope[i] >> check_intercept[i];
      // cout << check_slope[i] <<"\t" << check_intercept[i] << endl;
    }

  /*
  for(int i=0;i<48;i++)
    {
      bins[i] = Pixels.TeleP[1].Location[0][i].theta-0.5;
    }
  bins[48] =Pixels.TeleP[0].Location[0][47].theta*180./;
  for(int i=0;i<32;i++)
    {
      bins[i+49] =Pixels.TeleP[0].Location[0][i].theta-0.859375;
    }
  bins[81] =Pixels.TeleP[0].Location[0][31].theta+0.859375;

  Histo->Li7_AbsElasXS->SetBins(82,bins);
  */


}
//************************************************
telescope::~telescope()
{
  for(int ii =0;ii<34;ii++)
    {
      delete Loss[ii];
    }
  for(int ii=0;ii<32;ii++)
    {
      delete Pid[ii];
      delete blockCalCsid[ii];
      delete blockCalCsiA[ii];
    }
}
//***********************************************

int telescope::Block(int ring)
{


  if(id==0)
    {
      if(ring>=0&&ring<=3){return 1;}//cout << "really?" << endl;}
      if(ring>=4&&ring<=7){return 2;}//cout << "really?" << endl;}
      if(ring>=8&&ring<=11){return 3;}//cout << "really?" << endl;}
      if(ring>=12&&ring<=15){return 4;}//cout << "really?" << endl;}
      if(ring>=16&&ring<=19){return 5;}//cout << "really?" << endl;}
      if(ring>=20&&ring<=23){return 6;}//cout << "really?" << endl;}
      if(ring>=24&&ring<=27){return 7;}//cout << "really?" << endl;}
      if(ring>=28&&ring<=31){return 8;}
    }

  if(id==1)
    {
      if(ring>=0&&ring<=5){return 1;}//cout << "really?" << endl;}
      if(ring>=6&&ring<=11){return 2;}//cout << "really?" << endl;}
      if(ring>=12&&ring<=17){return 3;}//cout << "really?" << endl;}
      if(ring>=18&&ring<=23){return 4;}//cout << "really?" << endl;}
      if(ring>=24&&ring<=29){return 5;}//cout << "really?" << endl;}
      if(ring>=30&&ring<=35){return 6;}//cout << "really?" << endl;}
      if(ring>=36&&ring<=41){return 7;}//cout << "really?" << endl;}
      if(ring>=42&&ring<=47){return 7;} // not much data for calibration here, might consider throwing these rings out of analysis
    }


  return 0;


}



void telescope::analyze(int event)
{
  //  cout << "in telescope::analyze" << endl;
  Nsolution = 0;
  Np = 0;
  N6 = 0;
  penergy = 0;
  renergy = 0;

  Event = event;


  if(Pie.Nstore ==0 || Ring.Nstore == 0) return;
  
  phit = Pie.Order[0].strip;
  rhit = Ring.Order[0].strip;
  CsIhit = Csi.Order[0].strip + 16*id;

  penergy = Pie.Order[0].energy;
  renergy = Ring.Order[0].energy;

  if(id ==0)
    Histo->RusRingvPie->Fill(penergy,renergy);
  if(id ==1)
    Histo->S2RingvPie->Fill(penergy,renergy);

  if(id==1){if(rhit>=44)return;if(phit==15||phit==0)return;}

  /* //gets rid of all the pies on the outeredge of the CsI in the Russian
  if(id==0)
    {
      if(phit==1||phit==2||phit==5||phit==6||phit==9||phit==10||phit==13||phit==14||phit==17||phit==18)return;
      if(phit==21||phit==22||phit==25||phit==26||phit==29||phit==30||phit==33||phit==34||phit==37||phit==38)return;
      if(phit==41||phit==42||phit==45||phit==46||phit==49||phit==50||phit==53||phit==54||phit==57||phit==58)return;
      if(phit==61||phit==62)return;
    }
  */

  theta = Pixels.getAngle(id,phit,rhit);
  //theta = theta + Pixels.deltatheta*(ran->Rndm()-0.5);
  theta_use = theta*180./3.1415927;
  phi = Pixels.phi;
  phi = phi+Pixels.deltaphi*(ran->Rndm()-0.5); //correcting for pixels at center of sections, should put this somewhere else eventually
  
  // cout << phi << " " << theta << endl;
  //cout << Pixels.deltaphi << " " << Pixels.deltatheta << endl;
  xhit = Pixels.r*cos(phi)*sin(theta);
  yhit = Pixels.r*sin(phi)*sin(theta);


  //if(id==0 && yhit < 0.5 && yhit > -0.5 && xhit < 15 && xhit > -15) cout << xhit << endl;

  double z = Pixels.r*cos(theta);
  float rperp = sqrt(pow(xhit,2.)+pow(yhit,2.));

 

  if(rperp<15&&id==0){cout << "bzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz" << endl;}//getting rid of "dead" layer

  if(id==0&&rperp > 15.)
    {
      Histo->HitMap->Fill(-xhit,yhit);
      //if(rhit > 30)Histo->RusPhi_gated->Fill(phi*180./3.1415927);
    }
  if(id==1)
    {
      Histo->HitMap1->Fill(-xhit,yhit);
      Histo->RusPhi_gated->Fill(phi*180./3.1415927);
    }
  xhit = theta*180/3.1415927*cos(phi);
  yhit = theta*180/3.1415927*sin(phi);

  //Histo->HitMap_perspective->Fill(-xhit,yhit);

   //cout << rperp << endl;

  // cout << "Nrings = " << Ring.Nstore << endl;




  bool matched =0;
  

  /*
  if(Csi.Nstore!=1&&CsIhit==Csicheck)
    {
      for(int j =1;j<Csi.Nstore;j++)
	{
	  if((Csi.Order[j].strip+16*id)==(Csicheck+1)&&Block(Ring.Order[0].strip)==1)Histo->D2CsI->Fill(Csi.Order[0].energyR,Csi.Order[j].energyR);
	}

    }
  */

  if(Csi.Nstore ==1)
    {
      
      int icsi = Csi.Order[0].strip;
      int ipie = Pie.Order[0].strip;
      int iring = Ring.Order[0].strip;


      if(CsIMap[ipie]==icsi)
	{
	  matched = 1;
	}	
      if(!matched)return;

  
      float help;
      float help2;

      CsIenergy = Csi.Order[0].energyR;
      Sienergy = penergy;//renergy;//penergy;


      Histo->dEE[CsIhit]->Fill(CsIenergy,Sienergy);
      
      float energy = Csi.Order[0].energy;
   
      block = Block(iring);  
      if(block==0)cout << "block error" << endl;
      //remember this currently only goes over russian
      if(CsIhit==Csicheck)
	{
	  Histo->dEEtest->Fill(CsIenergy,Sienergy);
	  if(iring==1||iring==4)Histo->TSi->Fill(CsIenergy,Csi.Order[0].time);
	  if(block==1)//&&Pie.Order[0].time>4150&&Pie.Order[0].time<4250)
	    {Histo->ECsIblock_raw1->Fill(CsIenergy);}
	  if(block==2)
	    {Histo->ECsIblock_raw2->Fill(CsIenergy);}
	  if(block==3)
	    {Histo->ECsIblock_raw3->Fill(CsIenergy);}
	  if(block==4)
	    {Histo->ECsIblock_raw4->Fill(CsIenergy);}
	  if(block==5)
	    {Histo->ECsIblock_raw5->Fill(CsIenergy);}
	  if(block==6)
	    {Histo->ECsIblock_raw6->Fill(CsIenergy);}
	  if(block==7)
	    {Histo->ECsIblock_raw7->Fill(CsIenergy);}
	  if(block==8)
	    {Histo->ECsIblock_raw8->Fill(CsIenergy);}
	  // if(id==1)
	  Histo->CsIRings[iring]->Fill(CsIenergy);
	  Histo->ECsI_raw_all->Fill(CsIenergy);
	}




      
      //use this for testing calibrations
      if(CsIhit==Csicheck)
	{

	  //cout << phit << endl;
	  //if(block>4)  
	  energy = blockCalCsidRus->getE(CsIenergy+ran->Rndm(), CsIhit, theta_use);
	  //cout << "Calslope new = " << blockCalCsiARus->calSlope << "\t" << "CalInt new = " << blockCalCsiARus->calInt << endl;
	   
	  
	  //if(block<=4)
	  // energy = blockCalCsid[CsIhit]->getE(CsIenergy+ran->Rndm(), block);
	  //cout << "Calslope old = " << blockCalCsiA[CsIhit]->calSlope << "\t" << "CalInt old = " << blockCalCsiA[CsIhit]->calInt << endl;
	    
	  //if(block==1)energy = 0.04496*(CsIenergy+ran->Rndm())-5.05200;

	  if(block==1)Histo->ECsIblock1->Fill(energy);
	  if(block==2){Histo->ECsIblock2->Fill(energy);}
	  if(block==3){Histo->ECsIblock3->Fill(energy);}
	  if(block==4){Histo->ECsIblock4->Fill(energy);}
	  if(block==5){Histo->ECsIblock5->Fill(energy);}
	  if(block==6){Histo->ECsIblock6->Fill(energy);}
	  if(block==7){Histo->ECsIblock7->Fill(energy);}
	  if(block==8){Histo->ECsIblock8->Fill(energy);}

	  //  if(phit==20)
	  // {
	  if(block==1){Histo->ESiblock1->Fill(Sienergy);}
	  if(block==2){Histo->ESiblock2->Fill(Sienergy);}
	  if(block==3){Histo->ESiblock3->Fill(Sienergy);}
	  if(block==4){Histo->ESiblock4->Fill(Sienergy);}
	  if(block==5){Histo->ESiblock5->Fill(Sienergy);}
	  if(block==6){Histo->ESiblock6->Fill(Sienergy);}
	  if(block==7){Histo->ESiblock7->Fill(Sienergy);}
	  if(block==8){Histo->ESiblock8->Fill(Sienergy);}
	  // }
	}
      

      
      bool stat = Pid[CsIhit]->getPID(CsIenergy,Sienergy);
      if (!stat) return;
      
      int Z = 0;
      Z= Pid[CsIhit]->Z;
      int A = 0;
      A = Pid[CsIhit]->A;

      if(A==7)Histo->HitMap_perspective->Fill(-xhit,yhit);
      
      if (Z > 0 && A >0)
	{
	  
	  if(Z ==1)
	    { 
	      if(A ==2)
		{
		  energy = calCsid->getEnergy(0,CsIhit,CsIenergy+ran->Rndm());
		  energy = blockCalCsid[CsIhit]->getE(CsIenergy+ran->Rndm(),block);


		  if(id==0&&fitCal==1)
		  {
		   energy = blockCalCsidRus->getE(CsIenergy+ran->Rndm(),CsIhit,theta_use);
		  }
		}
	      //CsIenergy= raw channel number
	      if(A ==3)
		{
		  energy = calCsit->getEnergy(0,CsIhit,CsIenergy+ran->Rndm()); //fix this
		  //Histo->dEE_ta[CsIhit]->Fill(CsIenergy,Sienergy); //will want to get rid of this eventually		  
		  energy=triton_par1*blockCalCsid[CsIhit]->getE(CsIenergy+ran->Rndm(),block)-triton_shift+triton_quad*pow(CsIenergy,2.);
		  if(id==0&&fitCal==1)
		  {
		  energy=triton_par1*blockCalCsidRus->getE(CsIenergy+ran->Rndm(),CsIhit,theta_use)-triton_shift;
		  }

		}
	    }
	  else if(Z == 2)
	    {
	      if(A ==3)
		energy = calCsiHe3->getEnergy(0,CsIhit,CsIenergy+ran->Rndm());
	      else if(A ==4)
		{
		  energy = calCsiA->getEnergy(0,CsIhit,CsIenergy+ran->Rndm()); //fix this
		  // Histo->dEE_ta[CsIhit]->Fill(CsIenergy,Sienergy); //will want to get rid of this eventually (slows stuff down)
		  // energy=blockCalCsiA[CsIhit]->getE(CsIenergy+ran->Rndm(),block);
		  energy=blockCalCsiA[CsIhit]->getE(CsIenergy+ran->Rndm(),block)+alpha_quad*pow(CsIenergy,2.);		
		  if(id==0&&fitCal==1)
		  {
		    energy=blockCalCsiARus->getE(CsIenergy+ran->Rndm(),CsIhit,theta_use);
		  }
		}
		  else 
		{
		  cout << "found no calib for" << Z << " " << A << endl; 
		  abort();
		}
	    }
	  else if (Z == 3) 
	    {
	      energy = calCsiLi6->getEnergy(0,CsIhit,CsIenergy+ran->Rndm());
	    }

	  float sumEnergy = energy + Sienergy;


	  //energy loss in target
	  int ipid =0;
	  ipid = A;

	  /*
	  if(Z ==2)
	    ipid = A;
	  else 
	    ipid = A-1;
	  */

	  float light = Csi.Order[0].energy;
	  float thick;
	  float Ein=0;
	  float Ekin=0;
	  int please;


	  //can add back in Tantalum loss here with appropriate loss file

	  float thickTa = 13.8/cos(theta);
    	  float deadlayer = 0./cos(theta); // mg/cm^2 11.606 = 50 um dead layer of Si
	  // float thickSi = 121.863/cos(theta);

	  
	  if(id==0)
	    {
	      float thickSi = 121.863/cos(theta);
	      if(Z==1)
		{
		  sumEnergy = Loss[28+ipid]->getEin(energy,thickSi);
		  //		  sumEnergy = energy + Sienergy;
		}
	      if(Z==2)
		{
		  sumEnergy = Loss[29+ipid]->getEin(energy,thickSi);
		  //sumEnergy = energy + Sienergy;
		}
	    }

	  /*
	  if(id==1)
	    {
	      float thickSi = 121.863/cos(theta)
	      if(Z==1)
		{
		  sumEnergy = Loss[28+ipid]->getEin(energy,thickSi);
		  //		  sumEnergy = energy + Sienergy;
		}
	      if(Z==2)
		{
		  sumEnergy = Loss[29+ipid]->getEin(energy,thickSi);
		  //sumEnergy = energy + Sienergy;
		}
	    }
	  */
	  

	  if(CsIhit==Csicheck)
	    {
	
	      Histo->ECsI_sumenergy->Fill(sumEnergy);
	      // cout << sumEnergy/A << endl;
	      if(Z==1)
	      	{
		  // if(id==0)Ein = Loss[30]->getEin(sumEnergy,deadlayer); //change these depending on calibration beams
		  // cout << energy << endl;
		  Ein=sumEnergy;// energy = blockCalCsid[CsIhit]->getEnergy(0,CsIhit,CsIenergy+ran->Rndm());
		  Ein = Loss[25]->getEin(Ein,thickTa);
		  Histo->ECsI_Ta_loss->Fill(Ein-sumEnergy);
		  Histo->ECsI_all_inc->Fill(Ein);
		}
	      if(Z==2)
		{
		  // if(id==0)Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);
		  Ein=sumEnergy;
		  Ein = Loss[ipid+24]->getEin(Ein,thickTa);
		  Histo->ECsI_Ta_loss->Fill(Ein-sumEnergy);
		  //cout << Ein-sumEnergy << endl;
		  Histo->ECsI_all_inc->Fill(Ein);
		}
	    }




	  //if(Z==2 && A==4)
	  //	{
	  //	  cout << Ein << endl;
	  //	}

      if(target==1)
	{
	  thick = 9.472/2./cos(theta); // target thickness !!make target dependent 
	  // deal with loss target dependence too! something like (if  target==n then use Loss[ipid+(n-1)*7)

	      
	  if(Z==1)
	    {
	      if(id==0)
		{
		  Ein = Loss[ipid+28]->getEin(sumEnergy,deadlayer);		  
		  Ein = Loss[ipid+23]->getEin(Ein,thickTa);
		}
	      else
		{		 
		  Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		}

	      Ein = Loss[ipid-1]->getEin(Ein,thick);

	    }
	  if(Z==2)
	    {
	      if(id==0)
		{
		  Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);		  
		  Ein = Loss[ipid+24]->getEin(Ein,thickTa);
		}
	      else
		{		 
		  Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		}

		  
	      Ein = Loss[ipid]->getEin(Ein,thick);
	    }
	  if(Ein==0 && A!=7)cout << "ERROR in Target Energy Addition" << endl;
	  Ekin = Ein;
	}

      else if(target==2)
	{
	  thick = 9.6/2./cos(theta);
	  please  = ipid+(target-1)*8;

	  if(Z==1)
	    {
	      if(id==0)
		{
		  Ein = Loss[ipid+28]->getEin(sumEnergy,deadlayer);		  
		  Ein = Loss[ipid+23]->getEin(Ein,thickTa);
		      
		}
	      else
		{		 
		  Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		      
		}
	      Ein = Loss[please-1]->getEin(Ein,thick);
	    }
	  if(Z==2)
	    {
	      if(id==0)
		{
		  Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);		  
		  Ein = Loss[ipid+24]->getEin(Ein,thickTa);
		}
	      else
		{		 
		  Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		}

	      Ein = Loss[please]->getEin(Ein,thick);
	    }

	  Ekin = Ein;
	}

      else if(target==3)
	{
	  thick = 10.376/2./cos(theta);

	  please  = ipid+(target-1)*8;

	  if(Z==1)
	    {
	      if(id==0)
		{
		  Ein = Loss[ipid+28]->getEin(sumEnergy,deadlayer);		  
		  Ein = Loss[ipid+23]->getEin(Ein,thickTa);
		}
	      else
		{		 
		  Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		}

	      Ein = Loss[please-1]->getEin(Ein,thick);
	    }
	  if(Z==2)
	    {
	      if(id==0)
		{
		  Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);		  
		  Ein = Loss[ipid+24]->getEin(Ein,thickTa);
		}
	      else
		{		 
		  Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		}

	      Ein = Loss[please]->getEin(Ein,thick);
	    }

	  Ekin = Ein;
	}
      else 
	{
	  Ekin = Ein;
	}

      //if(A==7)Histo->Li7_AbsElasXS->Fill(theta);      
      
	  
          Solution[Nsolution].energy = energy;
	  Solution[Nsolution].energyR = Csi.Order[0].energyR;
	  Solution[Nsolution].denergy = Sienergy;
	  Solution[Nsolution].ipie = ipie;
	  Solution[Nsolution].iring = iring;
	  Solution[Nsolution].icsi = Csi.Order[0].strip;
	  Solution[Nsolution].itele = id;
	  Solution[Nsolution].iZ = Z;
	  Solution[Nsolution].iA = A;
	  Solution[Nsolution].mass = A;
	  Solution[Nsolution].Ekin = Ekin;
	  Solution[Nsolution].theta = theta;
	  Solution[Nsolution].phi = phi;
	  Nsolution=1;
	  
	  /*
	  cout << "id=" << id << endl;
	  cout << "ring= " << iring << endl;
	  cout << "csI= " << CsIhit << endl;
	  */
	  



	}
    }
    else if(Csi.Nstore>1 && Pie.Nstore >1 && Ring.Nstore >1)
      {
	//We will need to work on multihit at a later data KB
	//	cout << "HEY!!" << endl;
	multiHitCsi();
      }

  //  Addfake();
  getMomentum(); //Adds energytotal and momentum to Solutions
    

}
//*******************************************************
void telescope::reset()
{

  multPie = 0;
  multRing = 0;


  for(int i =0;i<Nsolution;i++)
    {
      Solution[i].energy = 0;
      Solution[i].energyR = 0;
      Solution[i].denergy = 0;
      Solution[i].ipie =0;
      Solution[i].iring = 0;
      Solution[i].icsi = 0;
      Solution[i].itele = 0;
      Solution[i].iZ = 0;
      Solution[i].iA = 0;
      Solution[i].mass = 0;
      Solution[i].Ekin = 0;
      Solution[i].theta = 0;
      Solution[i].phi = 0;
      Solution[i].energyTot = 0.;
      Solution[i].momentum = 0.;
      Solution[i].velocity = 0.;
      Solution[i].Mvect[0] = 0.;
      Solution[i].Mvect[1] = 0.;
      Solution[i].Mvect[2] = 0.;

    }

  Nsolution =0;

  Pie.reset();
  Ring.reset();
  Csi.reset();



}

//*************************************************************
int telescope::multiHitCsi()
{
  // find number of soultions ,i.e. back and front pairs of strips 
  // with the same energy 
   
  
  int isol = multiHit();
  if (isol <=0) return 0;


  //now assign each of these solutions a Csi detector location 
  int mult[16]={0};  //array for multipility of Si solution for each Csi
  int sil[16][10]; //contains a lits of solutions for each Csi
  for (int i=0;i<Nsolution;i++)
    {
     int ipie = Solution[i].ipie;
     int iring = Solution[i].iring;
     for (int icsi=0;icsi<16;icsi++)
       {
	 if(CsIMap[ipie]==icsi)
	   {
             sil[icsi][mult[icsi]] = i;
             mult[icsi]++;
             Solution[i].icsi = icsi;
             break;
	   }
       }
    }



  //make array of detect csi energies
  float energy[16]={-1.};
  float energyR[16]={-1};
  for (int i=0;i<Csi.mult;i++) 
    {
      energy[Csi.Order[i].strip] = Csi.Order[i].energy; //previously calibrated for proton // 12/9/2015 this is returning values of zero even when there is raw data
      energyR[Csi.Order[i].strip] = Csi.Order[i].energyR; //raw channel number
    }




  //loop over csi location
  for (int icsi = 0;icsi<16;icsi++)
    
   {
      // no solution for this location, ignore
      if (mult[icsi] == 0) continue;
      // more than one si solution for a single Csi location
      else if (mult[icsi] > 1)
	{
          for (int j=0;j<mult[icsi];j++) 
              Solution[sil[icsi][j]].ipid = 0;
          continue;
	}
      // only one si solution for this csi location
      else
	{
	  int ii = sil[icsi][0];
	  //now see if this csi fired 
	  if (energyR[icsi] <= 0.) // 12/9/2015 changed this from energy->energyR should at least count data now, this fixed the initial i ssue
	    {
	      //no csi recorded for this event
	      //stopped in silicon
	      Solution[ii].energy = 0.;
	      Solution[ii].iZ = 0;
	      continue;
	    }
	

	 	  
	  CsIhit = icsi+16*id;


	  Histo->dEE[CsIhit]->Fill(energyR[icsi],Solution[ii].denergy);


	  bool stat = Pid[CsIhit]->getPID(energyR[icsi],Solution[ii].denergy);
	  if(!stat)
	    {
	      Solution[ii].energy =0.;
	      Solution[ii].iZ =0;
	      continue;
	    }
	  //zline check

	  float help;
	  float help2;

	  help = Solution[ii].denergy;
	  help2 = check_intercept[CsIhit]-check_slope[CsIhit]*energyR[icsi];
	  //if(energyR[icsi]>1000&&help<help2){cout << "deuteron contamination in de_E=CsI-" << CsIhit << endl;}
      


	  int Z = Pid[CsIhit]->Z;
	  int A = Pid[CsIhit]->A;




	  int ipie = Solution[ii].ipie;
	  int iring = Solution[ii].iring;

	  block = Block(iring);

	  if(id==1){if(iring>=44)return 0;if(ipie==15||ipie==0)return 0;}	  

	  theta = Pixels.getAngle(id,ipie,iring);
	  //theta = theta+Pixels.deltatheta*(ran->Rndm()-0.5);
	  theta_use = theta*180./3.1415927;
	  phi = Pixels.phi;
	  phi = phi + Pixels.deltaphi*(ran->Rndm()-0.5);



	  xhit = Pixels.r*cos(phi)*sin(theta);
	  yhit = Pixels.r*sin(phi)*sin(theta);

	  double z = Pixels.r*cos(theta);
	  float rperp = sqrt(pow(xhit,2.)+pow(yhit,2.));


	   if(id==0)
	   Histo->HitMap->Fill(-xhit,yhit);
	   else
	   Histo->HitMap1->Fill(-xhit,yhit);

      //I want to look at position dependence, should only really worry about it with single hit data since calibrations are mostly single hit
	   
	   xhit = theta*180/3.1415927*cos(phi);
	   yhit = theta*180/3.1415927*sin(phi);

	   if(A==7)Histo->HitMap_perspective->Fill(-xhit,yhit);

	  if(Z >0 && A>0)
	    {

	      
	      if(Z ==1)
		{
		  if(A ==2)
		    {
		      energy[icsi] = calCsid->getEnergy(0,CsIhit,energyR[icsi]);
		      energy[icsi]=blockCalCsid[CsIhit]->getE(energyR[icsi]+ran->Rndm(),block);
		      if(id==0&&fitCal==1)
		      {
		      energy[icsi]=blockCalCsidRus->getE(energyR[icsi]+ran->Rndm(),CsIhit,theta_use);
		      }
		    }
	    
	    
	    
		  else if(A==3)
		    {
		      energy[icsi] = calCsit->getEnergy(0,CsIhit,energyR[icsi]);
		      //Histo->dEE_ta[CsIhit]->Fill(energyR[icsi],Solution[ii].denergy); // will want to get rid of this eventually (slows stuff down)
		      energy[icsi]=triton_par1*blockCalCsid[CsIhit]->getE(energyR[icsi]+ran->Rndm(),block)-triton_shift+triton_quad*pow(energyR[icsi],2.);
		      if(id==0&&fitCal==1)
		      {
		      energy[icsi]=triton_par1*blockCalCsidRus->getE(energyR[icsi]+ran->Rndm(),CsIhit,theta_use)-triton_shift;
		      }
		    }
		}	      
	      else if(Z == 2)
		{
		  if(A==3)
		    {    energy[icsi] = calCsiHe3->getEnergy(0,CsIhit,energyR[icsi]);

		    }
		  else if(A ==4)
		    {
		      energy[icsi] = calCsiA->getEnergy(0,CsIhit,energyR[icsi]);
		      //Histo->dEE_ta[CsIhit]->Fill(energyR[icsi],Solution[ii].denergy);
		      // energy[icsi]=blockCalCsiA[CsIhit]->getE(energyR[icsi]+ran->Rndm(),block);
		      energy[icsi]=blockCalCsiA[CsIhit]->getE(energyR[icsi]+ran->Rndm(),block)+alpha_quad*pow(energyR[icsi],2.);
		      if(id==0&&fitCal==1)
		      {
		      energy[icsi]=blockCalCsiARus->getE(energyR[icsi]+ran->Rndm(),CsIhit,theta_use);
		      }
		    }
		  else 
		    {
		      cout << "found no calib for Z= " << Z << " and A= " << A << endl; 
		      abort();
		    }
		}
	      else if (Z == 3) 
		{
		  energy[icsi] = calCsiLi6->getEnergy(0,CsIhit,energyR[icsi]);
		}

	      


	      float sumEnergy = energy[icsi] + Solution[ii].denergy;
	      //energy loss in target
	      int ipid =0;
	      ipid=A;

	      /*	      
	      if(Z ==1)
		ipid = A-1;
	      else if(Z==4)
		ipid = 9;
	      else if (Z==6)
		{
		  ipid = A+1;
		}
	      else if(Z==7)
		{
		  ipid =A+3;	      
		}
	      else if(Z==8)
		{
		  ipid = A+6;
		}
	      else 
		ipid = A;
	      */  
	   

	      float thick;
	      float Ein=0;
	      float Ekin=0;
	      int please;


	      float thickTa = 13.8/cos(theta);
	      float deadlayer = 5.8/cos(theta);//11.606/cos(theta); // mg/cm^2 11.606 = 50 um dead layer of Si
	      float thickSi = 121.863/cos(theta);
	      
	      
	      if(id==0)
		{
		  if(Z==1)
		    {
		      sumEnergy = Loss[28+ipid]->getEin(energy[icsi],thickSi);
		      // sumEnergy = energy[icsi]+Solution[ii].denergy;
		    }
		  if(Z==2)
		    {
		      sumEnergy = Loss[29+ipid]->getEin(energy[icsi],thickSi);
		      // sumEnergy = energy[icsi]+Solution[ii].denergy;
		    }
		}
	      


	  if(target==1)
	    {
	      thick = 9.472/2./cos(theta); // target thickness !!make target dependent 
	                                     // deal with loss target dependence too! something like (if                                                 target==n then use Loss[ipid+(n-1)*7)

	      
	      if(Z==1)
		{
		  if(id==0)
		    {
		      //Ein = Loss[ipid+28]->getEin(sumEnergy,deadlayer);		  
		      Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		    }
		  else
		    {		 
		      Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		    }

		  Ein = Loss[ipid-1]->getEin(Ein,thick);

		}
	      if(Z==2)
		{
		  if(id==0)
		    {
		      //Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);		  
		      Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		    }
		  else
		    {		 
		      Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		    }

		  
		  Ein = Loss[ipid]->getEin(Ein,thick);
		}
	      if(Ein==0 && A!=7)cout << "ERROR in Target Energy Addition" << endl;
	      Ekin = Ein;
	    }

	  if(target==2)
	    {
	      thick = 9.6/2./cos(theta);
	      please  = ipid+(target-1)*8;

	      if(Z==1)
		{
		  if(id==0)
		    {
		      Ein = Loss[ipid+28]->getEin(sumEnergy,deadlayer);		  
		      Ein = Loss[ipid+23]->getEin(Ein,thickTa);
		    }
		  else
		    {		 
		      Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		    }
		  Ein = Loss[please-1]->getEin(Ein,thick);
		}
	      if(Z==2)
		{
		  if(id==0)
		    {
		      Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);		  
		      Ein = Loss[ipid+24]->getEin(Ein,thickTa);
		    }
		  else
		    {		 
		      Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		    }

		  Ein = Loss[please]->getEin(Ein,thick);
		}

	      Ekin = Ein;
	    }

	  if(target==3)
	    {
	      thick = 10.376/2./cos(theta);

	      please  = ipid+(target-1)*8;

	      if(Z==1)
		{
		  if(id==0)
		    {
		      Ein = Loss[ipid+28]->getEin(sumEnergy,deadlayer);		  
		      Ein = Loss[ipid+23]->getEin(Ein,thickTa);
		    }
		  else
		    {		 
		      Ein = Loss[ipid+23]->getEin(sumEnergy,thickTa);
		    }

		  Ein = Loss[please-1]->getEin(Ein,thick);
		}
	      if(Z==2)
		{
		  if(id==0)
		    {
		      Ein = Loss[ipid+29]->getEin(sumEnergy,deadlayer);		  
		      Ein = Loss[ipid+24]->getEin(Ein,thickTa);
		    }
		  else
		    {		 
		      Ein = Loss[ipid+24]->getEin(sumEnergy,thickTa);
		    }

		  Ein = Loss[please]->getEin(Ein,thick);
		}

	      Ekin = Ein;
	    }

	   
	      Solution[ii].energy = energy[icsi];
	      Solution[ii].energyR = energyR[icsi];
	      Solution[ii].icsi = Csi.Order[ii].strip;
	      Solution[ii].iZ = Z;
	      Solution[ii].iA = A;
	      Solution[ii].mass = A;
	      Solution[ii].Ekin = Ekin;
	      Solution[ii].theta = theta;
	      Solution[ii].phi = phi;
	      Solution[ii].itele = id;


	    }
	}

   }
  
  return 1;
}
//****************************************************
//recursive subroutine  used for multihit subroutine
void telescope::loop(int depth)
{
  if (depth == NestDim )
    {
      // do stuff here
      int dstrip = 0;
      float de = 0.;
      for (int i=0;i<NestDim;i++)
	{
	  penergy = Pie.Order[NestArray[i]].energy;
	  renergy = Ring.Order[i].energy;
	  
          de += abs(renergy-penergy);
	}


      if (dstrip < dstripMin)
	{
          dstripMin = dstrip;
          for (int i=0;i<NestDim;i++) 
	  arrayD[i] = NestArray[i];
	}


      if (de < deMin)
	{
          deMin = de;
          for (int i=0;i<NestDim;i++) 
	  arrayB[i] = NestArray[i];
	}
      return;

    }

  for (int i=0;i<NestDim;i++)
    {
      NestArray[depth] = i;
      int leave = 0;
      for (int j=0;j<depth;j++)
	{
	  if (NestArray[j] == i)
	    {
	      leave =1;
              break; 
	    } 
	}
      if (leave) continue;
      loop(depth+1);
    }
}

//***************************************************
//extracts multiple particle from strip data 
int telescope::multiHit()
{
  int Ntries = min(Pie.Nstore,Ring.Nstore);
  if (Ntries > 4) Ntries =4;
  Nsolution = 0;
  if (Ntries <= 0) return 0;

  
  for (NestDim = Ntries;NestDim>0;NestDim--)
    {
      dstripMin = 1000;
      deMin = 10000.;
      
      //look for best solution
      loop(0);
      
      //check to see if best possible solution is reasonable
      int leave = 0;
      for (int i=0;i<NestDim;i++)
	{
	  penergy = Pie.Order[arrayB[i]].energy;
	  renergy = Ring.Order[i].energy;

	  float accept = 0.2;
	  if(penergy < 10.) accept =1.5/penergy;
	  
	  if (fabs(penergy-renergy) >penergy*accept)
	    {
	      leave = 1;
	      break;
	    }
	}
      
      if (leave) continue;
      // now load solution
      for (int i=0;i<NestDim;i++)
	{
	  penergy = Pie.Order[i].energy;
	  
	  Solution[i].denergy = penergy;
	  Solution[i].ipie = Pie.Order[i].strip;
	  Solution[i].iring = Ring.Order[arrayB[i]].strip;
          Solution[i].itele = id;
        }

      Nsolution = NestDim;
      
      break;
    }

  return Nsolution;
}
//***********************************************************
void telescope::load(int ipie,int icsi)
{

  CsIMap[ipie] = icsi;

}
//***********************************************************
void telescope::load(int F0low, int F1low,int F2low, int F3low,
                   int F0hi,  int F1hi, int F2hi,  int F3hi,
		   int B0low, int B1low,int B2low, int B3low,
		   int B0hi,  int B1hi, int B2hi,  int B3hi)
{
  FrontLow[0] = F0low;
  FrontLow[1] = F1low;
  FrontLow[2] = F2low;
  FrontLow[3] = F3low;

  FrontHigh[0] = F0hi;
  FrontHigh[1] = F1hi;
  FrontHigh[2] = F2hi;
  FrontHigh[3] = F3hi;

  BackLow[0] = B0low;
  BackLow[1] = B1low;
  BackLow[2] = B2low;
  BackLow[3] = B3low;

  BackHigh[0] = B0hi;
  BackHigh[1] = B1hi;
  BackHigh[2] = B2hi;
  BackHigh[3] = B3hi;

}
//***************************************************************
void telescope::Addfake()
{

  for(int i =0;i<Nsolution;i++)
    {
      if(Solution[i].iZ == 8 && Solution[i].iA ==15)
	{
	  // cout << "solution = " << Nsolution << endl;
	  Solution[Nsolution].ipie = Solution[i].ipie;
	  Solution[Nsolution].iring = Solution[i].iring;
	  Solution[Nsolution].icsi = Solution[i].icsi;
	  Solution[Nsolution].itele = Solution[i].itele;
	  Solution[Nsolution].theta = Solution[i].theta;
	  Solution[Nsolution].phi = Solution[i].phi;
	  Solution[Nsolution].denergy = Solution[i].denergy;

	  float CsI0 = Solution[i].energyR;
	  float fakeE = 0.;
	  int counter = 0;
	  float CsIE = 0.;
	  float dE0 = Solution[i].denergy;
	  float dE = 0.;
	  for(;;)
	    {
	      // CsIE = CsI0 - ran->Rndm()*min((float)500,CsI0);
	      CsIE = CsI0;

	      dE = dE0 - ran->Rndm()*min((float)50,dE0);
	      //dE = dE0;

	      bool stat = Pid[Solution[i].icsi]->getPID(CsIE,dE);
	      if(stat)
		{
		  int Z = Pid[Solution[i].icsi]->Z;
		  int A = Pid[Solution[i].icsi]->A;
		  if(Z == 8 && A ==14)
		    {
		      break;
		    }

		}
	      counter++;
	      if(counter ==100) break;
	    }

	  //if(counter > 30)cout << "counter " <<counter << endl;
	  int CsiHit = Solution[i].icsi + id*4;
	 
	  //fakeE  =calCsi->getEnergy(0,CsiHit,CsIE);
	  // fakeE = calCsiO14->getEnergy(0,CsiHit,fakeE);

	  fakeE = Solution[i].energy;

	  float sumEnergy = fakeE + dE;
	  float thick = 193./2./cos(theta);

	  float Ein = Loss[15]->getEin(sumEnergy,thick);
	  float Ekin = Ein;

	  Solution[Nsolution].Ekin =Ekin;
	  Solution[Nsolution].energy = fakeE;
	  Solution[Nsolution].denergy = dE;
	  Solution[Nsolution].energyR = CsIE;
	  Solution[Nsolution].iZ = 99;
	  Solution[Nsolution].iA = 99;
	  Solution[Nsolution].mass =14;
	  if(counter !=100) Nsolution++;

	}

    }


}

//**********************************************************
void telescope::getMomentum()
{
  for(int i = 0;i<Nsolution;i++)
    {	
      
      float theta = Solution[i].theta;
      float phi = Solution[i].phi;
      float momentum;
      if (relativity)
	{
	  Solution[i].mass *= 931.478;
	  momentum = Solution[i].Kinematics.
	    getMomentum(Solution[i].Ekin,Solution[i].mass);
	  Solution[i].energyTot = Solution[i].Ekin + Solution[i].mass;
	  
	}
      else
	{
	  momentum = sqrt(2.*Solution[i].mass*Solution[i].Ekin);
	  Solution[i].mass = 0.;
	}
      
      Solution[i].Mvect[0] = momentum*sin(theta)*cos(phi);
      Solution[i].Mvect[1] = momentum*sin(theta)*sin(phi);
      Solution[i].Mvect[2] = momentum*cos(theta); 
      Solution[i].momentum = sqrt(pow(Solution[i].Mvect[0],2)
				  +pow(Solution[i].Mvect[1],2)
				  +pow(Solution[i].Mvect[2],2));
      Solution[i].velocity = Solution[i].momentum/Solution[i].energyTot;
      
    }
  
}
