#include "correl.h"
#include <cmath>

using namespace std;
void correl::load(solution* particle)
{

  if (particle->iZ == 1 && particle->iA == 1)
    {
      proton[multProton] = particle;
      if (multProton < 7) multProton++;
      multPro[particle->itele]++;
      //proton[multProton]->A = 1;

    }
  else if (particle->iZ == 1 && particle->iA == 2)
    {
      H2[multH2] = particle;
      if (multH2 < 5) multH2++;
      itele = particle->itele;
      //H2[multH2]->A = 2;
      multh2[itele]++;
    }
  else if (particle->iZ == 1 && particle->iA == 3)
    {

      H3[multH3] = particle;
      if (multH3 < 5) multH3++;
      itele = particle->itele;
      //H3[multH3]->A = 3;
      multh3[itele]++;
    }
 else if (particle->iZ == 2 && particle->iA == 3)
    {
      He3[mult3He] = particle;
      if (mult3He < 5) mult3He++;
      itele = particle->itele;
      //He3[mult3He]->A = 3;
      mult3he[itele]++;
    }
  else if (particle->iZ == 2 && particle->iA ==4)
    {
      alpha[multAlpha] = particle;
      if (multAlpha < 5) multAlpha++;
      itele = particle->itele;
      //alpha[multAlpha]->A = 4;
      multAlp[itele]++;
    }

  else if (particle->iZ == 2 && particle->iA == 6)
    {
      He6[mult6He] = particle;
      if (mult6He < 5) mult6He++;
      itele = particle->itele;
     
      mult6he[itele]++;
    }
  else if (particle->iZ == 2 && particle->iA == 8)
    {
      He8[mult8He] = particle;
      if (mult8He < 5) mult8He++;
      itele = particle->itele;
     
      mult8he[itele]++;
    }
  else if (particle->iZ == 3 && particle->iA == 6)
    {
      Li6[mult6Li] = particle;
      if (mult6Li < 5) mult6Li++;
      itele = particle->itele;
     
      mult6li[itele]++;
    }
  else if (particle->iZ == 3 && particle->iA == 7)
    {
      Li7[mult7Li] = particle;
      if (mult7Li < 5) mult7Li++;
      itele = particle->itele;
     
      mult7li[itele]++;
    }
  else if (particle->iZ == 3 && particle->iA == 8)
    {
      Li8[mult8Li] = particle;
      if (mult8Li < 5) mult8Li++;
      itele = particle->itele;
     
      mult8li[itele]++;
    }
  else if (particle->iZ == 4 && particle->iA == 7)
    {
      Be7[mult7Be] = particle;
      if (mult7Be < 5) mult7Be++;
      itele = particle->itele;
     
      mult7be[itele]++;
    }
  else if (particle->iZ == 4 && particle->iA == 9)
    {
      Be9[mult9Be] = particle;
      if (mult9Be < 5) mult9Be++;
      itele = particle->itele;
     
      mult9be[itele]++;
    }
  else if (particle->iZ == 5 && particle->iA == 10)
    {
      B10[mult10B] = particle;
      if (mult10B < 5) mult10B++;
      itele = particle->itele;
     
      mult10b[itele]++;
    }
  else if (particle->iZ == 5 && particle->iA == 11)
    {
      B11[mult11B] = particle;
      if (mult11B < 5) mult11B++;
      itele = particle->itele;
 
      mult11b[itele]++;
    }
  else if (particle->iZ == 6 && particle->iA == 12)
    {
      C12[mult12C] = particle;
      if (mult12C < 5) mult12C++;
      itele = particle->itele;
     
      mult12c[itele]++;
    }
  else if (particle->iZ == 6 && particle->iA == 11)
    {
      C11[mult11C] = particle;
      if (mult11C < 5) mult11C++;
      itele = particle->itele;
     
      mult11c[itele]++;
    }
  else if (particle->iZ == 6 && particle->iA == 13)
    {
      C13[mult13C] = particle;
      if (mult13C < 5) mult13C++;
      itele = particle->itele;
     
      mult13c[itele]++;
    }
  else if (particle->iZ == 7 && particle->iA == 12)
    {
      N12[mult12N] = particle;
      if (mult12N < 5) mult12N++;
      itele = particle->itele;
     
      mult12n[itele]++;
    }
  else if (particle->iZ == 7 && particle->iA == 13)
    {
      N13[mult13N] = particle;
      if (mult13N < 5) mult13N++;
      itele = particle->itele;
     
      mult13n[itele]++;
    }
  else if (particle->iZ == 7 && particle->iA == 14)
    {
      N14[mult14N] = particle;
      if (mult14N < 5) mult14N++;
      itele = particle->itele;
     
      mult14n[itele]++;
    }
  else if (particle->iZ == 7 && particle->iA == 15)
    {
      N15[mult15N] = particle;
      if (mult15N < 5) mult15N++;
      itele = particle->itele;
     
      mult15n[itele]++;
    }
  else if (particle->iZ == 8 && particle->iA == 13)
    {
      O13[mult13O] = particle;
      if (mult13O < 5) mult13O++;
      itele = particle->itele;
     
      mult13o[itele]++;
    }
  else if (particle->iZ == 8 && particle->iA == 14)
    {
      O14[mult14O] = particle;
      if (mult14O < 5) mult14O++;
      itele = particle->itele;
     
      mult14o[itele]++;
    }
  else if (particle->iZ == 8 && particle->iA == 15)
    {
      O15[mult15O] = particle;
      if (mult15O < 5) mult15O++;
      itele = particle->itele;
     
      mult15o[itele]++;
    }
  else if (particle->iZ == 8 && particle->iA == 16)
    {
      O16[mult16O] = particle;
      if (mult16O < 5) mult16O++;
      itele = particle->itele;
     
      mult15o[itele]++;
    }
  else if (particle->iZ == 98 && particle->iA == 98)
    {
      fakeN13[mult13N] = particle;
      if (multfake13N < 5) multfake13N++;
      itele = particle->itele;
     
      multfake13n[itele]++;
    }
  else if (particle->iZ == 99 && particle->iA == 99)
    {
      fakeO14[mult14O] = particle;
      if (multfake14O < 5) multfake14O++;
      itele = particle->itele;
     
      multfake14o[itele]++;
    }
  /*
  else if (particle->ipid == 50)
    {
      gamma[multGamma] = particle;
      if (multGamma < 5) multGamma++;
      itele = particle->itele;
     
      mult_Gamma[itele]++;
    }
  */

}
//**********************************************************
void correl::setMask()
{
  for (int i=0;i<5;i++)
    {
      maskProton[i] = 1;
      maskAlpha[i] = 1;
      mask2H[i] = 1;
      mask3He[i] = 1;
      mask7Be[i] = 1;
      mask6He[i] = 1;
      mask6Li[i] = 1;
      mask7Li[i] = 1;
      mask7Be[i] = 1;
      mask9Be[i] = 1;
      mask10B[i] = 1;
      mask3H[i] = 1;
    }
}
//************************************************************
void correl::zeroMask()
{
  for (int i=0;i<5;i++)
    {
      maskProton[i] =0;
      maskAlpha[i] = 0;
      mask3He[i] = 0;
      mask2H[i] = 0;
      mask3H[i] = 0;
      mask6He[i] = 0;
      mask6Li[i] = 0;
      mask7Li[i] = 0;
      mask8Li[i] = 0;
      mask7Be[i] = 0;
      mask9Be[i] = 0;
      mask10B[i] = 0;
      mask11B[i] = 0;
      mask13O[i] = 0;
      mask14O[i] = 0;
      maskfake13N[i] = 0;
      maskfake14O[i] = 0;
      mask15O[i] = 0;
      mask16O[i] = 0;
      mask11C[i] = 0;
      mask12C[i] = 0;
      mask13C[i] = 0;
      mask12N[i] = 0;
      mask13N[i] = 0;
      mask14N[i] = 0;
      mask15N[i] = 0;
    }
}
//**********************************************************
void correl::zeroH2()
{
  for (int i=0;i<5;i++) mask2H[i] = 0;
}
//**********************************************************
void correl::zero3He()
{
  for (int i=0;i<5;i++) mask3He[i] = 0;
}
//*********************************************************
void correl::zeroProton()
{
  for (int i=0;i<5;i++) maskProton[i] = 0;
}
//************************************************
void correl::zeroAlpha()
{
  for (int i=0;i<5;i++) maskAlpha[i] = 0;
}
//*******************************************************
  /**
   * make an array of fragments for later analyisis
   */
void correl::makeArray(bool flagMask/*=0*/)
{

 N = 0;

 for (int j=0;j<multProton;j++) 
    {
    if (!flagMask || maskProton[j])
      {
	frag[N] = proton[j];
        N++;
      }
    }
  for (int j=0;j<multAlpha;j++)
    {
    if (!flagMask || maskAlpha[j]) 
      {
	frag[N] = alpha[j];
        N++;
      }
    }
  for (int j=0;j<multH2;j++)
    {
    if (!flagMask || mask2H[j]) 
      {
	frag[N] = H2[j];
	N++;
      }
    }
 for (int j=0;j<multH3;j++)
    {
    if (!flagMask || mask3H[j]) 
      {
	frag[N] = H3[j];
	N++;
      }
    }
 for (int j=0;j<mult3He;j++)
    {
    if (!flagMask || mask3He[j]) 
      {
	frag[N] = He3[j];
	N++;
      }
    } 
 for (int j=0;j<mult6He;j++)
    {
    if (!flagMask || mask6He[j]) 
      {
	frag[N] = He6[j];
	N++;
      }
    } 
 for (int j=0;j<mult6Li;j++)
    {
    if (!flagMask || mask6Li[j]) 
      {
	frag[N] = Li6[j];
	N++;
      }
    } 
 for (int j=0;j<mult7Li;j++)
    {
    if (!flagMask || mask7Li[j]) 
      {
	frag[N] = Li7[j];
	N++;
      }
    } 
 for (int j=0;j<mult8Li;j++)
    {
    if (!flagMask || mask8Li[j]) 
      {
	frag[N] = Li8[j];
	N++;
      }
    } 
 for (int j=0;j<mult7Be;j++)
    {
    if (!flagMask || mask7Be[j]) 
      {
	frag[N] = Be7[j];
	N++;
      }
    } 
 for (int j=0;j<mult11B;j++)
    {
    if (!flagMask || mask11B[j]) 
      {
	frag[N] = B11[j];
	N++;
      }
    } 
 for (int j=0;j<mult10B;j++)
    {
    if (!flagMask || mask10B[j]) 
      {
	frag[N] = B10[j];
	N++;
      }
    } 
 for (int j=0;j<mult11C;j++)
    {
    if (!flagMask || mask11C[j]) 
      {
	frag[N] = C11[j];
	N++;
      }
    } 
 for (int j=0;j<mult12C;j++)
    {
    if (!flagMask || mask12C[j]) 
      {
	frag[N] = C12[j];
	N++;
      }
    } 
 for (int j=0;j<mult13C;j++)
    {
    if (!flagMask || mask13C[j]) 
      {
	frag[N] = C13[j];
	N++;
      }
    } 
 for (int j=0;j<mult12N;j++)
    {
    if (!flagMask || mask12N[j]) 
      {
	frag[N] = N12[j];
	N++;
      }
    } 
 for (int j=0;j<mult13N;j++)
    {
    if (!flagMask || mask13N[j]) 
      {
	frag[N] = N13[j];
	N++;
      }
    } 
 for (int j=0;j<mult14N;j++)
    {
    if (!flagMask || mask14N[j]) 
      {
	frag[N] = N14[j];
	N++;
      }
    } 
 for (int j=0;j<mult15N;j++)
    {
    if (!flagMask || mask15N[j]) 
      {
	frag[N] = N15[j];
	N++;
      }
    } 
 for (int j=0;j<mult13O;j++)
    {
    if (!flagMask || mask13O[j]) 
      {
	frag[N] = O13[j];
	N++;
      }
    } 
 for (int j=0;j<mult14O;j++)
    {
    if (!flagMask || mask14O[j]) 
      {
	frag[N] = O14[j];
	N++;
      }
    } 
 for (int j=0;j<multfake13N;j++)
    {
    if (!flagMask || maskfake13N[j]) 
      {
	frag[N] = fakeN13[j];
	N++;
      }
    } 
 for (int j=0;j<multfake14O;j++)
    {
    if (!flagMask || maskfake14O[j]) 
      {
	frag[N] = fakeO14[j];
	N++;
      }
    } 
 for (int j=0;j<mult15O;j++)
    {
    if (!flagMask || mask15O[j]) 
      {
	frag[N] = O15[j];
	N++;
      }
    } 
 for (int j=0;j<mult16O;j++)
    {
    if (!flagMask || mask16O[j]) 
      {
	frag[N] = O16[j];
	N++;
      }
    } 
}



//***************************************
void correl::Reset()
{
  multProton = 0;
  multAlpha = 0;
  multH2 = 0;
  mult7Be = 0;
  mult6He = 0;
  mult3He = 0;
  mult9Be = 0;
  mult6Li = 0;
  mult7Li = 0;
  mult8Li = 0;
  mult10B = 0;
  mult11B = 0;
  multH3 = 0;
  mult11C =0;
  mult12C =0;
  mult13C =0;
  mult12N =0;
  mult13N =0;
  mult14N =0;
  mult15N =0;
  mult13O =0;
  mult14O =0;
  multfake14O =0;
  multfake13N =0;
  mult15O =0;
  mult16O =0;
  multGamma = 0;
  for (int itele=0;itele<14;itele++) 
    {
      multAlp[itele] = 0;
      multPro[itele] = 0;
      multh2[itele] = 0;
      mult7be[itele] = 0;
      mult6he[itele] = 0;
      mult3he[itele] = 0;
      mult6li[itele] = 0;
      mult7li[itele] = 0;
      mult8li[itele] = 0;
      mult9be[itele] = 0;
      mult10b[itele] = 0;
      mult11b[itele] = 0;
      multh3[itele] = 0;
      mult11c[itele] =0;
      mult12c[itele] =0;
      mult13c[itele] =0;
      mult12n[itele] =0;
      mult13n[itele] =0;
      mult14n[itele] =0;
      mult15n[itele] =0;
      mult13o[itele]=0;
      mult14o[itele]=0;
      multfake13n[itele]=0;
      multfake14o[itele]=0;
      mult15o[itele]=0;
      mult16o[itele]=0;
      mult_Gamma[itele] = 0;
    }
}

//**************************************************************
  /**
   * Finds the total kinetic energy of the fragments
   * in there center-of-mass frame.
   */
float correl::findErel(int target)
{
  
  //first find total momentum
  for (int i=0;i<3;i++) Mtot[i] = 0.;
  float energyTot = 0.;   // total energy for relativity, total mass for newton



  for (int i=0;i<N;i++) 
    {
      energyTot += frag[i]->energyTot;
      for (int j=0;j<3;j++)
        Mtot[j] += frag[i]->Mvect[j];
    }

  momentumCM = 0.;
  for (int j=0;j<3;j++) momentumCM += pow(Mtot[j],2);
  momentumCM = sqrt(momentumCM);

  //transform to average
  velC[0]=0.;//{0.,0.,6.677};
  velC[1]=0.;
  velC[2]=0.;
  //float velC[3]={0.,0.,10.8699}; // 7Be beam
  //float velC[3]={0.,0.,10.7685};// 9C
  //  float velC[3]={0.,0.,10.8381};// 9C //should inquire about this DH

 

  if(target==1)velC[2] = 6.66;//6.62038;//6.66214;//6.63562;//6.66214 //7Li beam after going through half the target...need to make this target dependent!
  if(target==2)velC[2] = 6.66053;
  if(target==3)velC[2] = 6.66187;


  Kinematics.transformMomentum(Mtot,velC,
        energyTot,momC);
  float mmc = 0.;
  for (int i=0;i<3;i++) mmc += pow(momC[i],2);
  mmc = sqrt(mmc);
  cosThetaC = momC[2]/mmc;
  PperpC = sqrt(pow(momC[0],2)+pow(momC[1],2));
  PparaC = momC[2];
  PtotC = sqrt(pow(momC[0],2)+pow(momC[1],2)+pow(momC[2],2));

  //velocity of DECAY center of mass i.e. velocity of particle in LAB frame before it breaks up
  velocityCM = momentumCM*Kinematics.c/energyTot;

  float velCM[3]={0.};
  for (int j=0;j<3;j++) velCM[j] = velocityCM/momentumCM*Mtot[j];
  thetaCM = acos(velCM[2]/velocityCM);
  phiCM = atan2(velCM[1],velCM[0]);

  
  

  float totalKE = 0.;
  for (int i=0;i<N;i++)
    {
      float eKinNew = 
      Kinematics.transformMomentum(frag[i]->Mvect,velCM,
        frag[i]->energyTot,frag[i]->MomCM);
      frag[i]->energyCM = eKinNew - Kinematics.scale*frag[i]->mass;
      totalKE += eKinNew - Kinematics.scale*frag[i]->mass;
    }


  if (N == 3)
    {
      float dot = 0.;
      float mm = 0.;
      for (int j=0;j<3;j++)
	{
         dot += frag[2]->MomCM[j]*momC[j];
         mm += pow(frag[2]->MomCM[j],2);
	}
      mm = sqrt(mm);
      cosAlphaQ = dot/mm/PtotC;
    }

  float mass_target=0.;
  if(target==1)mass_target=9.;
  if(target==2)mass_target=12.;
  if(target==3)mass_target=27.;

  float EPA0 = 24.;  
  float Pbeam2_mean = sqrt(pow((EPA0+931.478)*7.,2)-pow(7.*931.478,2));
  double Emane = (EPA0+931.478)*7 + mass_target*931.478;
  double Ecm = pow( pow(Emane,2.)-pow(Pbeam2_mean,2.) , 0.5); //Ecm = sqrt( (E1+E2)^2-(p1+p2)^2 )
  
  //reaction center of mass velocity
  // float Ecm = (EPA0+931.478)*7. + mass_target*931.478;
  float pcm = Pbeam2_mean;
  float vCM = Pbeam2_mean/Ecm*30. * mass_target/7.; //pcm = plab*m2/Ecm
                                                    //want velocity of 7Li so I divide by 7 (OK in newtonian limit)


  float velReactionCoM[3] ={0,0,velC[2]-vCM}; //velocity I want to transform by (in newtonian limit)
  float  * momCoM;
  momCoM = Kinematics.transformMom(Mtot,velReactionCoM,energyTot,momC);
  float momTot = sqrt(pow(momCoM[0],2.)+pow(momCoM[1],2.)+pow(momCoM[2],2.));
  //float mu = mas_target/(mass_target+7);
  //double thetaReactCoM2 = atan(pcm*sin(thetaCM)/(pcm*cos(thetaCM)-7*931.478*pcm/Emane));//acos(momCoM[2]/momTot);
  thetaReactCoM = acos(momCoM[2]/momTot);
  //cout << thetaReactCoM - thetaReactCoM2 << endl;
 
 // shout out to skisickness.com for a GREAT relativistic derivation of this relationship 
  // Should point out thetaCM is the angle of the reconstructed projectile in the LAB frame (CM refers to the reconstructed projectile) 
  // Confusing right?
  // I would use the other method but I believe my energy calibrations are screwing me again
  // If the total energy is off (which it is I'm sure) then the velocity transform to the CM frame will be wrong
  // Even if total energy is off the RELATIVE energies could be reasonably correct making my data a total wash
  // should try doing what I have simulation at some point 9/21/2016



  return totalKE;
}
//***********************************************************
void correl::getJacobi()
{

  for (int i=0;i<3;i++)
    {
      frag[i]->momentumCM = 0.;
      for (int k=0;k<3;k++) frag[i]->momentumCM +=
	pow(frag[i]->MomCM[k],2);
      frag[i]->momentumCM = sqrt(frag[i]->momentumCM);
    }


  //alpha is the  third fragment
  //first JacobiT
  float dot = 0.;
  float pp[3] = {0.};
  float PP = 0.;
  for (int k=0;k<3;k++)
    {
      pp[k] = frag[0]->MomCM[k] - frag[1]->MomCM[k];
      PP += pow(pp[k],2);
      dot += pp[k]*frag[2]->MomCM[k];
    }
  PP = sqrt(PP);
  cosThetaT = dot/PP/frag[2]->momentumCM;


  dot = 0;
  double PP1 = 0;
  double pp1[3]={0.};
  for (int k=0;k<3;k++)
    {
      pp1[k] = frag[0]->MomCM[k]/frag[0]->mass 
            - frag[2]->MomCM[k]/frag[2]->mass;
      PP1 += pow(pp1[k],2);
      dot += pp1[k]*frag[1]->MomCM[k];
    }  
  PP1 = sqrt(PP1);
  cosThetaY[0] = -dot/PP1/frag[1]->momentumCM;


  dot = 0;
  double PP2 = 0;
  double pp2[3]={0.};
  for (int k=0;k<3;k++)
    {
      pp2[k] = frag[1]->MomCM[k]/frag[1]->mass 
            - frag[2]->MomCM[k]/frag[2]->mass;
      PP2 += pow(pp2[k],2);
      dot += pp2[k]*frag[0]->MomCM[k];
    }  
  PP2 = sqrt(PP2);
  cosThetaY[1] = -dot/PP2/frag[0]->momentumCM;


  cosThetaV = (pp1[0]*pp2[0] + pp1[1]*pp2[1] + pp1[2]*pp2[2])/PP1/PP2;
}
//*****************************************************
void correl::printP()
{

  for (int i=0;i<5;i++)
    {
      frag[i]->momentumCM = 0.;
      for (int k=0;k<3;k++) frag[i]->momentumCM +=
	pow(frag[i]->MomCM[k],2);
      frag[i]->momentumCM = sqrt(frag[i]->momentumCM);
    }

  cout << frag[0]->MomCM[0] << " " << frag[0]->MomCM[1] << " " << 
          frag[0]->MomCM[2] <<endl;

  cout << frag[1]->MomCM[0] << " " << frag[1]->MomCM[1] << " " << 
          frag[1]->MomCM[2] <<endl;


  cout << frag[2]->MomCM[0] << " " << frag[2]->MomCM[1] << " " << 
          frag[2]->MomCM[2] <<endl;

  cout << frag[3]->MomCM[0] << " " << frag[3]->MomCM[1] << " " << 
          frag[3]->MomCM[2] <<endl;


  cout << frag[4]->MomCM[0] << " " << frag[4]->MomCM[1] << " " << 
    frag[4]->MomCM[2] << " " << frag[4]->momentumCM << " " <<
    Kinematics.getKE(frag[4]->momentumCM,frag[4]->mass) <<  endl;




  cout << endl;

    
}
//****************************************************************
float correl::getAlphaMom()
{
     frag[4]->momentumCM = 0.;
      for (int k=0;k<3;k++) frag[4]->momentumCM +=
	pow(frag[4]->MomCM[k],2);
      frag[4]->momentumCM = sqrt(frag[4]->momentumCM);
      return frag[4]->momentumCM;
}
//********************************************************************
void correl::rotate()
{
  float phi = atan2(frag[4]->MomCM[1],frag[4]->MomCM[0]);

  //rotate in x-y plane
  for (int j=0;j<5;j++)
    {
      frag[j]->MomRot[0] = frag[j]->MomCM[0]*cos(phi) + frag[j]->MomCM[1]*sin(phi);
      frag[j]->MomRot[1] = frag[j]->MomCM[1]*cos(phi) - frag[j]->MomCM[0]*sin(phi);
      frag[j]->MomRot[2] = frag[j]->MomCM[2];
    }

  float theta = acos(frag[4]->MomCM[2]/frag[4]->momentumCM);


  //rotate in x-y plane
  for (int j=0;j<5;j++)
    {
      frag[j]->MomRot2[2] = frag[j]->MomRot[2]*cos(theta) + frag[j]->MomRot[0]*sin(theta);
      frag[j]->MomRot2[0] = frag[j]->MomRot[0]*cos(theta) - frag[j]->MomRot[2]*sin(theta);
      frag[j]->MomRot2[1] = frag[j]->MomRot[1];
    }



}
void correl::find6Be(int i1,int i2,int j1, int j2)
{
  //first find total momentum
  float Mtot[3]={0.};
  float energyTot = 0.;   // total energy for relativity, total mass for newton

  energyTot += alpha[0]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += alpha[0]->Mvect[j];


  energyTot += proton[i1]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[i1]->Mvect[j];

  energyTot += proton[i2]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[i2]->Mvect[j];

  float momentum = 0.;
  for (int j=0;j<3;j++) momentum += pow(Mtot[j],2);
  momentum = sqrt(momentum);



  


  //velocity of 6Be
  float velocity6Be = momentum*Kinematics.c/energyTot;
  float vel6Be[3]={0.};
  for (int j=0;j<3;j++) 
    {
     vel6Be[j] = velocity6Be/momentum*Mtot[j];
    }


  float totalKE_6Be = 0.;
  float eKinNew =  Kinematics.transformMomentum(alpha[0]->Mvect,vel6Be,
        alpha[0]->energyTot,alpha[0]->MomCM);
  totalKE_6Be += eKinNew - Kinematics.scale*alpha[0]->mass;

  eKinNew =  Kinematics.transformMomentum(proton[i1]->Mvect,vel6Be,
        proton[i1]->energyTot,proton[i1]->MomCM);
  totalKE_6Be += eKinNew - Kinematics.scale*proton[i1]->mass;

  eKinNew =  Kinematics.transformMomentum(proton[i2]->Mvect,vel6Be,
        proton[i2]->energyTot,proton[i2]->MomCM);
  totalKE_6Be += eKinNew - Kinematics.scale*proton[i2]->mass;



  solution be6;
  be6.energyTot = 0.;
  be6.mass = 6.*Kinematics.nMass + totalKE_6Be;
  for (int j=0;j<3;j++) 
    {
     be6.Mvect[j] = Mtot[j];
     be6.energyTot += pow(be6.Mvect[j],2);
    }
  be6.energyTot = sqrt(be6.energyTot + pow(be6.mass,2));


  //now worry aboyt 8C->2p decat

  //first find total momentum
  for (int jj=0;jj<3;jj++)Mtot[jj]=0.;
  energyTot = 0.;   // total energy for relativity, total mass for newton

  energyTot += be6.energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += be6.Mvect[j];


  energyTot += proton[j1]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[j1]->Mvect[j];

  energyTot += proton[j2]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[j2]->Mvect[j];

  momentum = 0.;
  for (int j=0;j<3;j++) momentum += pow(Mtot[j],2);
  momentum = sqrt(momentum);


 //velocity of 8C
  float velocity8C = momentum*Kinematics.c/energyTot;
  float vel8C[3]={0.};
  for (int j=0;j<3;j++) 
    {
     vel8C[j] = velocity8C/momentum*Mtot[j];
    }


  float totalKE_8C = 0.;
  eKinNew =  Kinematics.transformMomentum(be6.Mvect,vel8C,
        be6.energyTot,be6.MomCM);
  totalKE_8C += eKinNew - Kinematics.scale*be6.mass;

  eKinNew =  Kinematics.transformMomentum(proton[j1]->Mvect,vel8C,
        proton[j1]->energyTot,proton[j1]->MomCM);
  totalKE_8C += eKinNew - Kinematics.scale*proton[j1]->mass;



  eKinNew =  Kinematics.transformMomentum(proton[j2]->Mvect,vel8C,
        proton[j2]->energyTot,proton[j2]->MomCM);
  totalKE_8C += eKinNew - Kinematics.scale*proton[j2]->mass;



  be6.momentumCM = 0.;
  proton[j1]->momentumCM = 0.;
  proton[j2]->momentumCM = 0.;
  for (int k=0;k<3;k++) 
    {
      be6.momentumCM += pow(be6.MomCM[k],2);
      proton[j1]->momentumCM += pow(proton[j1]->MomCM[k],2);
      proton[j2]->momentumCM += pow(proton[j2]->MomCM[k],2);
    }  
  be6.momentumCM  = sqrt(be6.momentumCM);
  proton[j1]->momentumCM = sqrt(proton[j1]->momentumCM);
  proton[j2]->momentumCM = sqrt(proton[j2]->momentumCM);

  float dot = 0.;
  float pp[3] = {0.};
  float PP = 0.;
  for (int k=0;k<3;k++)
    {
      pp[k] = proton[j1]->MomCM[k] - proton[j2]->MomCM[k];
      PP += pow(pp[k],2);
      dot += pp[k]*be6.MomCM[k];
    }
  PP = sqrt(PP);
  cosThetaT = dot/PP/be6.momentumCM;

   dot = 0;
  PP = 0;
  for (int k=0;k<3;k++)
    {
      pp[k] = proton[j1]->MomCM[k]*6. - be6.MomCM[k];
      PP += pow(pp[k],2);
      dot += pp[k]*proton[j2]->MomCM[k];
    }  
  PP = sqrt(PP);
  cosThetaY[0] = dot/PP/proton[j2]->momentumCM;


  dot = 0;
  PP = 0;
  for (int k=0;k<3;k++)
    {
      pp[k] = proton[j2]->MomCM[k]*6. - be6.MomCM[k];
      PP += pow(pp[k],2);
      dot += pp[k]*proton[j1]->MomCM[k];
    }  
  PP = sqrt(PP);
  cosThetaY[1] = dot/PP/proton[j1]->momentumCM;


  //first find total momentum
  for (int jj=0;jj<3;jj++)Mtot[jj]=0.;
  energyTot = 0.;   // total energy for relativity, total mass for newton

  energyTot += proton[j1]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[j1]->Mvect[j];

  energyTot += proton[j2]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[j2]->Mvect[j];

  momentum = 0.;
  for (int j=0;j<3;j++) momentum += pow(Mtot[j],2);
  momentum = sqrt(momentum);


 //velocity of pp
  float velocitypp = momentum*Kinematics.c/energyTot;
  float velpp[3]={0.};
  for (int j=0;j<3;j++) 
    {
     velpp[j] = velocitypp/momentum*Mtot[j];
    }


  float totalKE_pp = 0.;
  eKinNew =  Kinematics.transformMomentum(proton[j1]->Mvect,velpp,
        proton[j1]->energyTot,proton[j1]->MomCM);
  totalKE_pp += eKinNew - Kinematics.scale*proton[j1]->mass;

  eKinNew =  Kinematics.transformMomentum(proton[j2]->Mvect,velpp,
        proton[j2]->energyTot,proton[j2]->MomCM);
  totalKE_pp += eKinNew - Kinematics.scale*proton[j2]->mass;



  x_T = totalKE_pp/totalKE_8C;



  //first find total momentum
  for (int jj=0;jj<3;jj++)Mtot[jj]=0.;
  energyTot = 0.;   // total energy for relativity, total mass for newton

  energyTot += be6.energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += be6.Mvect[j];

  energyTot += proton[j2]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[j2]->Mvect[j];

  momentum = 0.;
  for (int j=0;j<3;j++) momentum += pow(Mtot[j],2);
  momentum = sqrt(momentum);


 //velocity of pbe
  float velocitypbe = momentum*Kinematics.c/energyTot;
  float velpbe[3]={0.};
  for (int j=0;j<3;j++) 
    {
     velpbe[j] = velocitypbe/momentum*Mtot[j];
    }


  float totalKE_pbe = 0.;
  eKinNew =  Kinematics.transformMomentum(be6.Mvect,velpbe,
        be6.energyTot,be6.MomCM);
  totalKE_pbe += eKinNew - Kinematics.scale*be6.mass;

  eKinNew =  Kinematics.transformMomentum(proton[j2]->Mvect,velpbe,
        proton[j2]->energyTot,proton[j2]->MomCM);
  totalKE_pbe += eKinNew - Kinematics.scale*proton[j2]->mass;



  x_Y[1] = totalKE_pbe/totalKE_8C;

  //first find total momentum
  for (int jj=0;jj<3;jj++)Mtot[jj]=0.;
  energyTot = 0.;   // total energy for relativity, total mass for newton

  energyTot += proton[j1]->energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += proton[j1]->Mvect[j];

  energyTot += be6.energyTot;
  for (int j=0;j<3;j++)
  Mtot[j] += be6.Mvect[j];

  momentum = 0.;
  for (int j=0;j<3;j++) momentum += pow(Mtot[j],2);
  momentum = sqrt(momentum);


 //velocity of pbe
  velocitypbe = momentum*Kinematics.c/energyTot;
  for (int jj=0;jj<3;jj++)velpbe[jj]=0.;
  for (int j=0;j<3;j++) 
    {
     velpbe[j] = velocitypbe/momentum*Mtot[j];
    }


  totalKE_pbe = 0.;
  eKinNew =  Kinematics.transformMomentum(proton[j1]->Mvect,velpbe,
        proton[j1]->energyTot,proton[j1]->MomCM);
  totalKE_pbe += eKinNew - Kinematics.scale*proton[j1]->mass;

  eKinNew =  Kinematics.transformMomentum(be6.Mvect,velpbe,
        be6.energyTot,be6.MomCM);
  totalKE_pbe += eKinNew - Kinematics.scale*be6.mass;


  x_Y[0] = totalKE_pbe/totalKE_8C;


}
