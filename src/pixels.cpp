#include "pixels.h"


pixels::pixels()
{

  targetdist[0] =  149;//150.; //in mm for Russian
  targetdist[1] =  354.;//358.; //in mm for S2


  //Zero the arrays
  for(int i=0;i<64;i++)
    {
      for(int j =0;j<48;j++)
	{
	  for(int k=0;k<2;k++)
	    {
	      TeleP[k].Location[i][j].theta=0;
	      TeleP[k].Location[i][j].phi=0;
	      TeleP[k].Location[i][j].x=0;
	      TeleP[k].Location[i][j].y=0;
	      TeleP[k].Location[i][j].z=0;
	      TeleP[k].Location[i][j].r=0;
	    }
	}
    }
    

  double theta,phi;
  //Get Theta and Phi for Russian Si
  for(int ipie =0;ipie<64;ipie++)
    {
      for(int iring=0;iring<32;iring++)
	{

	  double rrus = 15.4296875+0.859375*iring; //radius of each ring
	  TeleP[0].Location[ipie][iring].z = targetdist[0];
	  TeleP[0].Location[ipie][iring].theta = atan2(rrus,targetdist[0]);
	  TeleP[0].Location[ipie][iring].deltatheta = atan2(15+(iring+1.)*0.859375,targetdist[0])-atan2(15+iring*0.859375,targetdist[0]);
	  TeleP[0].Location[ipie][iring].steradian = 2*3.14159*(cos( atan2(rrus-0.859375,targetdist[0]) ) - cos( atan2(rrus+0.859375,targetdist[0]) ));

	  if(ipie <16)
	    {
	      TeleP[0].Location[ipie][iring].phi = (2.8125+5.625*(15-(ipie)))*acos(-1)/180.;
	      TeleP[0].Location[ipie][iring].deltaphi = 5.625*acos(-1)/180.;
	    }
	  else
	    {
	      TeleP[0].Location[ipie][iring].phi = (2.8125+5.625*(63-(ipie)+16))*acos(-1)/180.;
	      TeleP[0].Location[ipie][iring].deltaphi = 5.625*acos(-1)/180.;
	    }
	  
	  double z =TeleP[0].Location[ipie][iring].z;
	  TeleP[0].Location[ipie][iring].r = sqrt(z*z+rrus*rrus);
	  
	}
    }
  //Get Theta and Phi for S2 Si
  for(int ipie =0;ipie<16;ipie++)
    {
      for(int iring=0;iring<48;iring++)
	{

	  double rs2 = 11.25+0.5*iring; //radius of each ring
	  TeleP[1].Location[ipie][iring].z = targetdist[1];
	  TeleP[1].Location[ipie][iring].theta = atan2(rs2,targetdist[1]);
	  TeleP[1].Location[ipie][iring].deltatheta = atan2(11+(iring+1.)*0.859375,targetdist[1])-atan2(11+iring*0.859375,targetdist[1]);
	  TeleP[0].Location[ipie][iring].steradian = 2*3.14159*(cos( atan2(rs2-0.5,targetdist[1]) ) - cos( atan2(rs2+0.5,targetdist[1]) ));
	 
	  if(ipie >=12)
	    {
	      TeleP[1].Location[ipie][iring].phi = (11.25+22.5*(ipie-12))*acos(-1)/180.;
	      TeleP[1].Location[ipie][iring].deltaphi = 22.5*acos(-1)/180.;
	    }
	  else
	    {
	      TeleP[1].Location[ipie][iring].phi = (101.25+22.5*(ipie))*acos(-1)/180.;
	      TeleP[1].Location[ipie][iring].deltaphi = 22.5*acos(-1)/180.;
	    }
	  
	  double z =TeleP[1].Location[ipie][iring].z;
	  TeleP[1].Location[ipie][iring].r = sqrt(z*z+rs2*rs2);
	  
	}
    }


  //these are for position dependent calibrations used since they depend on the angle used
  double rr;
  for(int iring =0;iring<32;iring++)
    {
      rr = 15.4296875+0.859375*iring;
      thetaRus[iring] = atan2(rr,149);
    }
  for(int iring =0;iring<48;iring++)
    {
      rr = 11.25+0.5*iring;
      thetaS2[iring] = atan2(rr,358);
    }


}
//*********************************
location pixels::getCenter(int itele)
{
  float x = (TeleP[itele].Location[15][15].x 
          + TeleP[itele].Location[15][16].x
          + TeleP[itele].Location[16][15].x
          + TeleP[itele].Location[16][16].x)/4.;

  float y = (TeleP[itele].Location[15][15].y
          + TeleP[itele].Location[15][16].y
          + TeleP[itele].Location[16][15].y
          + TeleP[itele].Location[16][16].y)/4.;

  float z = (TeleP[itele].Location[15][15].z 
          + TeleP[itele].Location[15][16].z
          + TeleP[itele].Location[16][15].z
          + TeleP[itele].Location[16][16].z)/4.;
 
  location out;
  out.x = x;
  out.y = y;
  out.z = z;
  
  return out;
}
//***************************************
float pixels::getAngle(int itele, int ipie, int iring)
{
  r = TeleP[itele].Location[ipie][iring].r;
  phi = TeleP[itele].Location[ipie][iring].phi;
  deltaphi = TeleP[itele].Location[ipie][iring].deltaphi;
  deltatheta = TeleP[itele].Location[ipie][iring].deltatheta;
  return TeleP[itele].Location[ipie][iring].theta;
  
}
//*********************************
float pixels::getCsiCenter(int itele, int iCsi)
{
  float x,y,z;
  if (iCsi == 0)
    {
     x = (TeleP[itele].Location[7][7].x 
          + TeleP[itele].Location[7][8].x
          + TeleP[itele].Location[8][7].x
          + TeleP[itele].Location[8][8].x)/4.;

     y = (TeleP[itele].Location[7][7].y
          + TeleP[itele].Location[7][8].y
          + TeleP[itele].Location[8][7].y
          + TeleP[itele].Location[8][8].y)/4.;

     z = (TeleP[itele].Location[7][7].z 
          + TeleP[itele].Location[7][8].z
          + TeleP[itele].Location[8][7].z
          + TeleP[itele].Location[8][8].z)/4.;
    }
  if (iCsi == 1)
    {
     x = (TeleP[itele].Location[7][23].x 
          + TeleP[itele].Location[7][24].x
          + TeleP[itele].Location[8][23].x
          + TeleP[itele].Location[8][24].x)/4.;

     y = (TeleP[itele].Location[7][23].y
          + TeleP[itele].Location[7][24].y
          + TeleP[itele].Location[8][23].y
          + TeleP[itele].Location[8][24].y)/4.;

     z = (TeleP[itele].Location[7][23].z 
          + TeleP[itele].Location[7][24].z
          + TeleP[itele].Location[8][23].z
          + TeleP[itele].Location[8][24].z)/4.;
    }
  if (iCsi == 2)
    {
     x = (TeleP[itele].Location[23][23].x 
          + TeleP[itele].Location[23][24].x
          + TeleP[itele].Location[24][23].x
          + TeleP[itele].Location[24][24].x)/4.;

     y = (TeleP[itele].Location[23][23].y
          + TeleP[itele].Location[23][24].y
          + TeleP[itele].Location[24][23].y
          + TeleP[itele].Location[24][24].y)/4.;

     z = (TeleP[itele].Location[23][23].z 
          + TeleP[itele].Location[23][24].z
          + TeleP[itele].Location[24][23].z
          + TeleP[itele].Location[24][24].z)/4.;
    }
  if (iCsi == 3)
    {
     x = (TeleP[itele].Location[23][7].x 
          + TeleP[itele].Location[23][8].x
          + TeleP[itele].Location[24][7].x
          + TeleP[itele].Location[24][8].x)/4.;

     y = (TeleP[itele].Location[23][7].y
          + TeleP[itele].Location[23][8].y
          + TeleP[itele].Location[24][7].y
          + TeleP[itele].Location[24][8].y)/4.;

     z = (TeleP[itele].Location[23][7].z 
          + TeleP[itele].Location[23][8].z
          + TeleP[itele].Location[24][7].z
          + TeleP[itele].Location[24][8].z)/4.;
    }

  float r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  return acos(z/r)*180./3.14159;
}
