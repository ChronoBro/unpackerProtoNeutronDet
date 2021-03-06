//creates energy order (descending) lists of strips

#include "elist.h"
#include <algorithm>
#include <iostream>
using namespace std;

//places a new strip energy in an order list from max to min energy
void elist::Add(int StripNum, int underOver, float energy,  float time)
{

  //first find place in list
  int i = 0;
  for (;;)
    {
      if (i == Nstore) break;
      if (energy > Order[i].energy)break;
      i++;
    }
  if (i == nnn) return; // not enougth room in list 

  //new list length
  int N = min(nnn,Nstore+1);

  // move those in energy below the new value down the list 
  for (int j=N-1;j>i;j--) Order[j] = Order[j-1];

  //add present energy to list

  Order[i].energy = energy;
  Order[i].overflow = underOver;
  Order[i].strip = StripNum;
  Order[i].time = time;
  Order[i].energyRlow = 0.;  


  // increase list length
  Nstore = N;
  mult = N;
}
//***********************************************************************
//places a new strip energy in an order list from max to min energy
void elist::Add(int StripNum, int underOver, float energy,float energyR,float time)
{

  //first find place in list
  int i = 0;
  for (;;)
    {
      if (i == Nstore) break;
      if (energy > Order[i].energy)break;
      i++;
    }
  if (i == nnn) return; // not enougth room in list 

  //new list length
  int N = min(nnn,Nstore+1);

  // move those in energy below the new value down the list 
  for (int j=N-1;j>i;j--) Order[j] = Order[j-1];

  //add present energy to list

  Order[i].energy = energy;
  Order[i].energyR = energyR;
  Order[i].overflow = underOver;
  Order[i].strip = StripNum;
  Order[i].time = time;
  Order[i].energyRlow = 0.;
  Order[i].energylow = 0.;


  // increase list length
  Nstore = N;
  mult = N;
}

//**********************************************************************
//remove silicon strip cross talk from list 
int elist::Reduce(char*face)
{

  if (Nstore <= 1) return Nstore;
  for (int ii = Nstore-1;ii > 0;ii--) // ii location in list of interest
    {
     int xtalk = 0;
     if (Nstore <= 1) return Nstore;
     if (Nstore > nnn) 
       {
	 cout << "problem in Reduce" << endl;
         return 0;
       }
     if (Nstore < 0)
       {
	 cout << "problem in Reduce()" << endl;
	 return 0;
       }
     for (int i=ii-1;i>=0;i--)
       {
         //look for neighboring strips
         if (abs(Order[ii].strip-Order[i].strip) == 1)
	   {
	     if (*face == 'F')
	       {
                 if (Order[ii].energy < Order[i].energy*0.044
                      && Order[i].energy > 10.)
	           {
	             xtalk = 1; // signal xtalk
                     break;
	           }
	       }
	     else if (*face == 'B')
	       {
                 if (Order[ii].energy < Order[i].energy*0.30
                      && Order[i].energy > 10.)
	           {
	             xtalk = 1; // signal xtalk
                     break;
	           }
	       }
  	   }
        }
     if (xtalk == 0) continue;
     //remove from list
     if (ii+1 < Nstore) 
       {
       for (int j=ii;j<Nstore;j++)Order[j] = Order[j+1];
       }
     Nstore--;
    }   

  return Nstore;
}
//*********************************************************************
void elist::reset()
{
  Nstore = 0;
}
//*********************************************************************
  //looks for cross talks events
void elist::Neighbours(string side,float factor1,float factor2)
{

  if (Nstore <1) return; // nothing to look at 
  if (Nstore ==1) // no neighbors to look for
    {
      Order[0].energyMax = Order[0].energy;
      Order[0].neighbours = 0;
      return ;
    }
  int i=-1;
  for (;;)
    {
      i++;
      if (i >= Nstore)break;
      Order[i].energyMax = Order[i].energy;
      Order[i].neighbours = 0;

      int j = i;
      for (;;)
	{
          j++;
          if (j >= Nstore) break;
	  if (abs(Order[i].strip - Order[j].strip) == 1) //neighboring strips
	    {
              Order[i].energy += Order[j].energy; // add energy from
	                                             // adjacent strip
              Order[i].neighbours++;

              //remove this strip from the list
              if (j != Nstore-1)
		{
	         for (int k=j+1;k<Nstore;k++)Order[k-1] = Order[k];
		}
	      Nstore--;
              j--;
	    }
	}
    }
  for (int i=0;i<Nstore;i++) 
    {
     if (Order[i].overflow) 
         Order[i].energy = Order[i].energyMax;
     //to make on average the Zlines lie at the same location 
     //of events with no neighboring strips, make a small shift 
     //in the total energy
     else if (Order[i].energy != Order[i].energyMax) 
       {
         if (side=="Front")Order[i].energy -= .5;
           else if (side == "Back")
	     { 
     //  if (i!=2) Order[i].energy*=(int)pow(factor1,Order[i].neighbours);
      //else Order[i].energy*=factor2;
	     }
       }
    }
}
