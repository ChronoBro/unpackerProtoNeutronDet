#include <string>

using namespace std;

struct order
{
  float energy;  //high gain MeV
  float energyR;  // high gain channels
  float energyRlow;  //low gain channels
  float energylow; //low gain MeV
  float energyMax;
  int strip;
  int overflow;
  int neighbours;
  float time;
};

int const nnn=32;

/**
 * !\brief Energy ordered list
 *
 * This class creates an energy ordered list of the strips
 * read out from a strip detector, keeping track of the strip 
 * numbers that fired.
 */

class elist
{
 public:
  int Nstore; //number stored in list
  order Order[nnn];
  int mult;
  void Add(int,int,float,float);
  void Add(int,int,float,float,float);
  int  Reduce(char*);
  void reset();
  void Neighbours(string,float,float);
};
