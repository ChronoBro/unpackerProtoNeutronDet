// C++ file to read MSU event files 
//file numbers.beam contains runs to sort
//uses class hira to unpack hira data
//write out spectra in file sort.root 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "det.h"
#include <ctime>
#include "histo.h"
using namespace std;



int main(int argc,char *argv[])
{

  histo * Histo = new histo();
  unsigned short *point,*epoint;
  int unsigned words;
  int unsigned type;

  int physicsEvent = 0;
  int physicsEventGood = 0;
  int physicsEventCounter = 0;
  int scalerBuffer = 0;
  int Npauses = 0;
  int Nresumes = 0;
  int runno = 0;

 

  det Det(Histo);
  Det.Nalphat = 0;
  Det.N7Li_ta =0;

  Det.Silicon->Ncross=0;


  ostringstream outstring;
  int number;
  ifstream evtfile;
  FILE *pFile;
  bool fileProblem = false;
  bool endOfRun = false;
  bool first = true;  

  int raw_trig =0;
  int live_trig = 0;
  int Si_OR = 0;
  int CsI_OR = 0;
  int XLM_Complete =0;
  int Rus_Pie_OR = 0;
  int S2_Pie_OR = 0;
  int Rus_CsI_OR = 0;
  int S2_CsI_OR = 0;
  int FC = 0;

  //open file with run numbers to sort
  ifstream runFile;
  runFile.open("numbers.beam");
  //check if this file exists
  if (runFile.is_open() == 0 )
    {
      cout << " could open runfile " << endl;
      return 1;
    }
  for (;;)  // loop over run numbers
    {
      
      if (evtfile.is_open()) 
        cout << "problem previous file not closed" << endl;
      

      if(argc==1)
	{
	  runFile >> number;
      
	  //check to see if at end of file
	  if (runFile.eof())break;
	  if (runFile.bad())break; 
	}
      else if(argc==2)
	{
	  number = atoi(argv[1]);

	}

      for (int iExtra=0;iExtra<3;iExtra++) //loop over split evtfiles
	{
	  
	  //the following loop accounts for files that were split
	  endOfRun=false;
	  fileProblem = 0;
	  outstring.str("");
	  

	  if(number < 10)
	      outstring << "/data1/TAMU_aug2015/data/run-000" << number;
	  else if(number >=10 && number < 100)
	    outstring << "/data1/TAMU_aug2015/data/run-00" << number;
	  else 
	    outstring << "/data1/TAMU_aug2015/data/run-0"<<number;

	  	    
	  if (iExtra == 0)
	    outstring<<"-00.evt";
	  else 
	    outstring<<"_"<<iExtra<<"-13328.evt";
	   
         string name = outstring.str();

         //open evt file


         evtfile.clear();
         evtfile.open(name.c_str(),ios::binary);


	 //check to see if there are extra files     
         if (iExtra>0 && !fileProblem && !evtfile) 
	   {
             break;
	   }

         cout << '\n'<<name << endl;

         if (evtfile.bad()) cout << "bad" << endl;
         if (evtfile.fail()) cout << "fail" << endl;



         if (!evtfile)
           {
	     cout << "could not open event file" << endl;
             return 1;
           }
    /* 0=blank  1=9Be    2=natuarl C     3=27Al  4=Au  5=1/16" Al 6=1/8" Al*/


	 //	 Det.Silion.Telescope->variable = runno;

	 for(;;)  // loop over items in a evtfile
	   {
            int const hBufferWords = 4;
            int const hBufferBytes = hBufferWords*2;
            unsigned short hBuffer[hBufferWords];
   	    evtfile.read((char*)hBuffer,hBufferBytes);

	    if(evtfile.eof())
	      {
                cout << "eof found" << endl;
		fileProblem = true;
                break;
	      }
	    if(evtfile.bad())
	      {
                cout << " bad found" << endl;
                fileProblem = true;
                break;
	      }

	    point = hBuffer;
            int nbytes = *point++;
            int nbytes2 = *point++;
            int type = *point++;
            int type2 = *point;

            //cout << nbytes << " " << nbytes2 << " " << type << " " <<type2 << endl; 

            int dBufferBytes = nbytes - 8;
            int dBufferWords = dBufferBytes/2;

  	    //cout << "bytes= " << dBufferBytes << ", words " << dBufferWords << endl;
            unsigned short dBuffer[dBufferWords];
            evtfile.read((char*)dBuffer,dBufferBytes);
            point = dBuffer;

	    if (type == 1)
	      {

               runno = *point;
               cout << "run number = " << runno << endl; 
	      }
            else if (type == 30)
	      {
		physicsEvent++;
                if (physicsEvent%30000 == 0) 
                     cout << '\xd'<< physicsEvent << flush;

		//		cout << "events= " << physicsEvent << endl;
		
			
		/*	
		if(physicsEvent == 1)
		  {

		    cout << "here" << endl;
		    
		      unsigned short * pig = point;
		      for (int i=0;i<200;i++)
		      {
		      cout << dec<< *pig << " " << hex<< *pig << endl;
		      pig++;
		      }
		      abort();
		  }
		*/
		    
		/*	
		if(physicsEvent ==4535)
		  {
		    cout << "Bad Event found" << endl;
		    
		    unsigned short * pig = point;
		    unsigned short * pig2 = point;
		    for(int i=1;i<81;i++)
		      {
			cout << std::setfill('0') << setw(4)  <<hex << *pig2 << " " ;
			if(i%8 ==0)
			  {
			    cout << endl;
			  }
			pig2++;
		      }
		    cout << endl;
		    abort();
		    for (int i=0;i<80;i++)
		      {
		    	cout << dec<< *pig << " " << hex<< *pig << endl;
		     	pig++;
		      }
		    abort();
		    
		  }
		*/
		
		//	cout << "event = " << physicsEvent << endl;
		
		unsigned short * badpig = point;
                //cout << "in unpack " << endl;
		bool stat  =false;
		//if(physicsEvent == 625149)

		Det.Silicon->eventNum = physicsEvent;

		stat = Det.unpack(point,runno); //start the data unpack

		//cout << "out unpack " << endl;
                if (stat) 
		  {
		


		    
		    physicsEventGood++;
		    Det.analyze(physicsEvent);
		    
		  }
		/*		
		  else 
		  {
		  cout <<"Event " << physicsEvent << endl;
		  for (int i=0;i<200;i++)
		  {
		  cout << dec<< *badpig << " " << hex<< *badpig << endl;
		  badpig++;
		  }
		  abort();
		  
		  }
		*/
	      }
	    else if (type == 31)
	      {
                physicsEventCounter++;
	      }
            else if (type == 2)
	      {
                endOfRun = true;
                break;
	      }
	    else if (type == 20)
	      {
                scalerBuffer++;

		unsigned short * scaler = point;
		unsigned short size1;
		unsigned short size2;
		unsigned long size;
		scaler +=8; //skipping to the data

		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		live_trig += size;

		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		raw_trig += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;

		Si_OR += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		CsI_OR += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		XLM_Complete += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		Rus_Pie_OR += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		S2_Pie_OR += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		Rus_CsI_OR += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		S2_CsI_OR += size;
		size1 = *scaler++;
		size2 = *scaler++;
		size = (size2 << 16) | size1;
		FC += size;



		/*
		unsigned short * pig = point;
		for (int i=0;i<80;i++)
		  {
		    cout << dec<< *pig << " " << hex<< *pig << endl;
		    pig++;
		  }
		abort();
		*/


	      }
	    else if (type == 3) Npauses++;
            else if (type == 4) Nresumes++;


	   } //loop over items in a evtfile
       evtfile.close();
       evtfile.clear(); // clear event status in case we had a bad file
	
    } //end loop over file subsections
      if(argc ==2)
	{
	  break;
	}
    } //end loop of run file numbers

  cout << endl;
  cout <<"live_trig = " << live_trig << endl;
  cout <<"raw_trig = " << raw_trig << endl;
  cout <<"Si_OR = " << Si_OR << endl;
  cout <<"CsI_OR = " << CsI_OR << endl;
  cout <<"XLM Complete = " << XLM_Complete << endl;
  cout << "Russian Pie OR = " << Rus_Pie_OR << endl;
  cout << "S2 Pie OR = " << S2_Pie_OR << endl;
  cout << "Russian CsI OR = " << Rus_CsI_OR << endl;
  cout << "S2 CsI OR = " << S2_CsI_OR << endl;
  cout << "FC = " << FC << endl;
  
  cout << '\n'<<"physics Events = " << physicsEvent << endl;
  cout << "Good physics Events = " << physicsEventGood << endl;
  
  if (physicsEvent > 0)cout << "bad/total = " << 
			 (1.-(double)physicsEventGood/(double)physicsEvent)*100. << 
			 " %"<< endl;
  
  cout << "physics Event Counters = " << physicsEventCounter << endl;
  cout << "scaler buffers = " << scalerBuffer << endl;
  cout << "Numbers of pauses = " << Npauses << endl;
  cout << "Number of resumes = " << Nresumes << endl;
  //         cout << "number of 2p+6Li = " << Det.N2p6Li << endl;
  //cout << "number of 8B(IAS) = " << Det.N_IAS << endl;
  
  cout << "Number of cross-talk R21-R22 = " << Det.Silicon->Ncross << endl;
  
  cout << "Number of 7Li -> t + a = " << Det.N7Li_ta << endl;


  Histo->write(); // this forces the histrograms to be read out to file
  
}
