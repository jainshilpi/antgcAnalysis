#include <string>
//#include "TChain.h
#if defined(__CINT__) && !defined(__MAKECINT__)
#include "Analyse.C+"
#else
#include "Analyse.C"
#endif

#include <fstream>

void runAll()
{

  cout<<"inside runAll.C"<<endl;

  //string idType = "90"; ///70 percent, 80, 90 or 95 signal eff
  
  
  double lumi = 41*1000;  //pb-1
  cout<<"now making class object "<<endl;

  Analyse t;
  
  cout<<"done making class object "<<endl;
  
  
  cout<<"====Now running over samples===="<<endl;

  ifstream is("files.list");
  //ifstream is("files_data.list");
  string str;
  while(getline(is, str))
    {

      string fname, fnameout, idType;
      int isData, itype;
      double xsec;
      
      is >> fname >> fnameout >> isData >> xsec >> itype >> idType;
      int fstr = fname.find("#",0);
      if(fstr!=string::npos)
	{
	  continue;
	}
      
      cout<<str<<endl;
      
      if(!isData)
	xsec = xsec*lumi;
      
      t.Loop(fname, fnameout, isData, xsec, itype,idType);      
    }
  


  
}//void runAll()
