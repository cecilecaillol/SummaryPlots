#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "math.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "string.h"

using namespace std;

int main(int argc, char** argv) {
  string dir=*(argv + 1);
  string name=*(argv + 2);
  float m_min, m_max, step;
  if (argc > 1) {
        m_min = atof(argv[3]);
	m_max = atof(argv[4]);
        step = atof(argv[5]);
  }
  ofstream myfile;
  myfile.open (TString::Format("%s.txt",dir.c_str()));
  for(unsigned int i=0; i<((m_max-m_min)/step)+1; ++i){
    float mass = m_min+i*step;
    double value=-1.;
    TString fullpath;
    if(fmod(mass,1) == 0)
       fullpath=(TString::Format("%s/%d/%s%d.root", dir.c_str(), (int) mass, name.c_str(), (int) mass));
    else
       fullpath=(TString::Format("%s/%.1f/%s%.1f.root", dir.c_str(), mass, name.c_str(), mass));
    TFile* file = new TFile(fullpath);
    if(file->IsZombie()){
      std::cout << "INFO: file not found: " << fullpath  << std::endl; 
    }
    else{
      TTree* limit = (TTree*) file->Get("limit");
      double x; float y;
      limit->SetBranchAddress("limit", &x);
      limit->SetBranchAddress("quantileExpected", &y);
      int nevent = limit->GetEntries();
      for(int i=0; i<nevent; ++i){
          limit->GetEvent(i);
          if(y==-1)
            value = x;
      }
    }
    if (dir=="mmtt") value=value/100;
    cout<<mass<<" "<<value<<endl;
    myfile<<mass<<" "<<value<<endl;
    file->Close();
  }
  myfile.close();
  return 1;
}

