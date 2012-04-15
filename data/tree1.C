#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

void tree1w()
{
   //create a Tree file tree1.root
   
   //create the file, the Tree and a few branches
   TFile f("gauss2_uniform_2d_train.root","recreate");
   TTree t1("t1","a simple Tree with simple variables");
   Float_t px, py, pz;
   Int_t classLabel;
   t1.Branch("px",&px,"px/F");
   t1.Branch("py",&py,"py/F");
   t1.Branch("classLabel",&classLabel,"classLabel/I");
   
   // signal
   for (Int_t i=0;i<10000;i++) {
     gRandom->Rannor(px,py);
     classLabel = 1;
     t1.Fill();
  }

   // background
   for (Int_t i=0;i<10000;i++) {
     px = (gRandom->Rndm()-0.5) * 10;
     py = (gRandom->Rndm()-0.5) * 10;
     classLabel = 0;
     t1.Fill();
  }
  
  //save the Tree header. The file will be automatically closed
  //when going out of the function scope
  t1.Write();
}

void tree1() {
   tree1w();
}
