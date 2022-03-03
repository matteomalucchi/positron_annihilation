#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <stdio.h>
#include <cstdio>

#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>

#include "tools.h"

using namespace std;

void main_func (){
    gROOT->SetBatch(kFALSE);
    vector <TH1F*> histos;
    TH1F* histo_temp;
    TFile *f= new TFile("histograms.root", "RECREATE");
    list <string> names ={"pmt1_co_100", 
                        "pmt2_co_100",
                        "pmt1_na_100", 
                        "pmt2_na_100", 
                        "pmt1_cs_100", 
                        "pmt2_cs_100", 
                        "pmt1_bkg_100", 
                        "pmt2_bkg_100", 
                        "pmt1_null", 
                        "pmt2_null", 
                        /*"wave0_co_1", 
                        "wave0_co_2",
                        "wave0_na_1"*/};
                        
    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        cout << "Processing: " << *name << endl; 
        histos=make_histo("data/" + *name + ".txt", *name);

        histo_temp=histos[0];
        histo_temp->Write();
        histo_temp=histos[1];
        histo_temp->Write();
    }
    getchar();
    f->Close();
}

int main() {
    main_func();
}



// execute with:
// $ root -l
// $ .L tools.cpp
// $ .x main_func.cpp


// execute with:
// $ root -l
// $.L tools.cpp+
// $ root -l tools_cpp.so main_func.cpp

// $ g++ -Wall -Wextra -Werror -pedantic -std=c++14 main_func.cpp tools.cpp  `root-config --glibs --cflags` -o main_func
// ./main_func
