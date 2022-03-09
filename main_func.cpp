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
#include <TStopwatch.h>

#include "tools.h"

using namespace std;

void main_func (){
    gROOT->SetBatch(kFALSE);
    vector <TH1F*> histos;
    TFile *outfile= new TFile("histograms/histograms_NA_700_run2.root", "RECREATE");
    list <string> names ={"pmt1_NA_e6_700_or_run2",
                        "pmt2_NA_e6_700_or_run2"
                        /*"coincidence_pmt1_ext_run2",
                        "coincidence_pmt2_ext_run2",
                        "pmt1_NA_e6_700_or",
                        "pmt2_NA_e6_700_or"
                        "coincidence_pmt1_ext",
                        "coincidence_pmt2_ext",
                        "coincidence_pmt1"
                        "coincidence_pmt2"
                        "pmt1_NA_e6_100_690",
                        "pmt2_NA_e6_100_851",
                        "pmt1_co_100_690", 
                        "pmt2_co_100_851",
                        "pmt1_na_100_690", 
                        "pmt2_na_100_851", 
                        "pmt1_cs_100_690", 
                        "pmt2_cs_100_851",
                        "pmt1_bkg_100_690", 
                        "pmt2_bkg_100_851", 
                        "pmt1_null", 
                        "pmt2_null"*/};
                        
    for(list<string>::const_iterator name = names.begin(); name != names.end(); ++name){
        TStopwatch time;
        time.Start();
        cout << "Processing: " << *name << endl; 
        histos=make_histo("data/" + *name + ".txt", *name);
        histos[0]->Write();
        histos[1]->Write();
        time.Stop();
        time.Print();
    }
    //getchar();
    outfile->Close();
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
