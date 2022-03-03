#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <stdio.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>


#include "tools.h"

using namespace std;

void main_func (){
    gROOT->SetBatch(kFALSE);
    vector <TH1F*> histos;
    //TFile f("histograms.root", "RECREATE");
    list <string> names ={/*"pmt1_co_100", 
                        "pmt2_co_100",
                        "pmt1_na_100", 
                        "pmt2_na_100", 
                        "pmt1_cs_100", 
                        "pmt2_cs_100", 
                        "pmt1_bkg_100", 
                        "pmt2_bkg_100", 
                        "pmt1_null", 
                        "pmt2_null", 
                        "wave0_co_1", 
                        "wave0_co_2", */
                        "wave0_na_1"};
                        
    for(list<string>::const_iterator i = names.begin(); i != names.end(); ++i){
        histos=make_histo("data/" + *i + ".txt", *i);
        TCanvas *c_charge = new TCanvas(&(*i + "_charge")[0] , &(*i + "_charge")[0]);
        histos[0]->Draw();
        TCanvas *c_amp = new TCanvas(&(*i + "_amp")[0], &(*i + "_amp")[0]);
        histos[1]->Draw();

        /*histos[0]->SetName(&(*i + "_charge")[0]);
        histos[1]->SetName(&(*i + "_amp")[0]);
        histos[0]->Write();
        histos[1]->Write();*/
        c_charge->SaveAs(&("plots/" + *i + "_charge.png")[0]);
        c_amp->SaveAs(&("plots/" + *i + "_amp.png")[0]);

    }
    //f.Close();
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
// ./ someExecutable
