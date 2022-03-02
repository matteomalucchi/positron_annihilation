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

#include "tools.h"

using namespace std;

void main_func (){
    gROOT->SetBatch(kFALSE);
    vector <TH1F*> histos;

    //list <string> names ={"pmt1_co_100", "pmt2_co_100", "pmt1_na_100", "pmt2_na_100", "pmt1_cs_100", "pmt2_cs_100"};
    list <string> names ={"pmt2_cs_100"};
    for(list<string>::const_iterator i = names.begin(); i != names.end(); ++i){
        histos=make_histo("data/" + *i + ".txt", *i);
        TCanvas *c_charge = new TCanvas(&(*i + "_charge")[0] , &(*i + "_charge")[0]);
        histos[0]->Draw();
        TCanvas *c_amp = new TCanvas(&(*i + "_amp")[0], &(*i + "_amp")[0]);
        histos[1]->Draw();

        
        c_charge->SaveAs(&("plots/" + *i + "_charge.png")[0]);
        c_amp->SaveAs(&("plots/" + *i + "_amp.png")[0]);
    }
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

// $ g++ -Wall -Wextra -Werror -pedantic -std=c++14 main_func.cpp tools.cpp  `root-config --glibs --cflags` -o someExecutable
// ./ someExecutable
