#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
    
#include "make_histo.h"

using namespace std;

void plot (){

    TH1F *histo_co_1=  new TH1F("isto_co_1","co_1", 900, 200000, 500000);
    TCanvas *c_co_1 = new TCanvas("c_co_1","histo charge");
    histo_co_1 = make_histo("data/wave0_co_1.txt", histo_co_1);
    histo_co_1->Draw();

    TH1F *histo_co_2=  new TH1F("isto_co_2","co_2", 900, 200000, 500000);
    TCanvas *c_co_2 = new TCanvas("c_co_2","histo charge");
    histo_co_2 = make_histo("data/wave0_co_2.txt", histo_co_2);
    histo_co_2->Draw();

    TH1F *histo_na_1=  new TH1F("isto_na_1","na_1", 900, 200000, 500000);
    TCanvas *c_na_1 = new TCanvas("c_na_1","histo charge");
    histo_na_1 = make_histo("data/wave0_na_1.txt", histo_na_1);
    histo_na_1->Draw();
    
}

