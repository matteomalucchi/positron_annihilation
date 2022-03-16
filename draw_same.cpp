#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm> 
#include <iterator>
#include <list>
#include <float.h>
#include <tuple>
//#include<bits/stdc++.h>

#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>

using namespace std;

void draw_same(){

    map <string, string> pair_names ={
        {"pmt1_NA_e6_ext_run1_charge","pmt1_NA_e6_ext_run2_charge"},
        {"pmt2_NA_e6_ext_run1_charge","pmt2_NA_e6_ext_run2_charge"},
        {"pmt1_NA_e6_100_run1_charge", "pmt1_NA_e6_100_run2_charge"},
        {"pmt2_NA_e6_100_run1_charge", "pmt2_NA_e6_100_run2_charge"},
    };

    TFile *f = new TFile("histograms/histograms_calibration.root");
    TH1F *histo = nullptr
    f->GetObject(&(*sample + type)[0], histo);


}