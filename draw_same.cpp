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
        //{"pmt1_NA_e6_ext_run1_charge","pmt1_NA_e6_ext_run2_charge"},
        //{"pmt2_NA_e6_ext_run1_charge","pmt2_NA_e6_ext_run2_charge"},
        {"pmt1_NA_e6_100_run1_charge", "pmt1_NA_e6_100_run2_charge"},
        {"pmt2_NA_e6_100_run1_charge", "pmt2_NA_e6_100_run2_charge"}
    };

    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first;
        const auto name_b = pair_name.second;
/*
        string name = name_a.substr(5, name_a.size()-1);
        string run_a = name_a.substr(14, 18);
        string run_b = name_b.substr(14, 18);
*/
        TFile *f = new TFile("histograms/histograms_NA_run12.root");

        TH1F *histo_a = nullptr;
        TH1F *histo_b = nullptr;

        f->GetObject(&(name_a)[0], histo_a);
        f->GetObject(&(name_b)[0], histo_b);

        TCanvas *c= new TCanvas(&(name_a)[0], &(name_a)[0]);
        histo_a->SetLineColor(1);
        histo_b->SetLineColor(4);
        histo_a->Draw();
        histo_b->Draw("SAME");
    }


}