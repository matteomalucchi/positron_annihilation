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
#include<TNtuple.h>


using namespace std;

string type_of_file = "triple_lin";

void mean_life(){
    TFile *f= new TFile(&("triple/ntuple_" + type_of_file + ".root")[0]);
    vector<string> runs={"run3", "run9", "run10"};
    TNtuple *ntuple;
    int bins=500, min=-20,max=20;

    for(int i=0;i<runs.size();i++){
        f->GetObject("run3", ntuple);
        TCanvas *c = new TCanvas(&("histo_"+runs[i])[0], &("histo_"+runs[i])[0]);
        TH1F* h= new TH1F(&("histo_"+runs[i])[0], &("histo_"+runs[i])[0], bins, min, max);

        for (auto ievt: ROOT::TSeqI(ntuple->GetEntries())){// Loop over events
            ntuple->GetEntry(ievt);
            auto catm = ntuple->GetArgs(); // Get a row
            if (abs(catm[2]-catm[5])<13 && catm[6]>0.6 ){
                h->Fill(catm[2]-catm[8]);
            }
        }
        h->Draw();
        TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",-1,6.5);
        gaus1->SetParameters(1000, 5,5);
        h->Fit("gaus1","R", "SAME");
        gaus1->SetLineColor(kRed);
        gaus1->Draw("SAME");
        gStyle->SetOptFit(1111);

        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",-1,12);
        gaus2->FixParameter(0,gaus1->GetParameter(0));
        gaus2->FixParameter(1,gaus1->GetParameter(1));
        gaus2->FixParameter(2,gaus1->GetParameter(2));
        gaus2->SetLineColor(kGreen);
        gaus2->Draw("SAME");

        
        TH1F *h_sub = (TH1F*)h->Clone(&("histo_sub"+runs[i])[0]);
        for (int j=0;j<bins;j++){
            h_sub->SetBinContent(j,h->GetBinContent(j)-gaus2->Eval(h->GetXaxis()->GetBinCenter(j)));
        }
        TCanvas *c_sub = new TCanvas(&("histo_sub_"+runs[i])[0], &("histo_sub_"+runs[i])[0]);
        h_sub->Draw();

    }
}

