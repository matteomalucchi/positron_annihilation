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

/*    map<string, vector<float>> runs{
        {"run3",{-1,6, -1, 12, 7, 30,0,15,-1,0,8.28, -0.51}},
        {"run9",{-1,5, -1, 12, 5, 30,0,15,-1,0,7.45, -0.67}},
        {"run10",{-1,5, -1, 12, 5.5, 30,0,15,-1,0, 6.35, -0.52}},
    };
*/
    map<string, vector<float>> runs{
        {"run3",{-5,2.5, -5, 8, 4.5, 15,0,15,-1,0,8.28, -0.51}},
        {"run9",{-5,2.5, -5, 8, 2.5, 15,0,15,-1,0,7.45, -0.67}},
        {"run10",{-5,2.5, -5, 8, 3, 15,0,15,-1,0, 6.35, -0.52}},
    };
    TNtuple *ntuple;
    int bins=400, min=-10,max=20;

    for (const auto &run : runs) {
        const auto ntuple_name = run.first;
        const auto ranges = run.second;  
        cout << ntuple_name<< endl<<endl;
        f->GetObject(&(ntuple_name)[0], ntuple);
        TCanvas *c = new TCanvas(&("histo_"+ntuple_name)[0], &("histo_"+ntuple_name)[0]);
        TH1F* h= new TH1F(&("histo_"+ntuple_name)[0], &("histo_"+ntuple_name)[0], bins, min, max);

        for (auto ievt: ROOT::TSeqI(ntuple->GetEntries())){// Loop over events
            ntuple->GetEntry(ievt);
            auto catm = ntuple->GetArgs(); // Get a row
            if (abs(catm[2]-catm[5])<13 /*&& catm[6]>0.6*/ && catm[0]<0.6 && catm[3]<0.6 && catm[9]==1 && catm[7]>0.05){
                h->Fill((catm[2]+catm[5])/2-catm[8]);
            }
        }

        TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[0], ranges[1]);
        gaus1->SetParameters(100, 5,2);
        gaus1->SetLineColor(kRed);
        h->Fit("gaus1","R");

        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[2], ranges[3]);
        gaus2->FixParameter(0,gaus1->GetParameter(0));
        gaus2->FixParameter(1,gaus1->GetParameter(1));
        gaus2->FixParameter(2,gaus1->GetParameter(2));
        gaus2->SetLineColor(kYellow);
/*
        TF1 *expo1 = new TF1("expo1", "exp([0]+[1]*x)",ranges[4], ranges[5]);
        expo1->SetParameters(ranges[10], ranges[11]);
        expo1->SetLineColor(kGreen);
        h->Fit("expo1","R+");
        

        Double_t par[5];
        gaus1->GetParameters(&par[0]);
        expo1->GetParameters(&par[3]);
        */

        TF1 *gaus_exp = new TF1("gaus_exp","exp([0]+[1]*x)+[2]*exp(-0.5*pow(((x-[3])/[4]),2))",ranges[4], ranges[5]);
        gaus_exp->SetParameter(0,ranges[10]);
        gaus_exp->SetParameter(1,ranges[11]);
        
        gaus_exp->FixParameter(2,gaus1->GetParameter(0));
        gaus_exp->FixParameter(3,gaus1->GetParameter(1));
        gaus_exp->FixParameter(4,gaus1->GetParameter(2));
        gaus_exp->SetLineColor(kGreen);

        //TF1 *gaus_exp = new TF1("gaus_exp","[0]*exp(-0.5*pow(((x-[1])/[2]),2)) + exp([3]+[4]*x)",ranges[2], ranges[3]);
        //gaus_exp->SetParameters(100,5,2,ranges[10], ranges[11]);

        h->Fit("gaus_exp","R");

        TF1 *expo1 = new TF1("expo1","exp([0]+[1]*x)",ranges[4], ranges[5]);
        expo1->SetParameter(0,ranges[10]);
        expo1->SetParameter(1,ranges[11]);    
        expo1->SetLineColor(kBlack);
        h->Fit("expo1","R");

        TF1 *expo2 = new TF1("expo2","exp([0]+[1]*x)",ranges[4], ranges[5]);
        expo2->SetParameter(0,gaus_exp->GetParameter(0));
        expo2->SetParameter(1,gaus_exp->GetParameter(1));    
        expo2->SetLineColor(kGray);

        h->Draw();
        gaus1->Draw("SAME");
        gaus_exp->Draw("SAME");
        expo1->Draw("SAME");
        expo2->Draw("SAME");
        gaus2->Draw("same");

        c->BuildLegend();

        cout<<endl;
        cout << "mean life "<< ntuple_name<<"= "  <<1/abs(gaus_exp->GetParameter(1)) <<"+-"<< abs(gaus_exp->GetParError(1)/(gaus_exp->GetParameter(1)*gaus_exp->GetParameter(1)))<<endl;


/*
        gaus1->SetLineColor(kRed);
        gStyle->SetOptFit(1111);

        h->Fit("gaus1","R", "SAME");
        gaus1->SetLineColor(kRed);
        gStyle->SetOptFit(1111);

        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[2], ranges[3]);
        gaus2->FixParameter(0,gaus1->GetParameter(0));
        gaus2->FixParameter(1,gaus1->GetParameter(1));
        gaus2->FixParameter(2,gaus1->GetParameter(2));
        gaus2->SetLineColor(kGreen);
                
        h->Draw();
        gaus1->Draw("SAME");
        gaus2->Draw("SAME");

        TCanvas *c_sub = new TCanvas(&("histo_sub_"+ntuple_name)[0], &("histo_sub_"+ntuple_name)[0]);   
        TH1F *h_sub = (TH1F*)h->Clone(&("histo_sub"+ntuple_name)[0]);
        h_sub->Reset();
        for (int j=0;j<bins;j++){
            h_sub->SetBinContent(j,h->GetBinContent(j)-gaus2->Eval(h->GetXaxis()->GetBinCenter(j)));
        }
        TF1 *expo1 = new TF1("expo1", "expo(0)",ranges[4], ranges[5]);
        expo1->SetParameters(ranges[10], ranges[11]);
        expo1->SetParLimits(0,ranges[6], ranges[7]);
        expo1->SetParLimits(1, ranges[8], ranges[9]) ;
        
        h_sub->Fit("expo1","R");
        expo1->SetLineColor(kRed);
        gStyle->SetOptFit(1111);
        h_sub->Draw();
        expo1->Draw("same");
        cout<<endl;
        cout << "mean life "<< ntuple_name<<"= "  <<1/abs(expo1->GetParameter(1)) <<endl;
*/
    }
}

