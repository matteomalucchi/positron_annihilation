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
#include<TMath.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TStopwatch.h>
#include<TNtuple.h>


using namespace std;

string type_of_file = "triple_quad";

void chi_square_times(){
    TFile *f= new TFile(&("triple/ntuple_" + type_of_file + ".root")[0]);

    float err=1.8;
    string run = "run8";
    TNtuple *ntuple;
    int bins=400, min=-10,max=20;
    int num=3;
    float t0=0,weights_sum=0, t0_err, chi_square=0;
    Double_t p_value;
    f->GetObject(&(run)[0], ntuple);
    vector<float> times(num), weights(num);
    for (auto ievt: ROOT::TSeqI(ntuple->GetEntries())){// Loop over events
        ntuple->GetEntry(ievt);
        auto catm = ntuple->GetArgs(); // Get a row
        if (abs(catm[2]-catm[5])<13 && catm[6]<0.45 && catm[0]<0.45 && catm[3]<0.45 && catm[9]==1 && catm[7]>0.05 && abs(catm[2]/2+catm[5]/2-catm[8])<23 ){
            times[0]=catm[2];
            times[1]=catm[5];
            times[2]=catm[8];

            // weighted average
            for(int k=0; k<num;k++){
                weights[k]=1/(err*err);
                weights_sum+=weights[k];
            }

            for(int k=0; k<num;k++){
                t0+=times[k]*weights[k]/weights_sum;
            }
            t0_err=sqrt(1/weights_sum);
            cout << t0 <<"  "<< t0_err <<endl;
            for(int k=0; k<num;k++){
                chi_square+=pow((times[k]-t0)/err,2);
            }
            cout <<chi_square<<endl;

            p_value=TMath::Prob(chi_square, 2);
            cout <<p_value <<endl;
            cout<<endl;

            t0=0;
            chi_square=0;
            weights_sum=0;


        }
    }



}