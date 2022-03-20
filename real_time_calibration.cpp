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


using namespace std;



vector <TH1F*> make_histo(string name, float charge_min, float charge_max, float y_rescale, float peak){
    ifstream myfile;
    myfile.open("data/" + name + ".txt", ios::in | ios::out);  
    name= name + "_clear";
    vector <float> t(1030);
    vector <vector<float>> v;
    vector <TH1F*> histos;
    iota(begin(t), end(t), 0);


    if (myfile.is_open()){
        string tp;
        vector <float> w;
        int k=0, m, j=0;
        while(getline(myfile, tp)){ //read data from file object and put it into string.
            if (isalpha(tp[0]) == 0) {
                if (k>0) {
                    k=0;
                    w.clear();
                }
                m = stof(tp);
                w.push_back(m);
                j++;
            }
            else {
                if (j>0) {
                    j=0;
                    v.push_back(w);
                }                
                k++;
            }
        }
        //non prendo l'ultimo evento perchè è incompleto
        myfile.close(); //close the file object.
    }    
    
    vector <float> charge (v.size());
    vector <float> amp (v.size());
    float h, u, min;
    vector <float> idx_strange;
    for (long unsigned int i=0; i<v.size(); i++){
        h = 0;
        u=0;
        for (int j=0; j<100; j++){
            h += v[i][j] / 100;
            if (v[i][j]<14600 && u==0){
                idx_strange.push_back(i);
                u++;
            }
        }
        for (int j=0; j<1000; j++){
            charge[i] += h-v[i][j];
        }
        min =*min_element(v[i].begin(), v[i].end());
        amp[i]=h-min;
    }
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 1300, 0, 250000);
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 1400, 0, 3500);
    for (long unsigned int i=0; i< charge.size(); i++){
        if ((find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()) && charge[i]> charge_min && charge[i]<charge_max){
            histo_charge->Fill(charge[i]+peak/100);
            histo_amp->Fill(amp[i]);
        }
    }
    histo_charge->Scale(1/y_rescale);

    histo_charge->SetName(&(name + "_charge")[0]);
    histo_amp->SetName(&(name + "_amp")[0]);
    histos.push_back(histo_charge);
    histos.push_back(histo_amp);

    TCanvas *c_charge = new TCanvas(&(name + "_charge")[0] , &(name + "_charge")[0]);
    histo_charge->Draw();
    c_charge->SaveAs(&("real_time_calibration/" + name + "_charge.png")[0]);
/*
    TCanvas *c_amp = new TCanvas(&(name + "_amp")[0], &(name + "_amp")[0]);
    histo_amp->Draw();
    c_amp->SaveAs(&("real_time_calibration/" + name + "_amp.png")[0]);*/

    //cout << "numero di idx_strange  "<<idx_strange.size() <<endl;

    return histos;
}

vector<float> fit_gaus (TH1F* histo_sub, vector<float> ranges, string name, string type){

    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[0], ranges[1]);
    TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[2], ranges[3]);
    TF1 *gaus3 = new TF1("gaus3","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[7], ranges[8]);
    TF1 *gaus4 = new TF1("gaus4","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[9], ranges[10]);

    gaus1->SetParameters(700, abs(ranges[0]+ranges[1])/2, abs(ranges[0]-ranges[1]));
    gaus2->SetParameters(1000, abs(ranges[2]+ranges[3])/2, abs(ranges[2]-ranges[3]));
    gaus3->SetParameters(100, abs(ranges[7]+ranges[8])/2, abs(ranges[7]-ranges[8]));
    gaus4->SetParameters(100, abs(ranges[9]+ranges[10])/2, abs(ranges[9]-ranges[10]));

    TCanvas *c_sub = new TCanvas(&(name_b + type +"_sub")[0] ,&(name_b + type +"_sub")[0]);

    histo_sub->Fit("gaus1","R", "SAME");
    gaus1->Draw("SAME");
    gStyle->SetOptFit(1111);
    histo_sub->Fit("gaus2","R", "SAME"); 
    gaus2->Draw("SAME");
    gStyle->SetOptFit(1111);
    histo_sub->Fit("gaus3","R", "SAME"); 
    gaus3->Draw("SAME");
    gStyle->SetOptFit(1111);
    histo_sub->Fit("gaus4","R", "SAME"); 
    gaus4->Draw("SAME");
    gStyle->SetOptFit(1111);
    histo_sub->Draw("SAME");

    c_sub->SaveAs(&("real_time_calibration/"+name_b + type +"_sub.png")[0]);
    c_sub->Write();

    vector<float> gaus_params;

    gaus_params.push_back(gaus1->GetParameter(1));
    gaus_params.push_back(gaus1->GetParError(1));
    gaus_params.push_back(gaus2->GetParameter(1));
    gaus_params.push_back(gaus2->GetParError(1));        
    gaus_params.push_back(gaus3->GetParameter(1));
    gaus_params.push_back(gaus3->GetParError(1));
    gaus_params.push_back(gaus4->GetParameter(1));
    gaus_params.push_back(gaus4->GetParError(1));

    return gaus_params;
}

void real_time_calibration(){
    TFile *f = new TFile("histograms/histograms_RealTimeCalibration.root");
    TFile *outfile= new TFile("real_time_calibration/histograms_RealTimeCalibration_clear.root", "RECREATE");

    // primi due sono picco NA, poi picco cs, poi range da sottrare del cs e poi la posizione del picco del cs da sottrare
    // poi il range del primo picco co e poi del secondo picco co
    map <vector<string>, vector<float>> pair_names ={ 
        {{"pmt1_NA_cs_e6_500_run2", "pmt1_NA_cs_co_e6_500_run2"}, 
            {78000, 85000, 100000, 108000, 150000, 210000, 195000, 176000, 187000, 200000, 212000}}

        /*{"pmt1_NA+cs_e6_100_run1", 
                        "pmt2_NA+cs_e6_100_run1",
                        "pmt3_NA+cs_e6_30_run1",
                        "pmt1_NA+cs+co_e6_100_run1",
                        "pmt2_NA+cs+co_e6_100_run1",
                        "pmt3_NA+cs+co_e6_30_run1", */

                        /*"pmt2_NA+cs_e6_500_run2",
                        "pmt3_NA+cs_e6_100_run2",
                        "pmt2_NA+cs+co_e6_500_run2",
                        "pmt3_NA+cs+co_e6_100_run2"*/};
    
    vector <TH1F*> histos;

    string type= "_charge";
    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first[0];
        const auto name_b = pair_name.first[1];
        const auto ranges = pair_name.second;
        
        TH1F *histo_a = nullptr;
        f->GetObject(&(name_a + type)[0], histo_a);
        TH1F *histo_b = nullptr;
        f->GetObject(&(name_b + type)[0], histo_b);

        float int_NA_a=histo_a->Integral(histo_a->FindFixBin(ranges[0]),histo_a->FindFixBin(ranges[1]));
        float int_NA_b=histo_b->Integral(histo_b->FindFixBin(ranges[0]),histo_b->FindFixBin(ranges[1]));
        cout<< int_NA_a/int_NA_b << endl;

        float int_cs_a=histo_a->Integral(histo_a->FindFixBin(ranges[2]),histo_a->FindFixBin(ranges[3]));
        float int_cs_b=histo_b->Integral(histo_b->FindFixBin(ranges[2]),histo_b->FindFixBin(ranges[3]));
        cout<< int_cs_a/int_cs_b << endl;

        float y_rescale = (int_NA_a/int_NA_b + int_cs_a/int_cs_b )/2;

        cout << "Processing: " << name_a << endl; 
        histos=make_histo(name_a, ranges[4], ranges[5], y_rescale, ranges[6]);

        TH1F *histo_sub = (TH1F*)histo_b->Clone("histo_sub");
        histo_sub->Add(histos[0], -1);
        TCanvas *c_all = new TCanvas(&(name_b + type +"_all")[0] ,&(name_b + type +"_all")[0]);
        histo_sub->SetName(&(name_b + type +"_all")[0]);
        histo_b->SetLineColor(kRed);
        histo_sub->SetLineColor(kBlue);
        histos[0]->SetLineColor(kGreen);

        histo_sub->Draw();
        histo_b->Draw("SAME");
        histos[0]->Draw("SAME");
        c_all->SaveAs(&("real_time_calibration/"+name_b + type +"_all.png")[0]);



        fit_gaus(histo_sub, ranges, name_b, type);
        





        histo_a->Write();
        histo_b->Write();
        histos[0]->Write();
        //histos[1]->Write();
        histo_sub->Write();
        c_all->Write();



    }

    getchar();
    outfile->Close();

}