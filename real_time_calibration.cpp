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
    //TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 1400, 0, 3500);
    for (long unsigned int i=0; i< charge.size(); i++){
        if ((find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()) && charge[i]> charge_min && charge[i]<charge_max){
            histo_charge->Fill(charge[i]+peak/100);
            //histo_amp->Fill(amp[i]);
        }
    }
    histo_charge->Scale(1/y_rescale);

    histo_charge->SetName(&(name + "_charge")[0]);
    //histo_amp->SetName(&(name + "_amp")[0]);
    histos.push_back(histo_charge);
    //histos.push_back(histo_amp);

    TCanvas *c_charge = new TCanvas(&(name + "_charge")[0] , &(name + "_charge")[0]);
    histo_charge->SetLineColor(3);
    histo_charge->Draw();
    c_charge->SaveAs(&("real_time_calibration/" + name + "_charge.png")[0]);
/*
    TCanvas *c_amp = new TCanvas(&(name + "_amp")[0], &(name + "_amp")[0]);
    histo_amp->Draw();
    c_amp->SaveAs(&("real_time_calibration/" + name + "_amp.png")[0]);*/

    return histos;
}

vector<float> fit_gaus_NAb (TH1F* histo, vector<float> ranges, string name, string type){
    
    cout << endl;
    cout << "__________________________ Gaussian fit: " << name << type << " __________________________"<<endl; 
    cout << endl;

    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[6], ranges[7]);

    gaus1->SetParameters(700, abs(ranges[6]+ranges[7])/2, abs(ranges[6]-ranges[7]));


    TCanvas *c_NAb = new TCanvas(&(name + type +"_NAb_fit")[0] ,&(name + type +"_NAb_fit")[0]);

    histo->Fit("gaus1","R", "SAME");
    gaus1->Draw("SAME");
    gStyle->SetOptFit(1111);
    histo->Draw("SAME");

    c_NAb->SaveAs(&("real_time_calibration/"+name + type +"_NAb_fit.png")[0]);
    c_NAb->Write();

    vector<float> gaus_params;

    gaus_params.push_back(gaus1->GetParameter(1));
    gaus_params.push_back(gaus1->GetParError(1));

    return gaus_params;
}

vector<float> fit_gaus (TH1F* histo_sub, vector<float> ranges, string name, string type){
    
    cout << endl;
    cout << "__________________________ Gaussian fit: " << name << type << " __________________________"<<endl; 
    cout << endl;

    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[0], ranges[1]);
    TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[2], ranges[3]);
    TF1 *gaus3 = new TF1("gaus3","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[8], ranges[9]);
    TF1 *gaus4 = new TF1("gaus4","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges[10], ranges[11]);

    gaus1->SetParameters(700, abs(ranges[0]+ranges[1])/2, abs(ranges[0]-ranges[1]));
    gaus2->SetParameters(1000, abs(ranges[2]+ranges[3])/2, abs(ranges[2]-ranges[3]));
    gaus3->SetParameters(100, abs(ranges[8]+ranges[9])/2, abs(ranges[8]-ranges[9]));
    gaus4->SetParameters(100, abs(ranges[10]+ranges[11])/2, abs(ranges[10]-ranges[11]));

    TCanvas *c_sub = new TCanvas(&(name + type +"_sub_fit")[0] ,&(name + type +"_sub_fit")[0]);

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

    c_sub->SaveAs(&("real_time_calibration/"+name + type +"_sub_fit.png")[0]);
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

auto fit_lin(string name_final, string type, vector<float> gaus_params, TFile * outfile, vector<float> gaus_param_NAb){
    float x[4]={0.66, 1.17, 1.27, 1.33};
    float y[sizeof(x)/sizeof(x[0])]={gaus_params[2], gaus_params[4],gaus_param_NAb[0],  gaus_params[6]};
    //float y_err[sizeof(x)/sizeof(x[0])] = {gaus_params[2]/100, gaus_params[4]/100, gaus_params[6]/100};
    float y_err[sizeof(x)/sizeof(x[0])] = {gaus_params[3], gaus_params[5],gaus_param_NAb[1], gaus_params[7]};

    cout << endl;
    cout << "__________________________ Linear fit: " << name_final << type << " __________________________"<<endl; 
    cout << endl;
    gStyle->SetOptStat(0);

    TCanvas *c_fit_lin = new TCanvas(&(name_final + type + "_fit_lin")[0], &( name_final + type + "_fit_lin")[0]);
    TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.00001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    TGraphErrors* gr = new TGraphErrors(sizeof(x)/sizeof(x[0]),x,y,nullptr,y_err);
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.5);
    gr->SetMarkerStyle(1);
    gr->Draw("APE");
    TF1 *   linear  = new TF1("linear","[0]+[1]*x", 0, 1.4);
    linear->SetParNames ("Off-set","Calibration factor");
    gr->Fit("linear");
    linear->Draw("SAME");
    gStyle->SetOptFit(111);

    pad2->cd();
    float diff[sizeof(x)/sizeof(x[0])];
    float diff_norm[sizeof(x)/sizeof(x[0])];
    for (Int_t i=0;i<sizeof(x)/sizeof(x[0]);i++) {
        diff[i] =y[i]-linear->Eval(x[i]);
        diff_norm[i]=diff[i]/y_err[i];
    }  
    float one_array[sizeof(x)/sizeof(x[0])]={1,1,1,1};
    
    TGraphErrors* gr2 = new TGraphErrors(sizeof(x)/sizeof(x[0]),x,diff_norm,nullptr,one_array);
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.5);
    gr2->SetMarkerStyle(1);
    TF1 *   zero  = new TF1("zero","0*x", 0, 1.4);
    gr2->Draw("APE");
    zero->Draw("SAME");
    c_fit_lin->cd();
    c_fit_lin->SaveAs(&("real_time_calibration/" + name_final + type + "_fit_lin.png")[0]);
    outfile->cd();
    c_fit_lin->Write();

    vector<double> off_set{linear->GetParameter(0),linear->GetParError(0)};
    vector<double> cal_factor{linear->GetParameter(1),linear->GetParError(1)};
    vector<vector<double>> lin_params{off_set, cal_factor};

    return lin_params;
}

auto division(vector<double> off_set, vector<double> cal_factor, float peak, float peak_err){
    vector<float> p;
    p.push_back((peak-off_set[0])/cal_factor[0]);
    float e = sqrt(pow(peak_err/cal_factor[0], 2) +pow(off_set[1]/cal_factor[0], 2) +pow((peak-off_set[0])*cal_factor[1]/(cal_factor[0]*cal_factor[0]), 2));
    p.push_back(e);
    return p;
}

auto peak_energy(string name_final, string type, vector<vector<double>> lin_params, vector<float> gaus_params){
    vector<float> energy;
    vector<float> x(division(lin_params[0], lin_params[1], gaus_params[0], gaus_params[1]/*gaus_params[0]/100)*/));    
    //vector<float> y(division(lin_params[0], lin_params[1], gaus_params[10], gaus_params[11]));
    energy.insert(energy.end(), x.begin(), x.end());
    //energy.insert(energy.end(), y.begin(), y.end());

    string peak1 = "energy peak " + name_final + " " + type + " = " + to_string(energy[0]) + "+-" + to_string(energy[1]) + "\n";
    //string peak2 = "energy peak " + name_final + " " + type + " = " + to_string(energy[2]) + "+-" + to_string(energy[3])+ "\n";
    vector <string> energy_string;
    energy_string.push_back(peak1);
    //energy_string.push_back(peak2);

    return energy_string;

}


void real_time_calibration(){
    TFile *f = new TFile("histograms/histograms_RealTimeCalibration_rebin.root");
    TFile *outfile= new TFile("real_time_calibration/plots_RealTimeCalibration_fit_low.root", "RECREATE");
    ofstream out_file("real_time_calibration/peak_energy_low.txt");

    // primi due sono picco NA| poi picco cs| poi range da sottrare del cs |e poi il range del picco b del na da sottrarre|
    // poi il range del primo picco co |e poi del secondo picco co
    map <vector<string>, vector<float>> pair_names ={ 
        {{"pmt1_NA_cs_e6_100_run1", "pmt1_NA_cs_co_e6_100_run1"}, 
            {78000-2000, 85000-2000, 100000-2000, 108000-2000, 150000-2000, 210000-2000, 190000-2000, 204000-2000, 176000-2000, 187000-2000, 200000-2000, 212000-2000}},   
        /*{{"pmt1_NA_cs_e6_500_run2", "pmt1_NA_cs_co_e6_500_run2"}, 
            {78000, 85000, 100000, 108000, 150000, 210000, 190000, 204000, 176000, 187000, 200000, 212000}}, 

        {{"pmt2_NA_cs_e6_100_run1", "pmt2_NA_cs_co_e6_100_run1"}, 
            {69000, 76000, 90000, 108000-9000, 150000-9000, 210000-9000, 172000, 182000, 157000, 167000, 182000, 190000}},
        {{"pmt2_NA_cs_e6_500_run2", "pmt2_NA_cs_co_e6_500_run2"}, 
            {70000, 76000, 89000, 98000, 150000, 190000, 172000, 182000, 160000, 168000, 182000, 190000}},

        {{"pmt3_NA_cs_e6_30_run1", "pmt3_NA_cs_co_e6_30_run1"}, 
            {13000, 14700, 16700, 18800, 150000-9000, 210000-9000, 177000, 157000, 167000, 182000, 190000}},
        {{"pmt3_NA_cs_e6_100_run2", "pmt3_NA_cs_co_e6_100_run2"}, 
            {70000, 76000, 89000, 98000, 150000, 190000, 177000, 160000, 168000, 182000, 190000}}*/
        };
    
    vector <string> energy;

    string type= "_charge";

    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first[0];
        const auto name_b = pair_name.first[1];
        const auto ranges = pair_name.second;
        vector <TH1F*> histos;
        
        string pmt = name_a.substr(0, 4);
        string final = name_a.substr(10, name_a.size()-1);
        string name_final = pmt+final;

        cout << "Processing: " << name_final << endl; 

        TH1F *histo_a = nullptr;
        f->GetObject(&(name_a + type)[0], histo_a);
        TH1F *histo_b = nullptr;
        f->GetObject(&(name_b + type)[0], histo_b);
        TCanvas *c_a = new TCanvas(&(name_a)[0] ,&(name_a)[0]);
        histo_a->SetLineColor(1);
        histo_a->Draw();
        c_a->SaveAs(&("real_time_calibration/"+name_a+".png")[0]);

        TCanvas *c_b = new TCanvas(&(name_b)[0] ,&(name_b)[0]);
        histo_b->SetLineColor(2);
        histo_b->Draw();
        c_b->SaveAs(&("real_time_calibration/"+name_b+".png")[0]);

        float int_NA_a=histo_a->Integral(histo_a->FindFixBin(ranges[0]),histo_a->FindFixBin(ranges[1]));
        float int_NA_b=histo_b->Integral(histo_b->FindFixBin(ranges[0]),histo_b->FindFixBin(ranges[1]));
        cout<< int_NA_a/int_NA_b << endl;

        float int_cs_a=histo_a->Integral(histo_a->FindFixBin(ranges[2]),histo_a->FindFixBin(ranges[3]));
        float int_cs_b=histo_b->Integral(histo_b->FindFixBin(ranges[2]),histo_b->FindFixBin(ranges[3]));
        cout<< int_cs_a/int_cs_b << endl;

        float y_rescale = (int_NA_a/int_NA_b + int_cs_a/int_cs_b )/2;
        
        auto gaus_param_NAb=fit_gaus_NAb(histo_a, ranges, name_final, type);

        histos=make_histo(name_a, ranges[4], ranges[5], y_rescale, gaus_param_NAb[0]);

        auto gaus_param_NAb_clear=fit_gaus_NAb(histos[0], ranges, name_final, type);


        TH1F *histo_sub = (TH1F*)histo_b->Clone("histo_sub");

        histo_sub->Add(histos[0], -1);
        TCanvas *c_all = new TCanvas(&(name_final + type +"_all")[0] ,&(name_final + type +"_all")[0]);
        histo_sub->SetName(&(name_final + type +"_all")[0]);
        histo_sub->SetLineColor(4);
        histo_sub->Draw();
        histo_b->Draw("SAME");
        histo_a->Draw("SAME");
        histos[0]->Draw("SAME");
        c_all->SaveAs(&("real_time_calibration/"+name_final + type +"_all.png")[0]);
        
        histo_a->Write();
        histo_b->Write();
        histos[0]->Write();
        c_all->Write();

        vector<float> gaus_params;
        gaus_params=fit_gaus(histo_sub, ranges, name_final, type);
        auto lin_params = fit_lin(name_final, type, gaus_params, outfile, gaus_param_NAb_clear);

        vector<string> x= peak_energy(name_final, type, lin_params, gaus_params);
        energy.insert(energy.end(), x.begin(), x.end());
    }
    cout << endl;
    for (int i=0; i<energy.size(); i++){
        cout << energy[i];
        out_file << energy[i];
    }
    getchar();
    outfile->Close();

}