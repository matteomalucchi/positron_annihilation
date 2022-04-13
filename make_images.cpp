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

// Na spetrum in run1 e run2 discrepancies


void make_images2(){
    TFile *f = new TFile("histograms/histograms_csvsco_3.root");

    TH1F *histo1 = nullptr;
    TH1F *histo2 = nullptr;
    string name = "pmt3_Na(csvsco)_charge";
    string name_nice = "pmt3 Na charge";
    TCanvas *c = new TCanvas(&name[0],&name[0]);
    gStyle->SetOptStat("neou");


    f->GetObject("pmt3_NA_cs_e6_30_run1_charge", histo1);
    f->GetObject("pmt3_NA_cs_co_e6_30_run1_charge", histo2);

    histo1->SetNameTitle("pmt3 Na (no Co) charge", &name_nice[0]);
    histo2->SetNameTitle("pmt3 Na (Co) charge", &name_nice[0]);
    histo1->GetXaxis()->SetTitle("Charge [a.u.]");
    histo2->GetXaxis()->SetTitle("Charge [a.u.]");
    histo1->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo2->GetBinWidth(1))))[0]);
    histo2->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo2->GetBinWidth(1))))[0]);

    histo1->SetLineColor(kBlue);
    histo2->SetLineColor(kGreen);

    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",12500,15000);
    gaus1->SetParName(0, "const 1");
    gaus1->SetParName(1, "#mu 1");
    gaus1->SetParName(2, "#sigma 1");
    TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",12500,15000);
    gaus2->SetParName(0, "const 2");
    gaus2->SetParName(1, "#mu 2");
    gaus2->SetParName(2, "#sigma 2");

    gaus1->SetParameters(100, 14000, 3000);
    gaus2->SetParameters(100,  14000, 3000);

    histo2->Draw();
    TPaveStats *st = (TPaveStats*)histo2->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.4); //new x end position
    st->SetY1NDC(0.3); //new x start position
    st->SetY2NDC(0.6); //new x end position

    histo2->Fit("gaus2","R", "SAMES");
    gaus2->SetLineColor(kBlack);
    gaus2->Draw("SAME");
    gStyle->SetOptFit(111);
    //histo2->Fit("gaus2","R", "SAME"); 
    //gaus2->Draw("SAME");

    histo1->Fit("gaus1","R", "SAMES");
    gaus1->SetLineColor(kRed);
    gaus1->Draw("SAME");
    gStyle->SetOptFit(111);
    //histo1->Fit("gaus2","R", "SAME"); 
    //gaus2->Draw("SAME");
    histo1->Draw("SAME");
    TPaveStats *st1 = (TPaveStats*)histo1->FindObject("stats");
    st1->SetX1NDC(0.1); //new x start position
    st1->SetX2NDC(0.4); //new x end position
    st1->SetY1NDC(0.6); //new x start position
    st1->SetY2NDC(0.9); //new x end position

   auto legend = new TLegend(0.7,0.9,0.7,0.9);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(histo1,"pmt3 Na (no Co) charge","l");
   legend->AddEntry(histo2,"pmt3 Na (Co) charge","l");
   legend->AddEntry("gaus1","gaus 1","l");
   legend->AddEntry("gaus2","gaus 2","l");
   legend->Draw();

    c->SaveAs(&("final_images/"+ name+".pdf")[0]);

}


// pmt1 spectra of co cs na
void make_image2s(){
    TFile *f = new TFile("histograms/histograms_cocsna.root");

    TH1F *histo1 = nullptr;
    TH1F *histo2 = nullptr;
    TH1F *histo3 = nullptr;
    string name = "pmt1_amp";
    string name_nice = "pmt1 amplitude";
    TCanvas *c = new TCanvas(&name[0],&name[0]);


    f->GetObject("pmt1_co_100_amp", histo1);
    f->GetObject("pmt1_cs_100_amp", histo2);
    f->GetObject("pmt1_NA_e6_100_run2_amp", histo3);


    histo1->SetNameTitle("pmt1 Co amplitude", &name_nice[0]);
    histo2->SetNameTitle("pmt1 Cs amplitude", &name_nice[0]);
    histo3->SetNameTitle("pmt1 Na amplitude", &name_nice[0]);


    //gStyle->SetOptStat("neou");

    histo1->GetXaxis()->SetTitle("Amplitude [a.u.]");
    histo2->GetXaxis()->SetTitle("Amplitude [a.u.]");
    histo3->GetXaxis()->SetTitle("Amplitude [a.u.]");
    histo1->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo1->GetBinWidth(1))))[0]);
    histo2->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo2->GetBinWidth(1))))[0]);
    histo3->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo3->GetBinWidth(1))))[0]);

    histo1->SetLineColor(kBlue);
    histo2->SetLineColor(kGreen);
    histo3->SetLineColor(kRed);

    histo2->Draw();

    c->Update();
    TPaveStats *st2 = (TPaveStats*)histo2->FindObject("stats");
    st2->SetX1NDC(0.65); //new x start position
    st2->SetX2NDC(0.9); //new x end position
    st2->SetY1NDC(0.7); //new x start position
    st2->SetY2NDC(0.9); //new x end position

    histo3->Draw("sames");
    c->Update();
    TPaveStats *st = (TPaveStats*)histo3->FindObject("stats");
    st->SetX1NDC(0.65); //new x start position
    st->SetX2NDC(0.9); //new x end position
    st->SetY1NDC(0.5); //new x start position
    st->SetY2NDC(0.7); //new x end position

    histo1->Draw("SAMES");
    c->Update();
    TPaveStats *st1 = (TPaveStats*)histo1->FindObject("stats");
    st1->SetX1NDC(0.65); //new x start position
    st1->SetX2NDC(0.9); //new x end position
    st1->SetY1NDC(0.3); //new x start position
    st1->SetY2NDC(0.5); //new x end position

    auto legend = new TLegend(0.1,0.75,0.35,0.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(histo2,"pmt1 Cs amplitude","l");
    legend->AddEntry(histo3,"pmt1 Na amplitude","l");
    legend->AddEntry(histo1,"pmt1 Co amplitude","l");
    //legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
    legend->Draw();

    c->SaveAs(&("final_images/"+ name+".pdf")[0]);

}



void make_images(){
    TFile *f= new TFile(&("triple/ntuple_triple_quad.root")[0]);

    TNtuple *ntuple;
    int bins=150;
    float min=-20, max=30;


    f->GetObject("run7", ntuple);
    TCanvas *c = new TCanvas("c1", "c1");
    TH1F* histo1= new TH1F("t1-t2", "t1-t2", bins, min, max);
    gStyle->SetOptStat("neou");

    for (auto ievt: ROOT::TSeqI(ntuple->GetEntries())){// Loop over events
        ntuple->GetEntry(ievt);
        auto catm = ntuple->GetArgs(); // Get a row
        histo1->Fill(catm[2]-catm[5]);
    }
    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",-15,20);
    gaus1->SetParName(0, "const 1");
    gaus1->SetParName(1, "#mu 1");
    gaus1->SetParName(2, "#sigma 1");
    gaus1->SetParameters(100, 5, 5);

    histo1->Fit("gaus1");
    histo1->Draw();
    histo1->GetXaxis()->SetTitle("Time [ns]");
    histo1->GetYaxis()->SetTitle(&("Entries / "+to_string(histo1->GetBinWidth(1))+"[ns]")[0]);
    histo1->Draw();
    gStyle->SetOptFit(1111); 


    c->Update();
    TPaveStats *st1 = (TPaveStats*)histo1->FindObject("stats");
    st1->SetX1NDC(0.65); //new x start position
    st1->SetX2NDC(0.9); //new x end position
    st1->SetY1NDC(0.65); //new x start position
    st1->SetY2NDC(0.9); //new x end position



    c->SaveAs(&("final_images/t12_run7.pdf")[0]);

}

string type_of_file = "trasl";

vector <TH1F*> make_histo(string name, float charge_min, float charge_max,float amp_min, float amp_max, float y_rescale_charge,float y_rescale_amp, float x_rescale_charge, float x_rescale_amp){
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
        int k=0, m, j=0, f=0;
        while(getline(myfile, tp) /*&& f<50000*/){ //read data from file object and put it into string.
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
            f++;
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
        for (int j=0; j<20; j++){
            h += v[i][j] / 20;
            if (v[i][j]<14600 && u==0){
                idx_strange.push_back(i);
                u++;
            }
        }
        for (int j=300; j<950; j++){
            charge[i] += h-v[i][j];
        }
        min =*min_element(v[i].begin()+300, v[i].end()-80);
        amp[i]=h-min;
    }
    int range_charge=250000;
    int range_amp=3500; 
    if (name.find("pmt3")<name.length()){
        range_charge=50000;
        range_amp=600;
    }
    TH1F *histo_charge=  new TH1F(&(name + "_charge")[0],&(name + "_charge")[0], 950, 0, range_charge);
    TH1F *histo_amp=  new TH1F(&(name + "_amp")[0],&(name + "_amp")[0], 450, 0, range_amp);
    for (long unsigned int i=0; i< charge.size(); i++){
        if ((find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()) && charge[i]> charge_min && charge[i]<charge_max){
            if (type_of_file.find("trasl") < type_of_file.length()) histo_charge->Fill(charge[i]+x_rescale_charge);
            else histo_charge->Fill(charge[i]);
        }
    }
    for (long unsigned int i=0; i< charge.size(); i++){
        if ((find(idx_strange.begin(), idx_strange.end(), i) == idx_strange.end()) && amp[i]> amp_min && amp[i]<amp_max){
            if (type_of_file.find("trasl") < type_of_file.length()) histo_amp->Fill(amp[i]+x_rescale_amp);
            else histo_amp->Fill(amp[i]);        
        }
    }
    histo_charge->Scale(y_rescale_charge);
    histo_amp->Scale(y_rescale_amp);

    histo_charge->SetName(&(name + "_charge")[0]);
    histo_amp->SetName(&(name + "_amp")[0]);
    histos.push_back(histo_charge);
    histos.push_back(histo_amp);

    histo_charge->SetLineColor(3);
    histo_amp->SetLineColor(3);
    return histos;
}

vector<float> fit_gaus_single (TH1F* histo, float ranges1, float ranges2 , string name, string type){
    
    cout << endl;
    cout << "__________________________ Gaussian fit: " << name << type << " __________________________"<<endl; 
    cout << endl;

    int n_bin = histo->GetNbinsX();
    float x_max = histo->GetXaxis()->GetBinCenter(n_bin);


    TCanvas *c_Ne = new TCanvas(&(name + type +"_fit")[0] ,&(name + type +"_fit")[0]);
    TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
    pad1->SetBottomMargin(0.0001);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.00001);
    pad2->SetBottomMargin(0.1);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",ranges1, ranges2);
    gaus1->SetParameters(200, abs(ranges1+ranges2)/2, abs(ranges1-ranges2));
    histo->Fit("gaus1","R", "SAME");
    //gStyle->SetStatY(0.9);
    //gStyle->SetStatX(0.5);
    histo->GetYaxis()->SetTitle("Entries");
    //histo->GetXaxis()->SetTitle("Bin");
    histo->Draw("HIST");
    gaus1->Draw("SAME");
    //gStyle->SetOptFit(1111); 

    pad2->cd();
    TH1F *h2 = new TH1F("h2","residuals",n_bin,0,x_max);
    TF1 *zero = new TF1("zero","0*x",0, x_max);

    /*h2->GetXaxis()->SetLabelFont(63);
    h2->GetXaxis()->SetLabelSize(16);
    h2->GetXaxis()->SetTitle("Bin");
    h2->GetYaxis()->SetLabelFont(63);
    h2->GetYaxis()->SetLabelSize(16);*/
    for (int i=histo->FindFixBin(ranges1);i<=histo->FindFixBin(ranges2);i++) {
        float diff = (histo->GetBinContent(i)-gaus1->Eval(histo->GetBinCenter(i)))/histo->GetBinError(i);
        h2->SetBinContent(i,diff);
        h2->SetBinError(i,1);
    }
    h2->Draw("E"); //"P"
    zero->Draw("SAME");
    c_Ne->cd();

    
    vector<float> gaus_params;
    gaus_params.push_back(gaus1->GetParameter(1));
    gaus_params.push_back(gaus1->GetParError(1));

    return gaus_params;
}


void make_images1(){
    gROOT->SetBatch(kFALSE);
    ROOT::EnableImplicitMT();

    TFile *f = new TFile("histograms/histograms_realtimecalib.root");


    // primi due sono picco NA(solo cs)| poi picco cs(solo cs)| poi range da sottrare del NA (solo cs)|e poi il range da fittare del picco del Ne (solo cs)|
    // poi il range del primo picco co (cs+co)|e poi del secondo picco co(cs+co)
    map <vector<string>, vector<vector<float>>> pair_names ={ 
        {{"pmt1_NA_cs_e6_100_run1", "pmt1_NA_cs_co_e6_100_run1"}, 
            {{78000-2000, 85000-2000, 100000-2000, 108000-2000, 150000-2000, 210000-2000, 190000-2000, 204000-2000, 176000-2000, 187000-2000, 200000-2000, 212000-2000},
            {980, 1120, 1250, 1450, 1900, 2850, 2300, 2600, 2100, 2400, 2400, 2700}}}, 
/*
            {{"pmt1_NA_cs_e6_500_run2", "pmt1_NA_cs_co_e6_500_run2"}, 
            {{78000, 85000, 100000, 108000, 150000, 210000, 190000, 204000, 176000, 187000, 200000, 212000}, 
            {1025, 1150, 1275, 1475, 1900, 2850, 2350, 2650, 2225, 2425, 2425, 2800}}}, */
        };
    
    vector <float>   y_rescale, discrepancy;
    vector<vector<float>> energy_piccoNe(2), energy(2), gaus_param_Ne, gaus_param_Na_cs, gaus_param_Na_co, gaus_param_Cs_cs, gaus_param_Cs_co;
    vector<vector<vector<vector<double>>>> cal_params(2);

    string name_nice="pmt1 charge";

    list<string> types= {"_charge"};
    list<string> pmts= {"_pmt1", "_pmt2", "_pmt3"};


    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first[0];
        const auto name_b = pair_name.first[1];
        const auto ranges = pair_name.second;
        vector <TH1F*> histos;
        
        string pmt = name_a.substr(0, 4);
        string final = name_a.substr(10, name_a.size()-1);
        string name_final = pmt+final;

        cout << "Processing: " << name_final << endl; 
        int i=0;
        for(list<string>::const_iterator type = types.begin(); type != types.end(); ++type){
            TH1F *histo_a = nullptr;
            f->GetObject(&(name_a + *type)[0], histo_a);
            TH1F *histo_b = nullptr;
            f->GetObject(&(name_b + *type)[0], histo_b);
            histo_a->SetLineColor(1);
            histo_b->SetLineColor(7);

            float int_NA_a=histo_a->Integral(histo_a->FindFixBin(ranges[i][0]),histo_a->FindFixBin(ranges[i][1]));
            float int_NA_b=histo_b->Integral(histo_b->FindFixBin(ranges[i][0]),histo_b->FindFixBin(ranges[i][1]));
            //cout<< int_NA_a/int_NA_b << endl;

            float int_cs_a=histo_a->Integral(histo_a->FindFixBin(ranges[i][2]),histo_a->FindFixBin(ranges[i][3]));
            float int_cs_b=histo_b->Integral(histo_b->FindFixBin(ranges[i][2]),histo_b->FindFixBin(ranges[i][3]));
            //cout<< int_cs_a/int_cs_b << endl;

            y_rescale.push_back((int_NA_b/int_NA_a + int_cs_b/int_cs_a)/2);
            
            gaus_param_Ne.push_back(fit_gaus_single(histo_a, ranges[i][6],ranges[i][7] , name_final+"_Ne", *type));

            gaus_param_Na_cs.push_back(fit_gaus_single(histo_a, ranges[i][0],ranges[i][1] , name_final+"_Na_cs", *type));
            gaus_param_Na_co.push_back(fit_gaus_single(histo_b, ranges[i][0],ranges[i][1] , name_final+"_Na_co", *type));
            gaus_param_Cs_cs.push_back(fit_gaus_single(histo_a, ranges[i][2],ranges[i][3] , name_final+"_Cs_cs", *type));
            gaus_param_Cs_co.push_back(fit_gaus_single(histo_b, ranges[i][2],ranges[i][3] , name_final+"_Cs_co", *type));
            discrepancy.push_back(((gaus_param_Na_co[i][0]-gaus_param_Na_cs[i][0])+(gaus_param_Cs_co[i][0]-gaus_param_Cs_cs[i][0]))/2);

            i++;

        }

        histos=make_histo(name_a, ranges[0][4], ranges[0][5],ranges[1][4], ranges[1][5], 
                        y_rescale[0], y_rescale[1], discrepancy[0], discrepancy[1]);
        
        i=0;

        for(list<string>::const_iterator type = types.begin(); type != types.end(); ++type){
            TH1F *histo_a = nullptr;
            f->GetObject(&(name_a + *type)[0], histo_a);
            TH1F *histo_b = nullptr;
            f->GetObject(&(name_b + *type)[0], histo_b);

            //auto gaus_param_Ne_clear=fit_gaus_single(histos[i], ranges[i][6],ranges[i][7], name_final+"_clear", *type);
            TH1F *histo_sub = (TH1F*)histo_b->Clone(&("histo_sub"+*type)[0]);
            
            histo_sub->Add(histos[i], -1);
            TCanvas *c = new TCanvas(&(name_final + *type +"_all")[0] ,&(name_final + *type +"_all")[0]);

            TH1F *histos_a=  new TH1F("a", "a", 437, 105000, 220000);
            TH1F *histos_b=  new TH1F("a", "a", 437, 105000, 220000);
            TH1F *histoss=  new TH1F("a", "a", 437, 105000, 220000);
            TH1F *histos_sub=  new TH1F("a", "a", 437, 105000, 220000);

            int ua=0, oa=0, ub=0, ob=0, ui=0, oi=0, us=0, os=0;
            for (int y=0;y<950;y++){
                if (y<399){//underflow
                ua+=histo_a->GetBinContent(y);
                ub+=histo_b->GetBinContent(y);
                ui+=histos[i]->GetBinContent(y);
                us+=histo_sub->GetBinContent(y);
                }
                else if (y>836){//overflow
                oa+=histo_a->GetBinContent(y);
                ob+=histo_b->GetBinContent(y);
                oi+=histos[i]->GetBinContent(y);
                os+=histo_sub->GetBinContent(y);
                }
                else{
                    for(int p=0;p<histo_a->GetBinContent(y);p++){
                    histos_a->Fill(histo_a->GetBinCenter(y));
                    }
                    for(int p=0;p<histo_b->GetBinContent(y);p++){                    
                        histos_b->Fill(histo_b->GetBinCenter(y));
                    }
                    for(int p=0;p<histos[i]->GetBinContent(y);p++){
                        histoss->Fill(histos[i]->GetBinCenter(y));
                    }
                    for(int p=0;p<histo_sub->GetBinContent(y);p++){
                        histos_sub->Fill(histo_sub->GetBinCenter(y));
                    }                   
                }            
            }

            histos_a->SetBinContent(0, ua);
            histos_b->SetBinContent(0, ub);
            histoss->SetBinContent(0, ui);
            histos_sub->SetBinContent(0, us);
            
            histos_a->SetBinContent(229, oa);
            histos_b->SetBinContent(229, ob);
            histoss->SetBinContent(229, oi);
            histos_sub->SetBinContent(229, os);
                
                
            gStyle->SetOptStat("neou");


            histos_a->SetNameTitle("Na + Cs",  &name_nice[0]);
            histos_b->SetNameTitle("Na + Cs + Co", &name_nice[0]);
            histoss->SetNameTitle("Ne", &name_nice[0]);
            histos_sub->SetNameTitle("Na + Cs + Co + Ne", &name_nice[0]);


            histos_a->GetXaxis()->SetTitle("Charge [a.u.]");
            histos_b->GetXaxis()->SetTitle("Charge [a.u.]");
            histoss->GetXaxis()->SetTitle("Charge [a.u.]");
            histos_sub->GetXaxis()->SetTitle("Charge [a.u.]");
            histos_a->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histos_a->GetBinWidth(1))))[0]);
            histos_b->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histos_b->GetBinWidth(1))))[0]);
            histoss->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histoss->GetBinWidth(1))))[0]);
            histos_sub->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histos_sub->GetBinWidth(1))))[0]);

            histos_a->SetLineColor(2);
            histos_b->SetLineColor(7);
            histoss->SetLineColor(3);
            histos_sub->SetLineColor(4);


            histos_sub->Draw("HIST");
            /*histos_b->Draw("SAMES HIST");
            histos_a->Draw("SAMES HIST");
            histoss->Draw("SAMES HIST");*/

            c->Update();
            TPaveStats *st2 = (TPaveStats*)histos_sub->FindObject("stats");
            st2->SetX1NDC(0.65); //new x start position
            st2->SetX2NDC(0.9); //new x end position
            st2->SetY1NDC(0.75); //new x start position
            st2->SetY2NDC(0.9); //new x end position

            histoss->Draw("sames");
            c->Update();
            TPaveStats *st = (TPaveStats*)histoss->FindObject("stats");
            st->SetX1NDC(0.65); //new x start position
            st->SetX2NDC(0.9); //new x end position
            st->SetY1NDC(0.6); //new x start position
            st->SetY2NDC(0.75); //new x end position

            histos_a->Draw("SAMES");
            c->Update();
            TPaveStats *st1 = (TPaveStats*)histos_a->FindObject("stats");
            st1->SetX1NDC(0.4); //new x start position
            st1->SetX2NDC(0.65); //new x end position
            st1->SetY1NDC(0.75); //new x start position
            st1->SetY2NDC(0.9); //new x end position

            histos_b->Draw("SAMES");
            c->Update();
            TPaveStats *st3 = (TPaveStats*)histos_b->FindObject("stats");
            st3->SetX1NDC(0.4); //new x start position
            st3->SetX2NDC(0.65); //new x end position
            st3->SetY1NDC(0.6); //new x start position
            st3->SetY2NDC(0.75); //new x end position

            auto legend = new TLegend(0.2,0.75,0.4,0.9);
            //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
            legend->AddEntry(histos_b,"Na + Cs +Co","l");
            legend->AddEntry(histoss,"Ne","l");
            legend->AddEntry(histos_a,"Na + Cs","l");
            legend->AddEntry(histos_sub,"Na + Cs +Co + Ne","l");
            legend->Draw();

            c->SaveAs(&("final_images/pmt1_charge_all1.pdf")[0]);
        }
    }
}