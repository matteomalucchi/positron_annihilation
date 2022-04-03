#include <list>
#include <string>


#include <TFile.h>
#include <TPad.h>
#include <TGraphErrors.h>

using namespace std;

ofstream energy_file("triple/energy.txt");
ofstream params_file("triple/params.txt");



vector<vector<float>> fit_gaus (TFile *f ,string name_pmt, string type, vector<float> range, TFile *outfile){

    list <string> samples{
       /* name_pmt + "_NA_e6_ext_triple_90deg_run1",
        name_pmt + "_NA_l1_ext_triple_close_run2", 
        name_pmt + "_NA_l1_ext_triple_close_run3",
        name_pmt + "_NA_c6_ext_triple_merc_aero_run4",
        name_pmt + "_NA_c6_ext_coinc12_merc_metal_run5", */
        name_pmt + "_NA_c6_ext_coinc12_merc_metal_run6",
    };

    vector<vector<float>> gaus_params(samples.size());
    int i=0;

    for(list<string>::const_iterator sample = samples.begin(); sample != samples.end(); ++sample){
        TH1F *histo = nullptr;

        cout << endl;
        cout << "__________________________ Gaussian fit: " << *sample << type << " __________________________"<<endl; 
        cout << endl;
        
        f->GetObject(&(*sample + type)[0], histo);
        TCanvas *c = new TCanvas(&(*sample + type)[0], &(*sample + type)[0]);
        TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[0],range[1]);
        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[2],range[3]);
        gaus1->SetParameters(1000, abs(range[0]+range[1])/2, abs(range[1]-range[0]));
        gaus2->SetParameters(1000, abs(range[2]+range[3])/2, abs(range[3]-range[2]));
        histo->Fit("gaus1","R", "SAME");
        gaus1->Draw("SAME");
        gStyle->SetOptFit(1111);
        histo->Fit("gaus2","R", "SAME"); 
        gaus2->Draw("SAME");
        histo->Draw("SAME");
        gStyle->SetOptFit(1111);

        gaus_params[i].push_back(gaus1->GetParameter(1));
        //gaus_params.push_back(gaus1->GetParameter(2));
        gaus_params[i].push_back(gaus1->GetParError(1));
        gaus_params[i].push_back(gaus2->GetParameter(1));
        //gaus_params.push_back(gaus2->GetParameter(2));
        gaus_params[i].push_back(gaus2->GetParError(1));
        i++;
        c->SaveAs(&("calibration/" + *sample + type + "_calibration_low.png")[0]);
        outfile->cd();
        c->Write();
    }
    return gaus_params;
}

auto fit_lin(string name_pmt, string type, vector<vector<float>> gaus_params, TFile * outfile){
    float x[2]={0.511,/*1.061,*/ 1.274};
    vector<vector<double>> lin_params(gaus_params.size());

    for (int j=0;j<gaus_params.size();j++){
        float y[sizeof(x)/sizeof(x[0])]={gaus_params[j][0], gaus_params[j][2]};
        float y_err[sizeof(x)/sizeof(x[0])] = {gaus_params[j][1], gaus_params[j][3]};

        cout << endl;
        cout << "__________________________ Linear fit: " << name_pmt << type <<"_run"<< to_string(j+1) <<" __________________________"<<endl; 
        cout << endl;
        gStyle->SetOptStat(0);

        TCanvas *c_fit_lin = new TCanvas(&(name_pmt + type + "_run"+ to_string(j+1) +"_calibration")[0], &( name_pmt + type + "_run"+to_string(j+1)+"_calibration")[0]);
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
        //linear->SetParLimits(0,3000, 5000);
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
        float one_array[sizeof(x)/sizeof(x[0])]={1,1};

        TGraphErrors* gr2 = new TGraphErrors(sizeof(x)/sizeof(x[0]),x,diff_norm,nullptr,one_array);
        gStyle->SetStatY(0.9);
        gStyle->SetStatX(0.5);
        gr2->SetMarkerStyle(1);
        TF1 *   zero  = new TF1("zero","0*x", 0, 1.4);
        gr2->Draw("APE");
        zero->Draw("SAME");
        c_fit_lin->cd();
        c_fit_lin->SaveAs(&("calibration/" + name_pmt + type + "_run"+to_string(j+1)+ "_calibration.png")[0]);
        outfile->cd();
        c_fit_lin->Write();

        lin_params[j].push_back(linear->GetParameter(0));
        lin_params[j].push_back(linear->GetParError(0));
        lin_params[j].push_back(linear->GetParameter(1));
        lin_params[j].push_back(linear->GetParError(1));
    }
    return lin_params;
}
/*
auto division(vector<double> off_set, vector<double> cal_factor, float peak, float peak_err){
    vector<float> p;
    p.push_back((peak-off_set[0])/cal_factor[0]);
    float e = sqrt(pow(peak_err/cal_factor[0], 2) +pow(off_set[1]/cal_factor[0], 2) +pow((peak-off_set[0])*cal_factor[1]/(cal_factor[0]*cal_factor[0]), 2));
    p.push_back(e);
    return p;
}

auto peak_energy(string name_pmt, string type, vector<vector<double>> cal_params, vector<float> gaus_params){
    vector<float> energy;
    for (int j=0;j<cal_params.size();j++){

    vector<float> x(division(cal_params[j][0], cal_params[j][1], gaus_params[j][2], gaus_params[j][3]));    
    vector<float> y(division(cal_params[0], cal_params[1], gaus_params[10], gaus_params[11]));
    energy.insert(energy.end(), x.begin(), x.end());
    energy.insert(energy.end(), y.begin(), y.end());

    string peak1 = "energy peak " + name_pmt + " " + type + " = " + to_string(energy[0]) + "+-" + to_string(energy[1]) + "\n";
    string peak2 = "energy peak " + name_pmt + " " + type + " = " + to_string(energy[2]) + "+-" + to_string(energy[3])+ "\n";
    vector <string> energy_string;
    energy_string.push_back(peak1);
    energy_string.push_back(peak2);

    return energy_string;

}*/


void calibration_triple(){
    // ordine del vector:
    // primo è per charge e secondo è per ampiezza
    // prima i 2 del NA |poi i 2 punti del compton edge Ne | e i 2 del picco del Ne
    map<string, vector<vector<float>>> pmts{
        {"pmt1",{{70000, 79000, /* 140000, 152000,*/ 175000, 192000},
                {900, 1100, /*1830, 1920,*/ 2200, 2540}}},
        {"pmt2",{{63000, 73000,/* 130000, 137000,*/ 157000, 176000},
                {880, 1070, /*1760, 1900,*/ 2080, 2440}}},
        {"pmt3",{{11000, 14000, /*23500,26500,*/ 28000,34000},
                {160, 210,/* 350, 380,*/ 400, 470}}}
        };
    TFile *outfile= new TFile("triple/calibration_plots.root", "RECREATE");
    TFile *f = new TFile("histograms/histograms_triple.root");

    vector <string> energy;

    for (const auto &pmt : pmts) {
        const auto name_pmt = pmt.first;
        const auto ranges = pmt.second;  

        cout << endl; 
        cout << "========================= Analyzing: " << name_pmt << " ========================="<<endl;    
        cout << endl;

        auto gaus_charge=fit_gaus(f, name_pmt, "_charge", ranges[0], outfile);
        auto gaus_amp=fit_gaus(f, name_pmt, "_amp", ranges[1], outfile);
        auto cal_charge = fit_lin(name_pmt, "_charge", gaus_charge, outfile);
        auto cal_amp = fit_lin(name_pmt, "_amp", gaus_amp, outfile);
        for (int r=0;r<gaus_charge.size();r++){
            params_file<< "off_set "<<name_pmt<< "_run"<< to_string(r+1)<< "_charge = "<< cal_charge[r][0] << " +- " <<  cal_charge[r][1]<<"\n"; 
            params_file<< "cal_factor "<<name_pmt<< "_run"<< to_string(r+1)<< "_charge = "<< cal_charge[r][2] << " +- " <<  cal_charge[r][3]<<"\n"; 
            params_file<< "off_set "<<name_pmt<< "_run"<<to_string(r+1)<< "_amp = "<< cal_amp[r][0] << " +- " <<  cal_amp[r][1]<<"\n"; 
            params_file<< "cal_factor "<<name_pmt<< "_run"<<to_string(r+1)<< "_amp = "<< cal_amp[r][2] << " +- " <<  cal_amp[r][3]<<"\n"; 
            params_file <<"\n";
        }
        /*
        vector<string> x= peak_energy(name_pmt, "_charge", cal_charge, gaus_charge);
        vector<string> y= peak_energy(name_pmt, "_amp", cal_amp, gaus_amp);
        energy.insert(energy.end(), x.begin(), x.end());
        energy.insert(energy.end(), y.begin(), y.end());
*/

    }
    /*
    for (int i=0; i<energy.size(); i++){
        cout << energy[i];
        energy_file << energy[i];
    }
    */
    outfile->Close();
}

