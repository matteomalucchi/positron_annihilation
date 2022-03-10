#include <list>
#include <string>


#include <TFile.h>
#include <TGraphErrors.h>

using namespace std;

vector<float> fit_gaus (TFile *f ,string name_pmt, string type, vector<float> range){

    list <string> samples{
        name_pmt + "_co_100",
        name_pmt + "_cs_100"
    };

    vector<float> fit_params;
    int i=0;

    for(list<string>::const_iterator sample = samples.begin(); sample != samples.end(); ++sample){
        TH1F *histo = nullptr;

        f->GetObject(&(*sample + type)[0], histo);
        TCanvas *c = new TCanvas(&(*sample + type)[0], &(*sample + type)[0]);
        TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[i+0],range[i+1]);
        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[i+2],range[i+3]);
        gaus1->SetParameters(2000, abs(range[i+0]+range[i+1])/2, abs(range[i+1]-range[i+0]));
        gaus2->SetParameters(2000, abs(range[i+2]+range[i+3])/2, abs(range[i+3]-range[i+2]));
        histo->Fit("gaus1","R", "SAME");
        gaus1->Draw("SAME");
        histo->Fit("gaus2","R", "SAME"); 
        gaus2->Draw("SAME");
        histo->Draw("SAME");
        gStyle->SetOptFit(111111);
        fit_params.push_back(gaus1->GetParameter(1));
        fit_params.push_back(gaus1->GetParameter(2));
        fit_params.push_back(gaus2->GetParameter(1));
        fit_params.push_back(gaus2->GetParameter(2));
        i +=4;
    }
    return fit_params;
}
/*
void fit_lin(vector<vector<float>> fit_params){
    double x[4]={0.18, 0.66, 1.17, 1.33};
    double y[4]={0.31, 1.01, 1.75, 1.98};
    double y_err[4] = {0.07, 0.03, 0.04, 0.04};

    TCanvas *c_fit_lin = new TCanvas("linear fit", "linear fit");
    TGraphErrors* gr = new TGraphErrors(4,x,y,nullptr,y_err);
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.5);
    gr->SetMarkerStyle(1);
    gr->Draw("APE");
    TF1 *   linear  = new TF1("linear","[0]+[1]*x", 0, 1.4);
    linear->SetParNames ("x_0 constant","Conversion factor");
    gr->Fit("linear");
    linear->Draw("SAME");
    gStyle->SetOptFit(111);
    c_fit_lin->SaveAs("images/linear_fit.png");
}
*/
void calibration(){
    map<string, vector<vector<float>>> pmts{
        {"pmt1",{{165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000},
            {165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000}}},
        {"pmt2",{{165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000},
            {165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000}}}
    };
    
    TFile *f = new TFile("histograms/histograms_calibration.root");


    for (const auto &pmt : pmts) {
        const auto name_pmt = pmt.first;
        const auto ranges = pmt.second;        
        vector<float> fit_charge=fit_gaus(f, name_pmt, "_charge", ranges[0]);
        cout << fit_charge[0]<< endl;
        //fit_lin(fit_params);
    }

}