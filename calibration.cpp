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
        fit_params.push_back(gaus1->GetParError(1));
        fit_params.push_back(gaus2->GetParameter(1));
        fit_params.push_back(gaus2->GetParError(1));
        i +=4;
    }
    return fit_params;
}

void fit_lin(string name_pmt, string type, vector<float> fit_params){
    double x[4]={0.18, 0.66, 1.17, 1.33};
    double y[4]={fit_params[4], fit_params[6], fit_params[0], fit_params[2]};
    double y_err[4] = {fit_params[5], fit_params[7], fit_params[1], fit_params[3]};

    TCanvas *c_fit_lin = new TCanvas(&("images/" + name_pmt + type + "_calibration")[0], &("images/" + name_pmt + type + "_calibration")[0]);
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
    c_fit_lin->SaveAs(&("images/" + name_pmt + type + "_calibration.png")[0]);
}

void calibration(){
    // ordine del vector:
    // primo è per charge e secondo è per ampiezza
    // prima i 4 punti del co e poi i 4 del cs (e poi i 4 del NA)
    map<string, vector<vector<float>>> pmts{
        {"pmt1",{{165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000},
            {2300, 2350, 2350, 2800, 300, 600, 1150, 1600}}}
        /*{"pmt2",{{165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000},
            {165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000}}},
        {"pmt3",{{165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000},
            {165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000}}}*/
    };
    
    TFile *f = new TFile("histograms/histograms_calibration.root");


    for (const auto &pmt : pmts) {
        const auto name_pmt = pmt.first;
        const auto ranges = pmt.second;        
        vector<float> fit_charge=fit_gaus(f, name_pmt, "_charge", ranges[0]);
        for (int i=0; i<fit_charge.size(); i++){
        cout << fit_charge[i]<< endl;            
        }
        vector<float> fit_amp=fit_gaus(f, name_pmt, "_amp", ranges[1]);
        //cout << fit_amp[0]<< endl;
        fit_lin(name_pmt, "_charge", fit_charge);
        fit_lin(name_pmt, "_amp", fit_amp);
    }

}

// aggiungi funzione che converte picco NA da digit a ev e propaga errore