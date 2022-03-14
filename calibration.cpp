#include <list>
#include <string>


#include <TFile.h>
#include <TGraphErrors.h>

using namespace std;

vector<float> fit_gaus (TFile *f ,string name_pmt, string type, vector<float> range){

    list <string> samples{
        name_pmt + "_co_100",
        name_pmt + "_cs_100", 
        name_pmt + "_NA_e6_100"
    };

    vector<float> gaus_params;
    int i=0;

    for(list<string>::const_iterator sample = samples.begin(); sample != samples.end(); ++sample){
        TH1F *histo = nullptr;

        cout << endl;
        cout << "__________________________ Gaussian fit: " << *sample << type << " __________________________"<<endl; 
        cout << endl;
        
        f->GetObject(&(*sample + type)[0], histo);
        TCanvas *c = new TCanvas(&(*sample + type)[0], &(*sample + type)[0]);
        TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[i+0],range[i+1]);
        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[i+2],range[i+3]);
        gaus1->SetParameters(1000, abs(range[i+0]+range[i+1])/2, abs(range[i+1]-range[i+0]));
        gaus2->SetParameters(1500, abs(range[i+2]+range[i+3])/2, abs(range[i+3]-range[i+2]));
        histo->Fit("gaus1","R", "SAME");
        gaus1->Draw("SAME");
        gStyle->SetOptFit(1111);
        histo->Fit("gaus2","R", "SAME"); 
        gaus2->Draw("SAME");
        histo->Draw("SAME");
        gStyle->SetOptFit(1111);

        gaus_params.push_back(gaus1->GetParameter(1));
        gaus_params.push_back(gaus1->GetParError(1)*10);
        gaus_params.push_back(gaus2->GetParameter(1));
        gaus_params.push_back(gaus2->GetParError(1)*10);
        i +=4;
        c->SaveAs(&("calibration/" + *sample + type + "_calibration.png")[0]);

    }
    return gaus_params;
}

auto fit_lin(string name_pmt, string type, vector<float> gaus_params){
    float x[4]={0.18, 0.66, 1.17, 1.33};
    float y[4]={gaus_params[4], gaus_params[6], gaus_params[0], gaus_params[2]};
    float y_err[4] = {gaus_params[5], gaus_params[7], gaus_params[1], gaus_params[3]};

    cout << endl;
    cout << "__________________________ Linear fit: " << name_pmt << type << " __________________________"<<endl; 
    cout << endl;

    TCanvas *c_fit_lin = new TCanvas(&(name_pmt + type + "_calibration")[0], &( name_pmt + type + "_calibration")[0]);
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
    c_fit_lin->SaveAs(&("calibration/" + name_pmt + type + "_calibration.png")[0]);

    vector<double> cal_factor{linear->GetParameter(1),linear->GetParError(1)};
    return cal_factor;
}

auto division(vector<double> cal_factor, float peak, float peak_err){
    vector<float> p;
    p.push_back(peak/cal_factor[0]);
    float e = sqrt(pow(peak_err/cal_factor[0], 2) + pow(peak*cal_factor[1]/(cal_factor[0]*cal_factor[0]), 2));
    p.push_back(e);
    return p;
}

auto peak_energy(string name_pmt, string type, vector<double> cal_factor, vector<float> gaus_params){
    vector<float> energy;
    vector<float> x(division(cal_factor, gaus_params[8], gaus_params[9]));    
    vector<float> y(division(cal_factor, gaus_params[10], gaus_params[11]));
    energy.insert(energy.end(), x.begin(), x.end());
    energy.insert(energy.end(), y.begin(), y.end());

    string peak1 = "energy peak " + name_pmt + " " + type + " = " + to_string(energy[0]) + "+-" + to_string(energy[1]) + "\n";
    string peak2 = "energy peak " + name_pmt + " " + type + " = " + to_string(energy[2]) + "+-" + to_string(energy[3])+ "\n";
    vector <string> energy_string;
    energy_string.push_back(peak1);
    energy_string.push_back(peak2);
    //cout<< "energy peak "<< name_pmt <<" "<< type << " = " << energy[0] <<"+-" << energy[1] <<endl;
    //cout<< "energy peak "<< name_pmt <<" "<< type << " = " << energy[2] <<"+-" << energy[3] <<endl;
    return energy_string;

}


void calibration(){
    // ordine del vector:
    // primo è per charge e secondo è per ampiezza
    // prima i 4 punti del co e poi i 4 del cs (e poi i 4 del NA)
    map<string, vector<vector<float>>> pmts{
        {"pmt1",{{165000, 185000, 185000, 210000, 26000, 37000, 93000, 113000, 65000, 90000, 175000, 210000},
            {2050, 2350, 2350, 2800, 300, 600, 1150, 1600, 850, 1300, 2200, 2700}}},
        {"pmt2",{{150000, 172000, 172000, 196000, 20000, 37000, 80000, 110000, 60000,90000, 165000, 200000},
            {1900, 2300, 2300, 2800, 250, 600, 1050, 1600, 800, 1200, 2100, 2600}}}
        /*{"pmt3",{{28000,33000,33000,38000,22000, 40000, 90000, 120000},
            {}}}*/
    };
    ofstream out_file("calibration/peak_energy.txt");
    TFile *f = new TFile("histograms/histograms_calibration.root");

    vector <string> energy;

    for (const auto &pmt : pmts) {
        const auto name_pmt = pmt.first;
        const auto ranges = pmt.second;  

        cout << endl; 
        cout << "========================= Analyzing: " << name_pmt << " ========================="<<endl;    
        cout << endl;

        vector<float> gaus_charge=fit_gaus(f, name_pmt, "_charge", ranges[0]);
        vector<float> gaus_amp=fit_gaus(f, name_pmt, "_amp", ranges[1]);
        auto cal_charge = fit_lin(name_pmt, "_charge", gaus_charge);
        auto cal_amp = fit_lin(name_pmt, "_amp", gaus_amp);
        vector<string> x= peak_energy(name_pmt, "_charge", cal_charge, gaus_charge);
        vector<string> y= peak_energy(name_pmt, "_amp", cal_amp, gaus_amp);
        energy.insert(energy.end(), x.begin(), x.end());
        energy.insert(energy.end(), y.begin(), y.end());

    }
    cout << endl;
    for (int i=0; i<energy.size(); i++){
        cout << energy[i];
        out_file << energy[i];
    }


}

