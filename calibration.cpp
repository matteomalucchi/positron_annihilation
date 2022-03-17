#include <list>
#include <string>


#include <TFile.h>
#include <TPad.h>
#include <TGraphErrors.h>

using namespace std;

vector<float> fit_gaus (TFile *f ,string name_pmt, string type, vector<float> range, TFile *outfile){

    list <string> samples{
        name_pmt + "_co_100",
        name_pmt + "_cs_100", 
        name_pmt + "_NA_e6_100_run2",
        //name_pmt + "_NA_e6_ext_run2"
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
        //gaus_params.push_back(gaus1->GetParameter(2));
        gaus_params.push_back(gaus1->GetParError(1));
        gaus_params.push_back(gaus2->GetParameter(1));
        //gaus_params.push_back(gaus2->GetParameter(2));
        gaus_params.push_back(gaus2->GetParError(1));
        i +=4;
        c->SaveAs(&("calibration/" + *sample + type + "_calibration_low.png")[0]);
        outfile->cd();
        c->Write();


    }
    return gaus_params;
}

auto fit_lin(string name_pmt, string type, vector<float> gaus_params, TFile * outfile){
    float x[4]={0.18, 0.66, 1.17, 1.33};
    float y[sizeof(x)/sizeof(x[0])]={gaus_params[4], gaus_params[6], gaus_params[0], gaus_params[2]};
    float y_err[sizeof(x)/sizeof(x[0])] = {gaus_params[5], gaus_params[7], gaus_params[1], gaus_params[3]};

    cout << endl;
    cout << "__________________________ Linear fit: " << name_pmt << type << " __________________________"<<endl; 
    cout << endl;
    gStyle->SetOptStat(0);

    TCanvas *c_fit_lin = new TCanvas(&(name_pmt + type + "_calibration")[0], &( name_pmt + type + "_calibration")[0]);
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

    TGraphErrors* gr = new TGraphErrors(4,x,y,nullptr,y_err);
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.5);
    gr->SetMarkerStyle(1);
    gr->Draw("APE");
    TF1 *   linear  = new TF1("linear","[0]+[1]*x", -1, 1.4);
    linear->SetParNames ("Off-set","Calibration factor");
    gr->Fit("linear");
    linear->Draw("SAME");
    gStyle->SetOptFit(111);

    pad2->cd();
    float diff[sizeof(x)/sizeof(x[0])];
    float diff_err[sizeof(x)/sizeof(x[0])];
    for (Int_t i=0;i<sizeof(x)/sizeof(x[0]);i++) {
        diff[i] =y[i]-linear->Eval(x[i]);
        diff_err[i]=diff[i]/y_err[i];
    }  

    TGraphErrors* gr2 = new TGraphErrors(4,x,diff_err,nullptr,y_err);
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.5);
    gr2->SetMarkerStyle(1);
    TF1 *   zero  = new TF1("zero","0*x", -0.1, 1.4);
    gr2->Draw("APE");
    zero->Draw("SAME");
    c_fit_lin->cd();
    c_fit_lin->SaveAs(&("calibration/" + name_pmt + type + "_calibration_low.png")[0]);
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

auto peak_energy(string name_pmt, string type, vector<vector<double>> cal_params, vector<float> gaus_params){
    vector<float> energy;
    vector<float> x(division(cal_params[0], cal_params[1], gaus_params[8], gaus_params[9]));    
    vector<float> y(division(cal_params[0], cal_params[1], gaus_params[10], gaus_params[11]));
    energy.insert(energy.end(), x.begin(), x.end());
    energy.insert(energy.end(), y.begin(), y.end());

    string peak1 = "energy peak " + name_pmt + " " + type + " = " + to_string(energy[0]) + "+-" + to_string(energy[1]) + "\n";
    string peak2 = "energy peak " + name_pmt + " " + type + " = " + to_string(energy[2]) + "+-" + to_string(energy[3])+ "\n";
    vector <string> energy_string;
    energy_string.push_back(peak1);
    energy_string.push_back(peak2);

    return energy_string;

}


void calibration(){
    // ordine del vector:
    // primo è per charge e secondo è per ampiezza
    // prima i 4 punti del co, poi i 4 del cs e poi i 4 del NA
    map<string, vector<vector<float>>> pmts{
        {"pmt1",{{170000, 180000, 192000, 204000,
                 26000, 35000, 97000, 106000,
                 65000, 90000, 175000, 210000,
                 //72000, 80000, 182000, 192000
                 },
                {2100, 2350, 2350, 2600,
                 370, 500, 1230, 1400, 
                 850, 1300, 2200, 2700,
                 //900,1100, 2250, 2450
                 }}},
        /*{"pmt2",{{150000, 172000, 172000, 196000, 
                 20000, 37000, 80000, 110000, 
                 60000,90000, 165000, 200000},
                {1900, 2300, 2300, 2800, 
                 250, 600, 1050, 1600,
                 800, 1200, 2100, 2600}}}*/
        /*{"pmt3",{{28000,33000,33000,38000,22000, 40000, 90000, 120000},
            {}}}*/
    };
    ofstream out_file("calibration/peak_energy_low.txt");
    TFile *outfile= new TFile("calibration/calibration_plots.root", "RECREATE");
    TFile *f = new TFile("histograms/histograms.root");

    vector <string> energy;

    for (const auto &pmt : pmts) {
        const auto name_pmt = pmt.first;
        const auto ranges = pmt.second;  

        cout << endl; 
        cout << "========================= Analyzing: " << name_pmt << " ========================="<<endl;    
        cout << endl;

        vector<float> gaus_charge=fit_gaus(f, name_pmt, "_charge", ranges[0], outfile);
        vector<float> gaus_amp=fit_gaus(f, name_pmt, "_amp", ranges[1], outfile);
        auto cal_charge = fit_lin(name_pmt, "_charge", gaus_charge, outfile);
        auto cal_amp = fit_lin(name_pmt, "_amp", gaus_amp, outfile);
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
    outfile->Close();



}

