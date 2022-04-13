#include <list>
#include <string>


#include <TFile.h>
#include <TPad.h>
#include <TGraphErrors.h>

using namespace std;

ofstream params_file("triple_calib/quad_params.txt");
ofstream chi_file("triple_calib/chi_square.txt");
ofstream gaus_file("triple_calib/gaus_params.txt");

vector<float> fit_gaus (TFile *f ,string name_pmt, string type, vector<float> range, TFile *outfile){

    vector<float> gaus_params;

    TH1F *histo = nullptr;

    cout << endl;
    cout << "__________________________ Gaussian fit: " << name_pmt << type << " __________________________"<<endl; 
    cout << endl;
    
    f->GetObject(&(name_pmt + type)[0], histo);
    TCanvas *c = new TCanvas(&(name_pmt + type)[0], &(name_pmt + type)[0]);
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

    gaus_params.push_back(gaus1->GetParameter(1));
    gaus_params.push_back(gaus1->GetParError(1));
    gaus_params.push_back(gaus2->GetParameter(1));
    gaus_params.push_back(gaus2->GetParError(1));
    c->SaveAs(&("triple_calib/" + name_pmt + type + "_calibration_low.png")[0]);
    outfile->cd();
    c->Write();
    chi_file<< name_pmt + type << "   chi2/ndof  NA   =  "<< gaus1->GetChisquare() <<"/"<< gaus1->GetNDF() <<"\n";
    gaus_file<< name_pmt + type << "   peak NA  =   "<< gaus1->GetParameter(1) <<"+-"<< gaus1->GetParError(1) <<"\n";
    chi_file<< name_pmt + type << "   chi2/ndof  Ne   =  "<< gaus2->GetChisquare() <<"/"<< gaus2->GetNDF() <<"\n";
    gaus_file<< name_pmt + type << "   peak Ne  =   "<< gaus2->GetParameter(1) <<"+-"<< gaus2->GetParError(1) <<"\n";
    return gaus_params;
}

auto fit_quad(string name_final, string type, vector<float> gaus_params, TFile * outfile){
    // na + ne
    float x[2]={0.5106, 1.274};
    float y[sizeof(x)/sizeof(x[0])]={gaus_params[0], gaus_params[2]};
    float y_err[sizeof(x)/sizeof(x[0])] = {gaus_params[1], gaus_params[3]};

    cout << endl;
    cout << "__________________________ Quadratic fit: " << name_final << type << " __________________________"<<endl; 
    cout << endl;
    gStyle->SetOptStat(0);

    TCanvas *c_fit_quad = new TCanvas(&(name_final + type + "_fit_quad")[0], &( name_final + type + "_fit_quad")[0]);
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
    gr->GetYaxis()->SetTitle("Peak [a.u.]");
    gr->GetXaxis()->SetTitle("Energy [MeV]");
    gr->GetXaxis()->SetTitleSize(0.15);
    gr->GetXaxis()->SetLimits(0, 1.4);

    gr->Draw("APE");
    TF1 *   quadratic  = new TF1("quadratic","[0]+[1]*x+[2]*x*x", 0, 1.4);
    quadratic->SetParNames ("Off-set","Linear factor", "Quadratic factor");
    quadratic->FixParameter(0,0);
    
    gr->Fit("quadratic");
    quadratic->Draw("SAME");
    gStyle->SetOptFit(111);

    pad2->cd();
    float diff[sizeof(x)/sizeof(x[0])];
    float diff_norm[sizeof(x)/sizeof(x[0])];
    for (Int_t i=0;i<sizeof(x)/sizeof(x[0]);i++) {
        diff[i] =y[i]-quadratic->Eval(x[i]);
        diff_norm[i]=diff[i]/y_err[i];
    }  
    float one_array[sizeof(x)/sizeof(x[0])]={1};
    
    TGraphErrors* gr2 = new TGraphErrors(sizeof(x)/sizeof(x[0]),x,diff_norm,nullptr,one_array);
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.5);
    gr2->SetMarkerStyle(1);
    gr2->GetXaxis()->SetLimits(0, 1.4);

    TF1 *   zero  = new TF1("zero","0*x", 0, 1.4);
    gr2->Draw("APE");
    zero->Draw("SAME");
    c_fit_quad->cd();
    //c_fit_quad->SaveAs(&("real_time_calibration/" + name_final + type + "_fit_quad.png")[0]);
    outfile->cd();
    c_fit_quad->Write();

    vector<double> off_set{quadratic->GetParameter(0),quadratic->GetParError(0)};
    vector<double> lin_factor{quadratic->GetParameter(1),quadratic->GetParError(1)};
    vector<double> quad_factor{quadratic->GetParameter(2),quadratic->GetParError(2)};
    vector<vector<double>> quad_params{off_set, lin_factor, quad_factor};

    return quad_params;
}



void calibration_triple(){
    // ordine del vector:
    // primo è per charge e secondo è per ampiezza
    // prima i 2 del NA |poi i 2 punti del compton edge Ne | e i 2 del picco del Ne
    map<string, vector<vector<float>>> pmts{
        {"pmt1_NA_c6_ext_coinc12_merc_metal_block_run8",{{70000, 79000, /* 140000, 152000,*/ 175000, 192000},
                {900, 1100, /*1830, 1920,*/ 2200, 2540}}},
        {"pmt2_NA_c6_ext_coinc12_merc_metal_block_run8",{{63000, 73000,/* 130000, 137000,*/ 157000, 176000},
                {880, 1070, /*1760, 1900,*/ 2080, 2440}}},
        {"pmt3_NA_c6_ext_coinc12_merc_metal_run6",{{11000, 14000, /*23500,26500,*/ 28000,34000},
                {160, 210,/* 350, 380,*/ 400, 470}}}
        };
    TFile *outfile= new TFile("triple_calib/calibration_plots.root", "RECREATE");
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
        auto cal_charge = fit_quad(name_pmt, "_charge", gaus_charge, outfile);
        auto cal_amp = fit_quad(name_pmt, "_amp", gaus_amp, outfile);

        for (int r=0;r<gaus_charge.size();r++){
            params_file<< "off_set "<<name_pmt<< "_run"<< to_string(r+1)<< "_charge = "<< cal_charge[0][0] << " +- " <<  cal_charge[0][1]<<"\n"; 
            params_file<< "cal_factor "<<name_pmt<< "_run"<< to_string(r+1)<< "_charge = "<< cal_charge[1][0] << " +- " <<  cal_charge[1][1]<<"\n"; 
            params_file<< "quad_factor "<<name_pmt<< "_run"<< to_string(r+1)<< "_charge = "<< cal_charge[2][0] << " +- " <<  cal_charge[2][1]<<"\n"; 
            
            params_file<< "off_set "<<name_pmt<< "_run"<< to_string(r+1)<< "_amp = "<< cal_amp[0][0] << " +- " <<  cal_amp[0][1]<<"\n"; 
            params_file<< "cal_factor "<<name_pmt<< "_run"<< to_string(r+1)<< "_amp = "<< cal_amp[1][0] << " +- " <<  cal_amp[1][1]<<"\n"; 
            params_file<< "quad_factor "<<name_pmt<< "_run"<< to_string(r+1)<< "_amp = "<< cal_amp[2][0] << " +- " <<  cal_amp[2][1]<<"\n"; 
           
            params_file <<"\n";
        } 
    }
    outfile->Close();
}

