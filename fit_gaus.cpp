#include <list>
#include <string>


#include <TFile.h>


void fit_gaus (){
    TFile *f = new TFile("histograms.root");

    std::map<std::string, vector<float>> samples = {
    {"pmt1_co_100_690", {165000, 185000, 185000, 210000}}/*,
    {"pmt2_co_100_851", 0.077 / 1499064.0 * scaleFactorZZTo4l * integratedLuminosity},
    {"pmt1_na_100_690", 0.077 / 1499093.0 * scaleFactorZZTo4l * integratedLuminosity},
    {"pmt2_na_100_851", 0.18 / 1497445.0 * scaleFactorZZTo4l * integratedLuminosity},
    {"pmt1_cs_100_690", 1.0},
    {"pmt2_cs_100_851", 1.0}*/
    };
    TH1F *h_charge = nullptr;
    TH1F *h_amp = nullptr;

    for (const auto &sample : samples) {
        const auto name = sample.first;
        const auto range = sample.second;
        std::cout << range[0] << std::endl;

        f->GetObject(&(name + "_charge")[0], h_charge);
        TCanvas *c_charge = new TCanvas(&(name + "_charge")[0], &(name + "_charge")[0]);
        TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[0],range[1]);
        TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",range[2],range[3]);
        //gaus1->SetParNames ("const","Conversion factor");
        gaus1->SetParameters(00, abs(range[0]+range[1])/2, abs(range[1]-range[0]));
        gaus2->SetParameters(00, abs(range[2]+range[3])/2, abs(range[3]-range[2]));
        h_charge->Fit("gaus1","R", "SAME");
        gaus1->Draw("SAME");
        h_charge->Fit("gaus2","R", "SAME"); 
        //gaus2->Draw("SAME");
        h_charge->Draw("SAME");
        gStyle->SetOptFit(111111);

        //f->GetObject(&(name + "_amp")[0], h_amp);
        //h_amp->Draw();
    }
}