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


// Na spetrum in run1 e run2 discrepancies

/*
void make_images(){
    TFile *f = new TFile("histograms/histograms_NA1.root");

    TH1F *histo1 = nullptr;
    TH1F *histo2 = nullptr;
    string name = "pmt1_Na_charge";
    string name_nice = "pmt1 Na charge";
    TCanvas *c = new TCanvas(&name[0],&name[0]);
    gStyle->SetOptStat("neou");


    f->GetObject("pmt1_NA_e6_100_run1_charge", histo1);
    f->GetObject("pmt1_NA_e6_100_run2_charge", histo2);

    histo1->SetNameTitle("pmt1 Na charge run1", &name_nice[0]);
    histo2->SetNameTitle("pmt1 Na charge run2", &name_nice[0]);
    histo1->GetXaxis()->SetTitle("Charge [a.u.]");
    histo2->GetXaxis()->SetTitle("Charge [a.u.]");
    histo1->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo2->GetBinWidth(1))))[0]);
    histo2->GetYaxis()->SetTitle(&("Entries / "+to_string(int(histo2->GetBinWidth(1))))[0]);

    histo1->SetLineColor(kBlue);
    histo2->SetLineColor(kGreen);

    TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",73000,84000);
    gaus1->SetParName(0, "const 1");
    gaus1->SetParName(1, "#mu 1");
    gaus1->SetParName(2, "#sigma 1");
    TF1 *gaus2 = new TF1("gaus2","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",73000,84000);
    gaus2->SetParName(0, "const 2");
    gaus2->SetParName(1, "#mu 2");
    gaus2->SetParName(2, "#sigma 2");

    gaus1->SetParameters(1000, 78000, 3000);
    gaus2->SetParameters(1000,  78000, 3000);

    histo2->Draw();
    TPaveStats *st = (TPaveStats*)histo2->FindObject("stats");
    st->SetX1NDC(0.61); //new x start position
    st->SetX2NDC(0.9); //new x end position
    st->SetY1NDC(0.6); //new x start position
    st->SetY2NDC(0.9); //new x end position

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
    st1->SetX2NDC(0.39); //new x end position
    st1->SetY1NDC(0.6); //new x start position
    st1->SetY2NDC(0.9); //new x end position

   auto legend = new TLegend(0.1,0.4,0.3,0.6);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(histo1,"pmt1 Na charge run1","l");
   legend->AddEntry(histo2,"pmt1 Na charge run2","l");
   legend->AddEntry("gaus1","gaus 1","l");
   legend->AddEntry("gaus2","gaus 2","l");
   legend->Draw();

    c->SaveAs(&("final_images/"+ name+".pdf")[0]);

}


// pmt1 spectra of co cs na
void make_images(){
    TFile *f = new TFile("histograms/histograms_cocsna.root");

    TH1F *histo1 = nullptr;
    TH1F *histo2 = nullptr;
    TH1F *histo3 = nullptr;
    string name = "pmt1_charge";
    string name_nice = "pmt1 charge";
    TCanvas *c = new TCanvas(&name[0],&name[0]);


    f->GetObject("pmt1_co_100_charge", histo1);
    f->GetObject("pmt1_cs_100_charge", histo2);
    f->GetObject("pmt1_NA_e6_100_run2_charge", histo3);


    histo1->SetNameTitle("pmt1 Co charge", &name_nice[0]);
    histo2->SetNameTitle("pmt1 Cs charge", &name_nice[0]);
    histo3->SetNameTitle("pmt1 Na charge", &name_nice[0]);


    //gStyle->SetOptStat("neou");

    histo1->GetXaxis()->SetTitle("Charge [a.u.]");
    histo2->GetXaxis()->SetTitle("Charge [a.u.]");
    histo3->GetXaxis()->SetTitle("Charge [a.u.]");
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
    legend->AddEntry(histo2,"pmt1 Cs charge","l");
    legend->AddEntry(histo3,"pmt1 Na charge","l");
    legend->AddEntry(histo1,"pmt1 Co charge","l");
    //legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
    legend->Draw();

    c->SaveAs(&("final_images/"+ name+".pdf")[0]);

}
*/


void make_images(){
    TFile *f = new TFile("histograms/histograms_cocsna.root");

    TH1F *histo1 = nullptr;
    TH1F *histo2 = nullptr;
    TH1F *histo3 = nullptr;
    string name = "pmt1_charge";
    string name_nice = "pmt1 charge";
    TCanvas *c = new TCanvas(&name[0],&name[0]);


    f->GetObject("pmt1_co_100_charge", histo1);
    f->GetObject("pmt1_cs_100_charge", histo2);
    f->GetObject("pmt1_NA_e6_100_run2_charge", histo3);


    histo1->SetNameTitle("pmt1 Co charge", &name_nice[0]);
    histo2->SetNameTitle("pmt1 Cs charge", &name_nice[0]);
    histo3->SetNameTitle("pmt1 Na charge", &name_nice[0]);


    //gStyle->SetOptStat("neou");

    histo1->GetXaxis()->SetTitle("Charge [a.u.]");
    histo2->GetXaxis()->SetTitle("Charge [a.u.]");
    histo3->GetXaxis()->SetTitle("Charge [a.u.]");
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
    legend->AddEntry(histo2,"pmt1 Cs charge","l");
    legend->AddEntry(histo3,"pmt1 Na charge","l");
    legend->AddEntry(histo1,"pmt1 Co charge","l");
    //legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
    legend->Draw();

    c->SaveAs(&("final_images/"+ name+".pdf")[0]);

}