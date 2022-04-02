void p(TH1F* h1) {
   //gStyle->SetOptStat(0);
   auto c1 = new TCanvas("c1", "fit residual simple");
   TF1 *gaus1 = new TF1("gaus1","[0]*exp(-0.5*pow(((x-[1])/[2]),2))",-1,1);
    gaus1->SetParameters(200,0,1);

   h1->Fit("gaus1","R", "SAME");
   //h1->GetXaxis()->SetTitle("x");
   c1->Clear(); // Fit does not draw into correct pad
   auto rp1 = new TRatioPlot(h1);
   rp1->Draw();
   rp1->GetLowerRefYaxis()->SetTitle("ratio");
   rp1->GetUpperRefYaxis()->SetTitle("entries");
   c1->Update();
}
void ratioplot2(){
   auto h = new TH1F("a", "a", 50, -6, 6);
   h->FillRandom("gaus", 2000);

    p(h);
}