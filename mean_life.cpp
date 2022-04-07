#include<TFile.h>
#include<TNtuple.h>

void mean_life(){
    TFile *f= new TFile("triple/ntuple_triple_lin.root");
    TNtuple* n_run3;
    f->GetObject("run3", n_run3);
    TH1F* h= new TH1F("h", "h", 64, 0, 2);
    for (auto ievt: ROOT::TSeqI(n_run3->GetEntries())) {
        n_run3->GetEntry(ievt);
        auto catm = n_run3->GetArgs(); // Get a row
        h->Fill(catm[0]);
    }
    h->Draw();
}