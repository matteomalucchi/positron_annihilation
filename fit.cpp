#include <TGraphErrors.h>

void fit(){
    double x[4]={0.18, 0.66, 1.17, 1.33};
    double y[4]={0.31, 1.01, 1.75, 1.98};
    double y_err[4] = {0.07, 0.03, 0.04, 0.04};
    TGraphErrors* gr = new TGraphErrors(5,x,y,nullptr,y_err);
    gr->SetMarkerStyle(1);
    gr->Draw("APE");
    TF1 *   linear  = new TF1("linear","[0]+[1]*x",0, 1.4);
    linear->SetParameter(0, 2800);
    linear->SetParNames ("x_0 constant","Conversion factor");
    gr->Fit("linear");
    linear->Draw("SAME");
    gStyle->SetOptFit(111);
}

