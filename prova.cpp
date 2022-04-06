void prova(){
    /*
    float m_ave=0.5109, m_ave_err=0.003;
    int dimension =6;
    float y[dimension], y_err[dimension];
    for (int d=0;d<dimension;d++){
        y[d]=d+1;
        y_err[d]=1;
    }

    TCanvas *c_mass = new TCanvas("c", "c",200,10,600,400);

    float ave[dimension], ave_err[dimension];
    for(int p=0;p<dimension;p++){
        ave[p]=m_ave;
        ave_err[p]=m_ave_err;
    }


    auto ge = new TGraphErrors(dimension, ave, y, ave_err);
    ge->GetXaxis()->SetLimits(0.4,0.6);
    ge->SetFillColor(4);
    ge->SetFillStyle(3010);
    ge->Draw("a4");

    auto meanline=new TLine(m_ave, 0.5, m_ave, dimension+0.5);
    meanline->SetLineColor(4);
        meanline->SetLineStyle(9);

    //meanline->Draw("SAME");*/
    
   auto c41 = new TCanvas("c41","c41",200,10,600,400);
   double x[] = {0, 1, 2, 3, 4};
   double y[] = {0, 2, 4, 1, 3};
   double ex[] = {0.1, 0.2, 0.3, 0.4, 0.5};
   double ey[] = {1, 0.5, 1, 0.5, 1};
   auto ge = new TGraphErrors(5, x, y, ex, ey);
   ge->SetFillColor(4);
   ge->SetFillStyle(3010);
   ge->Draw("a3");


}