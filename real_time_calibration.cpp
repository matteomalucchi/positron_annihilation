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


using namespace std;

void real_time_calibration(){
    TFile *f = new TFile("histograms/histograms_RealTimeCalibration.root");

    map <string, string> pair_names ={ 
        /*{"pmt1_NA+cs_e6_100_run1", 
                        "pmt2_NA+cs_e6_100_run1",
                        "pmt3_NA+cs_e6_30_run1",
                        "pmt1_NA+cs+co_e6_100_run1",
                        "pmt2_NA+cs+co_e6_100_run1",
                        "pmt3_NA+cs+co_e6_30_run1", */
                        {"pmt1_NA+cs_e6_500_run2", "pmt1_NA+cs+co_e6_500_run2"}

                        /*"pmt2_NA+cs_e6_500_run2",
                        "pmt3_NA+cs_e6_100_run2",
                        "pmt2_NA+cs+co_e6_500_run2",
                        "pmt3_NA+cs+co_e6_100_run2"*/};

    for (const auto &pair_name : pair_names){  
        const auto name_a = pair_name.first;
        const auto name_b = pair_name.second;
        



}