// include headers
#include "libreries.h"

using namespace std;


void BcgEstimation()
{
    //TFile* data = TFile::Open("Data/BcgEstimation.root", "read");
    //TH1D *hPtMiss = (TH1D*)data->Get("PtMiss");

    TFile* data = TFile::Open("Data/nonExBackground.root", "read");
    TH2D *hTmp = (TH2D*)data->Get("InvMassVsPtMissPion");
    //cout<<hTmp<<"\n";
    TH1D *hPtMiss = (TH1D*)hTmp->ProjectionY();

    gStyle->SetOptStat("1111");
    gStyle->SetOptFit(1); 

    TF1 *f1 = new TF1("f1", "[0]*x + [1]*x*x",0.0,0.5);
    f1->SetParameters(10,5);
    hPtMiss->Fit("f1","","",0.3,0.5);
    hPtMiss->SetTitle(" ; p_{T}^{miss} [GeV];Number of events");
    hPtMiss->Draw();
    f1->Draw("same");

    TCanvas *c2 = new TCanvas("cCanvas2","cCanvas2",800,700);
    TH1F* hSigToBcg = new TH1F("SigToBcg", "Signal to Background ratio", 280, 59, 204);
    double bcg, signal, ratio, cut;
    for (int i = 1; i <= hSigToBcg->GetXaxis()->GetNbins(); ++i)
    {
        cut = (60 + i*0.5)/1000;
        bcg = f1->Integral(0.0, cut);
        signal = hPtMiss->Integral(1, i) - bcg;
        ratio = signal / bcg;
        //cout<<"S: "<<signal<<" B: "<<bcg<<" R: "<<ratio<<"\n";
        hSigToBcg->SetBinContent(i,ratio);
    }
    hSigToBcg->Draw();
    cout<<"Max signal to background ratio: "<< hSigToBcg->GetMaximum()<<" for cut "<< hSigToBcg->GetXaxis()->GetBinCenter(hSigToBcg->GetMaximumBin())<<"\n";
}



