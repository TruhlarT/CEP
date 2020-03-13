 // c++ headers
#include <iostream>
#include <utility>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath> 
#include <cstdlib>
#include <unistd.h>

// ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include <TH2.h> 
#include <TF1.h> 
#include <TF2.h> 
#include <THStack.h> 
#include <TStyle.h> 
#include <TGraph.h> 
#include <TGraph2D.h> 
#include <TGraphErrors.h> 
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TGaxis.h> 
#include <TString.h> 
#include <TColor.h> 
#include <TLine.h> 
#include <TExec.h> 
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TLatex.h> 
#include <TMath.h>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TDirectory.h"


using namespace std;


TFile* data;

void savePlots()
{
	TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/part1New/StRP.root";
    TString path = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/part1New/Figs/";

	data = TFile::Open(input, "read");
	if (!data)
	{
		cout<<"Error: cannot open "<<input<<endl;
		return;
	}


    TString name[] = {TString{"cuts"}, TString{"ptMiss"}, TString{"qualityCut"}, TString{"zVertex"}, 
                    TString{"NHitsFit"}, TString{"NHits"}, TString{"DCAz"}, TString{"DCAxy"}, 
                    TString{"eta"}, TString{"t"}, TString{"pid"}, TString{"m2"}, TString{"pionsUncorr"}, 
                    TString{"pions"}, TString{"pionsUncorrEl"}, TString{"pionsEl"}, TString{"pionsUncorrInel"}, 
                    TString{"pionsInel"}, TString{"kaonsUncorr"}, TString{"kaons"}, TString{"kaonsUncorrEl"}, 
                    TString{"kaonsEl"}, TString{"kaonsUncorrInel"}, TString{"kaonsInel"}, TString{"protonsUncorr"}, 
                    TString{"protons"}, TString{"protonsUncorrEl"}, TString{"protonsEl"}, 
                    TString{"protonsUncorrInel"}, TString{"protonsInel"},
                    TString{"ratioPlotPionEl+Inel"}, TString{"ratioPlotPionEl"},
                    TString{"ratioPlotPionInel"}, TString{"ratioPlotKaonEl+Inel"},
                    TString{"ratioPlotKaonEl"}, TString{"ratioPlotKaonInel"},
                    TString{"ratioPlotProtonEl+Inel"}, TString{"ratioPlotProtonEl"},
                    TString{"ratioPlotProtonInel"},
                    TString{"granInvMassEl+InelPion"}, 
                    TString{"granInvMassElPion"}, 
                    TString{"granInvMassInelPion"}, 
                    TString{"granInvMassEl+InelKaon"}, 
                    TString{"granInvMassElKaon"}, 
                    TString{"granInvMassInelKaon"}, 
                    TString{"granInvMassEl+InelProton"}, 
                    TString{"granInvMassElProton"}, 
                    TString{"granInvMassInelProton"},
                    TString{"granPairRapEl+InelPion"}, 
                    TString{"granPairRapElPion"}, 
                    TString{"granPairRapInelPion"}, 
                    TString{"granPairRapEl+InelKaon"}, 
                    TString{"granPairRapElKaon"}, 
                    TString{"granPairRapInelKaon"}, 
                    TString{"granPairRapEl+InelProton"}, 
                    TString{"granPairRapElProton"}, 
                    TString{"granPairRapInelProton"},
                    TString{"granPhiEl+InelPion"}, 
                    TString{"granPhiElPion"}, 
                    TString{"granPhiInelPion"}, 
                    TString{"granPhiEl+InelKaon"}, 
                    TString{"granPhiElKaon"}, 
                    TString{"granPhiInelKaon"}, 
                    TString{"granPhiEl+InelProton"}, 
                    TString{"granPhiElProton"}, 
                    TString{"granPhiInelProton"}, 
                    TString{"grantSumEl+InelPion"}, 
                    TString{"grantSumElPion"}, 
                    TString{"grantSumInelPion"}, 
                    TString{"grantSumEl+InelKaon"}, 
                    TString{"grantSumElKaon"}, 
                    TString{"grantSumInelKaon"}, 
                    TString{"grantSumEl+InelProton"}, 
                    TString{"grantSumElProton"}, 
                    TString{"grantSumInelProton"}, 
                    TString{"tSumEl+InelPion"}, 
                    TString{"tSumElPion"}, 
                    TString{"tSumInelPion"}, 
                    TString{"tSumEl+InelKaon"}, 
                    TString{"tSumElKaon"}, 
                    TString{"tSumInelKaon"}, 
                    TString{"tSumEl+InelProton"}, 
                    TString{"tSumElProton"}, 
                    TString{"tSumInelProton"}, 
                    TString{"PhiEl+InelPion"}, 
                    TString{"PhiElPion"}, 
                    TString{"PhiInelPion"}, 
                    TString{"PhiEl+InelKaon"}, 
                    TString{"PhiElKaon"}, 
                    TString{"PhiInelKaon"}, 
                    TString{"PhiEl+InelProton"}, 
                    TString{"PhiElProton"}, 
                    TString{"PhiInelProton"},
                    TString{"PairRapEl+InelPion"}, 
                    TString{"PairRapElPion"}, 
                    TString{"PairRapInelPion"}, 
                    TString{"PairRapEl+InelKaon"}, 
                    TString{"PairRapElKaon"}, 
                    TString{"PairRapInelKaon"}, 
                    TString{"PairRapEl+InelProton"}, 
                    TString{"PairRapElProton"}, 
                    TString{"PairRapInelProton"}
                };
    TString nameSource[] = {TString{"BasicPlots/CutsFlow"}, TString{"BasicPlots/hMissingPt"}, 
                    TString{"Cuts/CutsFlow"}, TString{"El+Inel/trackQuality/zVertex"}, 
                    TString{"Cuts/trckQ_0/NhitsFit"}, TString{"Cuts/trckQ_1/NhitsDEdx"}, 
                    TString{"Cuts/trckQ_2/DcaZ"}, TString{"Cuts/trckQ_3/DcaXY"}, 
                    TString{"Cuts/trckQ_4/Eta"}, TString{"Cuts/trckQ_5/t"}, 
                    TString{"GoldenEvents/El+Inel/PID/mSquaredVschiPairPion"}, 
                    TString{"GoldenEvents/El+Inel/PID/mSquared"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassPionEl+Inel"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassPionEl+Inel"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassPionEl"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassPionEl"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassPionInel"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassPionInel"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassKaonEl+Inel"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassKaonEl+Inel"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassKaonEl"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassKaonEl"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassKaonInel"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassKaonInel"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassProtonEl+Inel"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassProtonEl+Inel"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassProtonEl"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassProtonEl"}, 
                    TString{"MoneyPlots/Uncorrected/uncorrInvMassProtonInel"}, 
                    TString{"MoneyPlots/Corrected/corrInvMassProtonInel"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotPionEl+Inel"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotPionEl"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotPionInel"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotKaonEl+Inel"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotKaonEl"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotKaonInel"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotProtonEl+Inel"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotProtonEl"},
                    TString{"MoneyPlots/RatioPlots/ratioPlotProtonInel"},
                    TString{"Graniitti/El+Inel/invMass_statePion"}, 
                    TString{"Graniitti/El/invMass_statePion"}, 
                    TString{"Graniitti/Inel/invMass_statePion"}, 
                    TString{"Graniitti/El+Inel/invMass_stateKaon"}, 
                    TString{"Graniitti/El/invMass_stateKaon"}, 
                    TString{"Graniitti/Inel/invMass_stateKaon"}, 
                    TString{"Graniitti/El+Inel/invMass_stateProton"}, 
                    TString{"Graniitti/El/invMass_stateProton"}, 
                    TString{"Graniitti/Inel/invMass_stateProton"},
                    TString{"Graniitti/El+Inel/rap_statePion"}, 
                    TString{"Graniitti/El/rap_statePion"}, 
                    TString{"Graniitti/Inel/rap_statePion"}, 
                    TString{"Graniitti/El+Inel/rap_stateKaon"}, 
                    TString{"Graniitti/El/rap_stateKaon"}, 
                    TString{"Graniitti/Inel/rap_stateKaon"}, 
                    TString{"Graniitti/El+Inel/rap_stateProton"}, 
                    TString{"Graniitti/El/rap_stateProton"}, 
                    TString{"Graniitti/Inel/rap_stateProton"}, 
                    TString{"Graniitti/El+Inel/phiPion"}, 
                    TString{"Graniitti/El/phiPion"}, 
                    TString{"Graniitti/Inel/phiPion"}, 
                    TString{"Graniitti/El+Inel/phiKaon"}, 
                    TString{"Graniitti/El/phiKaon"}, 
                    TString{"Graniitti/Inel/phiKaon"}, 
                    TString{"Graniitti/El+Inel/phiProton"}, 
                    TString{"Graniitti/El/phiProton"}, 
                    TString{"Graniitti/Inel/phiProton"}, 
                    TString{"Graniitti/El+Inel/tSumPion"}, 
                    TString{"Graniitti/El/tSumPion"}, 
                    TString{"Graniitti/Inel/tSumPion"}, 
                    TString{"Graniitti/El+Inel/tSumKaon"}, 
                    TString{"Graniitti/El/tSumKaon"}, 
                    TString{"Graniitti/Inel/tSumKaon"}, 
                    TString{"Graniitti/El+Inel/tSumProton"}, 
                    TString{"Graniitti/El/tSumProton"}, 
                    TString{"Graniitti/Inel/tSumProton"}, 
                    TString{"MoneyPlots/Uncorrected/tSumPionEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/tSumPionEl"},
                    TString{"MoneyPlots/Uncorrected/tSumPionInel"},
                    TString{"MoneyPlots/Uncorrected/tSumKaonEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/tSumKaonEl"},
                    TString{"MoneyPlots/Uncorrected/tSumKaonInel"},
                    TString{"MoneyPlots/Uncorrected/tSumProtonEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/tSumProtonEl"},
                    TString{"MoneyPlots/Uncorrected/tSumProtonInel"},
                    TString{"MoneyPlots/Uncorrected/phiRPPionEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/phiRPPionEl"},
                    TString{"MoneyPlots/Uncorrected/phiRPPionInel"},
                    TString{"MoneyPlots/Uncorrected/phiRPKaonEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/phiRPKaonEl"},
                    TString{"MoneyPlots/Uncorrected/phiRPKaonInel"},
                    TString{"MoneyPlots/Uncorrected/phiRPProtonEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/phiRPProtonEl"},
                    TString{"MoneyPlots/Uncorrected/phiRPProtonInel"},
                    TString{"MoneyPlots/Uncorrected/pairRapPionEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/pairRapPionEl"},
                    TString{"MoneyPlots/Uncorrected/pairRapPionInel"},
                    TString{"MoneyPlots/Uncorrected/pairRapKaonEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/pairRapKaonEl"},
                    TString{"MoneyPlots/Uncorrected/pairRapKaonInel"},
                    TString{"MoneyPlots/Uncorrected/pairRapProtonEl+Inel"},
                    TString{"MoneyPlots/Uncorrected/pairRapProtonEl"},
                    TString{"MoneyPlots/Uncorrected/pairRapProtonInel"},
                }; 

    TCanvas *cCanvas; 

    int size = (sizeof(name)/sizeof(*name));
    int sizeSource = (sizeof(nameSource)/sizeof(*nameSource));

    if(size != sizeSource)
    {
        cout<<"Size of names is different from size of name source..."<<endl;
        cout<<"Size / sizeSource: "<<size<<" / "<<sizeSource<<endl;
        cout<<"Exiting with error"<<endl;
        return;
    }

    for (int i = 0; i < size; ++i)
    {
        cCanvas = (TCanvas*)data -> Get(nameSource[i]);
        if(!cCanvas)
        {
            cout<<"Error: cannot get "<<nameSource[i]<<endl;
            cout<<"Exiting with error"<<endl;
            return;
        }
        cCanvas->SaveAs(path + name[i] + ".pdf");
    }
    data->Close();
    
}
