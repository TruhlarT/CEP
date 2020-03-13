#ifndef trackQuality_h
#define trackQuality_h


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
#include "TEllipse.h"

using namespace std;

class Plot;

class trackQuality{


public:

	trackQuality(TFile* dataInput, TFile* fileOut, TString outnam="/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/", bool text = false, TString inputCuts="");
	~trackQuality();

	void PlotHistogram();



	bool TEXT;

	TH1F* tmpHist;
	TH1F* tmpHist2;
	TH1F* tmpHist3;

	TH1F* histSignal;
	TH1F* histBackground;

	TH2F* tmp2DHist;
	TH2F* tmp2DHist2;
	TH2F* tmp2DHist3;

	TString output;
	TString cuts, cutsWithPrefix;
	TString dataLabel;

	TFile* data;
	TFile* fout;


};

#endif