#ifndef RPplots_h
#define RPplots_h


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

class RPplots{

public:

	RPplots(TFile* dataInput, TFile* fileOut, TString outnam="/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/", bool text = false, TString inputCuts="");
	~RPplots();

	void PlotHistogram();


	bool TEXT;

	TH1F* histSignal;
	TH1F* histBackground;
	TH1F* tmpHist;

	TH2F* tmp2DHist;
	TH2F* hist2DSignal;
	TH2F* hist2DBackground;


	TLegend* leg1;

	TString output;
	TString cuts;

	TFile* data;
	TFile* fout;

};


#endif