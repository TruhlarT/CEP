#ifndef PID_h
#define PID_h

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

class Plot;

class PID{


public:

	PID(TFile* dataInput, TFile* fileOut, TString outnam="/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/", bool text = false, TString inputCuts="");
	~PID();

	void PlotHistogram();
	void SetBichsel(TH1* hist, Int_t color = 1, Float_t xMin = 0.2, Float_t xMax = 6, Float_t width = 4.0);
	void DrawBichsel();

	bool TEXT;


	TString output;
	TString cuts, cutsWithPrefix;

	TFile* data;
	TFile *fout;

	TH1F* tmpHist;
	TH1F* tmpHist2;
	TH1F* tmpHist3;
	TH1F* tmpHist4;

	TH2F* tmp2DHist;
	TH2F* tmp2DHist2;
	TH2F* tmp2DHist3;
	TH2F* tmp2DHist4;
	TH2F* tmp2DHist5;
	TH2F* tmp2DHist6;


};



#endif