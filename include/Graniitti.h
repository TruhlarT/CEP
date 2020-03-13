#ifndef Graniitti_h
#define Graniitti_h

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

class Graniitti{


public:

	Graniitti(TFile* dataInput, TFile* graniittiInput, TFile* fileOut, TString outnam="/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/", bool text = false, TString inputCuts="");
	~Graniitti();

	void PlotHistogram();
	void PlotComparison();

	bool TEXT;
	enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};

	TString output;
	TString cuts, cutsWithPrefix;
	TString dataLabel;
	
	TFile* data;
	TFile* graniitti;
	TFile *fout;

	TH1F* tmpHist;
	TH1F* tmpHist2;

};



#endif