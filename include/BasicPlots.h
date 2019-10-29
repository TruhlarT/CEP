#ifndef BasicPlots_h
#define BasicPlots_h


#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TTree.h"
#include "TLatex.h"
#include "TF1.h"
#include "TFile.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include <fstream>

using namespace std;

class Plot;

class BasicPlots{

public:

	BasicPlots(TFile* dataInput, TFile* fileOut, TString outnam="/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/", bool text = false, TString inputCuts="");
	~BasicPlots();

	void PlotHistogram();

	TList *mHistList; // list of output histograms

	TH1D* hMissingPtTPC;
	TH1D* hMissingPtTOF;
	TH1D* hMissingPtQ0;
	TH1D* hMissingPtExc;
	TH1F* hCuts;
	TH1F* tmpHist;

	TLegend* leg1;

	TString output;
	TString cuts, cutsWithPrefix;

	bool TEXT;
	bool bcground;
	
	TFile* data;
	TFile* fout;

	Float_t siz;
	Int_t font;

};


#endif