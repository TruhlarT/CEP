#ifndef Plot_h
#define Plot_h


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

class Plot{

public:

	Plot();
	~Plot();



	void SetGraphStyle(TH1* hist, Int_t markColor = 4, Int_t markStyle = 20, Int_t markSize = 1, 
						Int_t lineColor = 4, Int_t lineStyle = 1, Float_t lineWidth = 1, Float_t xOffSet = 0.9, Float_t yOffSet = 1.1);

	void SetMarkerStyle(TH1* hist, Int_t markColor = 4, Int_t markStyle = 20, Int_t markSize = 1, 
							Int_t lineColor = 4, Int_t lineStyle = 1, Float_t lineWidth = 1);

	void DrawText(TH1* hist, Bool_t pair = false, Float_t xMin = 0.74, Float_t yMin = 0.82, Float_t xMax = 0.81, Float_t yMax = 0.96);

	void SetLegendStyle(TLegend* leg1);
	void SetTextStyle(TPaveText* text);
	void SetLineStyle(TLine* line, Int_t style = 1, Int_t color = 1, Int_t width = 4);

	Float_t siz;
	Int_t font;

};


#endif
