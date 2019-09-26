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

#include "BasicPlots.h"
#include "PID.h"
#include "trackQuality.h"

using namespace std;

//_____________________________________________________________________________
int main(int argc, char** argv) {
	//open output file

	bool text = false;

	if(argc != 2){
		cout<<"Error: wrong input"<<endl;
		cout<<"You should do: ./Analysis DataSet"<<endl;
		return 1;
	}

	//const string name = "StRP_production_0000.root";
	const string& DataSet = argv[1];
	const string name = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/" + DataSet + ".root";
	const string output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/" + DataSet + "/";
	

	BasicPlots Plots(name , output, text);
	Plots.PlotHistogram();


	PID PIDPlots(name, output, text);
	PIDPlots.PlotHistogram();

	trackQuality TrackPlots(name, output, text);
	TrackPlots.PlotHistogram();

	cout<<"Ending Analysis... GOOD BYE!"<<endl;
	return 0;
}//main

