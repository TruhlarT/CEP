#include "RPplots.h"
#include "Plot.h"


//_____________________________________________________________________________
RPplots::RPplots(TFile* dataInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	output = outnam;
	cuts = inputCuts;

	data = dataInput;
	fout = fileOut;

	TEXT = text;

	cout << "trackQuality::trackQuality() called" << endl;
}//RPplots

//_____________________________________________________________________________
RPplots::~RPplots()
{
  //destructor

  cout << "RPplots::~RPplots() destructor called" << endl;

}//~RPplots


//_____________________________________________________________________________
void RPplots::PlotHistogram() {

	if (!data){
		cout<<"Error: cannot open input file"<<endl;
		return;
    }

	TTree* tree = dynamic_cast<TTree*>( data->Get("recTree") );
	TTree* treeBack = dynamic_cast<TTree*>( data->Get("Background") );
	
	if (!tree || !treeBack){
		cout<<"Error: cannot open one of the TTree"<<endl;
		return;
    }

	if (!fout){
		cout<<"Error: cannot open output file"<<endl;
		return;
    }

    TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
    gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
	gPad->SetTickx();
	gPad->SetTicky();  
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetLineWidth(2);      //axis line
	gStyle->SetFrameLineWidth(2); //frame line

	tree->UseCurrentStyle();
	treeBack->UseCurrentStyle();

	Plot tool;
	TString variable, variable2, cutsOption;
	Int_t nBins, nBins2, nInputs;
	Float_t min, max, min2, max2;

//////////////////////////////////////////////////
// Plot XYEastCor vertex
    gPad->SetLogz();
    variable = "yCorrelationsEast";
    variable2 = "xCorrelationsEast";
    cutsOption = cuts;
    nBins = 100;
    min = -1;
    max = 1.4;
    nBins2 = 100;
    min2 = -1;
    max2 = 1.4;

    tree->Draw(variable + ":" + variable2 + ">>Sig(" + nBins "," + min + "," + max + "," + nBins2 "," + min2 + "," + max2 +")",cuts,"colz");
    hist2DSignal = (TH2F*)gPad->GetPrimitive("Sig");
    hist2DSignal->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");
    tool.SetGraphStyle(hist2DSignal,4,20,1,4,1,1,0.9,0.7);
    hist2DSignal->Draw("colz");
    tool.DrawText(hist2DSignal,0,true,0.62,0.81,0.76,0.94);
    tool.DrawTextStar(hist2DSignal);
    if(TEXT){
        const Int_t n = 200;
        Double_t x[n], y[n];
        Double_t tmp;
        for(int i = 0; i < n; ++i){
            x[i] = -0.036 + (0.152*i)/n;
            tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
            y[i] = sqrt(abs(tmp));
            if(tmp<0)
                y[i] = - y[i];
        }
        TGraph* gr = new TGraph(n,x,y);
        gr->SetLineWidth(4);
        gr->Draw("same");

        for(int i = 0; i < n; ++i){
            x[i] = -0.22 + (0.185*i)/n;
            tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
            y[i] = sqrt(abs(tmp));
            if(tmp<0)
                y[i] = - y[i];
        }
        TGraph* gr1 = new TGraph(n,x,y);
        gr1->SetLineWidth(4);
        gr1->Draw("same");

        TLine *left02 = new TLine(-0.219,-0.2,0.116,-0.2);
        tool.SetLineStyle(left02);
        left02->Draw("same");


        TLine *left01 = new TLine(-0.219,0.2,0.116,0.2);
        tool.SetLineStyle(left01);
        left01->Draw("same");

        for(int i = 0; i < n; ++i){
            x[i] = -0.036 + (0.152*i)/n;
            tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
            y[i] = sqrt(abs(tmp));
        }
        TGraph* gr2 = new TGraph(n,x,y);
        gr2->SetLineWidth(4);
        gr2->Draw("same");

        for(int i = 0; i < n; ++i){
            x[i] = -0.22 + (0.185*i)/n;
            tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
            y[i] = sqrt(abs(tmp));
        }
        TGraph* gr3 = new TGraph(n,x,y);
        gr3->SetLineWidth(4);
        gr3->Draw("same");
    }
    cCanvas2D->Update();
    cCanvas2D->Write(variable + "VS" + variable2);
//////////////////////////////////////////
//////////////////////////////////////////////////
// Plot XYEastCor vertex
    gPad->SetLogz();
    variable = "yCorrelationsWest";
    variable2 = "xCorrelationsWest";
    cutsOption = cuts;
    nBins = 100;
    min = -1;
    max = 1.4;
    nBins2 = 100;
    min2 = -1;
    max2 = 1.4;

    tree->Draw(variable + ":" + variable2 + ">>Sig(" + nBins "," + min + "," + max + "," + nBins2 "," + min2 + "," + max2 +")",cuts,"colz");
    hist2DSignal = (TH2F*)gPad->GetPrimitive("Sig");
    hist2DSignal->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");
    tool.SetGraphStyle(hist2DSignal,4,20,1,4,1,1,0.9,0.7);
    hist2DSignal->Draw("colz");
    tool.DrawText(hist2DSignal,0,true,0.62,0.81,0.76,0.94);
    tool.DrawTextStar(hist2DSignal);
    if(TEXT){
        const Int_t n = 200;
        Double_t x[n], y[n];
        Double_t tmp;
        for(int i = 0; i < n; ++i){
            x[i] = -0.036 + (0.152*i)/n;
            tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
            y[i] = sqrt(abs(tmp));
            if(tmp<0)
                y[i] = - y[i];
        }
        TGraph* gr = new TGraph(n,x,y);
        gr->SetLineWidth(4);
        gr->Draw("same");

        for(int i = 0; i < n; ++i){
            x[i] = -0.22 + (0.185*i)/n;
            tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
            y[i] = sqrt(abs(tmp));
            if(tmp<0)
                y[i] = - y[i];
        }
        TGraph* gr1 = new TGraph(n,x,y);
        gr1->SetLineWidth(4);
        gr1->Draw("same");

        TLine *left02 = new TLine(-0.219,-0.2,0.116,-0.2);
        tool.SetLineStyle(left02);
        left02->Draw("same");


        TLine *left01 = new TLine(-0.219,0.2,0.116,0.2);
        tool.SetLineStyle(left01);
        left01->Draw("same");

        for(int i = 0; i < n; ++i){
            x[i] = -0.036 + (0.152*i)/n;
            tmp = (x[i] + 0.3)*(x[i]+0.3) - 0.46*0.46;
            y[i] = sqrt(abs(tmp));
        }
        TGraph* gr2 = new TGraph(n,x,y);
        gr2->SetLineWidth(4);
        gr2->Draw("same");

        for(int i = 0; i < n; ++i){
            x[i] = -0.22 + (0.185*i)/n;
            tmp = (x[i] + -0.15)*(x[i]-0.15) - 0.42*0.42;
            y[i] = sqrt(abs(tmp));
        }
        TGraph* gr3 = new TGraph(n,x,y);
        gr3->SetLineWidth(4);
        gr3->Draw("same");
    }
    cCanvas2D->Update();
    cCanvas2D->Write(variable + "VS" + variable2);
//////////////////////////////////////////


}//RPplots::PlotHistogram


