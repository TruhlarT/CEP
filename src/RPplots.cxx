#include "RPplots.h"
#include "Plot.h"


//_____________________________________________________________________________
RPplots::RPplots(TFile* dataInput, TFile* fileOut, TString outnam, bool text, TString inputCuts)
{
	//constructor
	output = outnam;
	cuts = inputCuts;
    if(inputCuts != "")
    {
        cutsWithPrefix = " && " + inputCuts;
    }
    else
    {
        cutsWithPrefix="";
    }

	data = dataInput;
	fout = fileOut;

	TEXT = text;

    if(inputCuts.Contains("!elastic"))
    {
        dataLabel = "Inel";
    }
    else if(inputCuts.Contains("elastic"))
    {
        dataLabel = "El";
    }
    else
    {
        dataLabel = "El+Inel";
    }

	cout << "RPplots::RPplots() called" << endl;
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

    TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas",800,700);
    gPad->SetMargin(0.12,0.02,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky();  
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line

	Plot tool;
	TString variable, variable2, cutsOption;
	Int_t nBins, nBins2, nInputs;
	Float_t min, max, min2, max2;

//////////////////////////////////////////////////
// Plot t 
    variable = "t";
    nBins = 100;
    min = -2.0;
    max = 0.0;
    treeBack->Draw(variable +"East>>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",cuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    treeBack->Draw(variable +"West>>" + variable +"Bcg2(" + nBins + "," + min + "," + max + ")",cuts);
    tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg2");
    histBackground->Add(tmpHist);
    tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);


    tree->Draw(variable +"East>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",cuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
    tree->Draw(variable +"West>>" + variable +"Sig2(" + nBins + "," + min + "," + max + ")",cuts);
    tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
    histSignal->Add(tmpHist);   
    histSignal->SetTitle(" ; t [GeV^{2}]; Number of tracks");
    tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
    tool.SetMarkerStyle(histSignal);
    histSignal->Draw("E");
    tool.DrawText(histSignal,0,false,0.15,0.74,0.28,0.9,12);
    tool.DrawTextStar(histSignal,3);
    histBackground->Draw("ESAME");

    if(TEXT){
        TLine *left4 = new TLine(-1.0,0,-1.0,tmpHist->GetMaximum()/2);
        tool.SetLineStyle(left4,10,1,4);
        left4->Draw("same");

        TLine *left5 = new TLine(-0.12,0,-0.12  ,tmpHist->GetMaximum()/2);
        tool.SetLineStyle(left5,10,1,4);
        left5->Draw("same");     
    }
                                
    TLegend* leg1 = new TLegend(0.15, 0.65, 0.28, 0.74);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(histSignal,dataLabel + " (unlike-sign pairs)","p");
    leg1->AddEntry(histBackground,dataLabel + " (like-sign pairs)","p");
    leg1->Draw("same");


    cCanvas->Update();
    cCanvas->Write(variable);
//////////////////////////////////////
// Plot t 
    variable = "t";
    nBins = 100;
    min = 0.0;
    max = 4.0;
    treeBack->Draw("TMath::Abs( " + variable +"East + " + variable +"West)>>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",cuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

    tree->Draw("TMath::Abs( " + variable +"East + " + variable +"West)>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",cuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
    histSignal->SetTitle(" ; |t_{1} + t_{2}| [GeV^{2}]; Number of tracks");
    tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
    tool.SetMarkerStyle(histSignal);
    histSignal->Draw("E");
    tool.DrawText(histSignal,0,false,0.68, 0.75, 0.9, 0.88);
    tool.DrawTextStar(histSignal,2);
    histBackground->Draw("ESAME");


    leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(histSignal,dataLabel + " (unlike-sign pairs)","p");
    leg1->AddEntry(histBackground,dataLabel + " (like-sign pairs)","p");
    leg1->Draw("same");


    cCanvas->Update();
    cCanvas->Write(variable + "Sum");    
//////////////////////////////////////////////////
// Plot t 
    variable = "phiRp";
    nBins = 100;
    min = -3.0;
    max = 4.5;
    treeBack->Draw(variable +"East>>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",cuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    treeBack->Draw(variable +"West>>" + variable +"Bcg2(" + nBins + "," + min + "," + max + ")",cuts);
    tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Bcg2");
    histBackground->Add(tmpHist);
    tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);


    tree->Draw(variable +"East>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",cuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");
    tree->Draw(variable +"West>>" + variable +"Sig2(" + nBins + "," + min + "," + max + ")",cuts);
    tmpHist = (TH1F*)gPad->GetPrimitive(variable +"Sig2");
    histSignal->Add(tmpHist);   
    histSignal->SetTitle(" ; #phi [rad]; Number of tracks");
    tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
    tool.SetMarkerStyle(histSignal);
    histSignal->Draw("E");
    tool.DrawText(histSignal,0,false,0.68, 0.75, 0.9, 0.88);
    tool.DrawTextStar(histSignal,2);
    histBackground->Draw("ESAME");

                                
    leg1 = new TLegend(0.15, 0.88, 0.28, 0.97);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(histSignal,dataLabel + " (unlike-sign pairs)","p");
    leg1->AddEntry(histBackground,dataLabel + " (like-sign pairs)","p");
    leg1->Draw("same");


    cCanvas->Update();
    cCanvas->Write(variable);
/////////////////////////////////////////////
// Plot phi 
 /*   variable = "deltaPhi";
    nBins = 100;
    min = 0.0;
    max = 200.0;
    treeBack->Draw("TMath::Abs( " + variable +"*57.2957795)>>" + variable +"Bcg1(" + nBins + "," + min + "," + max + ")",cuts);
    histBackground = (TH1F*)gPad->GetPrimitive(variable +"Bcg1");
    tool.SetMarkerStyle(histBackground,2,20,1,2,1,1);

    tree->Draw("TMath::Abs( " + variable +"*57.2957795)>>" + variable +"Sig1(" + nBins + "," + min + "," + max + ")",cuts);
    histSignal = (TH1F*)gPad->GetPrimitive(variable +"Sig1");   
    histSignal->SetTitle(" ; |#Delta #varphi| [deg]; Number of tracks");
    tool.SetGraphStyle(histSignal,4,20,1,4,1,1,0.9,1.4);
    tool.SetMarkerStyle(histSignal);
    histSignal->Draw("E");
    tool.DrawText(histSignal,0,false,0.68, 0.75, 0.9, 0.88);
    tool.DrawTextStar(histSignal,2);
    histBackground->Draw("ESAME");


    leg1 = new TLegend(0.6, 0.65, 0.78, 0.74);
    tool.SetLegendStyle(leg1);
    leg1->AddEntry(histSignal,dataLabel + " (unlike-sign pairs)","p");
    leg1->AddEntry(histBackground,dataLabel + " (like-sign pairs)","p");
    leg1->Draw("same");


    cCanvas->Update();
    cCanvas->Write("Delta" + variable); */
//////////////////////////////////////////////////
    TCanvas *cCanvas2D = new TCanvas("cCanvas2D","cCanvas2D",800,700);
    gPad->SetMargin(0.09,0.13,0.1,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gStyle->SetPalette(1);
    gPad->SetTickx();
    gPad->SetTicky(); 

// Plot XYEastCor vertex
    gPad->SetLogz();
    variable = "yCorrelationsRpEast";
    variable2 = "xCorrelationsRpEast";
    cutsOption = cuts;
    nBins = 100;
    min = -1;
    max = 1.4;
    nBins2 = 100;
    min2 = -1;
    max2 = 1.4;

    tree->Draw(variable + ":" + variable2 + ">>Sig(" + nBins + "," + min + "," + max + "," + nBins2 + "," + min2 + "," + max2 +")",cuts,"colz");
    hist2DSignal = (TH2F*)gPad->GetPrimitive("Sig");
    hist2DSignal->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");
    tool.SetGraphStyle(hist2DSignal,4,20,1,4,1,1,0.9,0.7);
    hist2DSignal->Draw("colz");
    tool.DrawText(hist2DSignal,0,true,0.61,0.75,0.76,0.9);
    tool.DrawTextStar(hist2DSignal, 1);
    /*if(TEXT){
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
    }*/
    cCanvas2D->Update();
    cCanvas2D->Write(variable + "VS" + variable2);
//////////////////////////////////////////
//////////////////////////////////////////////////
// Plot XYEastCor vertex
    gPad->SetLogz();
    variable = "yCorrelationsRpWest";
    variable2 = "xCorrelationsRpWest";
    cutsOption = cuts;
    nBins = 100;
    min = -1;
    max = 1.4;
    nBins2 = 100;
    min2 = -1;
    max2 = 1.4;

    tree->Draw(variable + ":" + variable2 + ">>Sig(" + nBins + "," + min + "," + max + "," + nBins2 + "," + min2 + "," + max2 +")",cuts,"colz");
    hist2DSignal = (TH2F*)gPad->GetPrimitive("Sig");
    hist2DSignal->SetTitle(" ; p_{x} [GeV/c]; p_{y} [GeV/c]");
    tool.SetGraphStyle(hist2DSignal,4,20,1,4,1,1,0.9,0.7);
    hist2DSignal->Draw("colz");
    tool.DrawText(hist2DSignal,0,true,0.61,0.75,0.76,0.9);
    tool.DrawTextStar(hist2DSignal, 1);
    /*if(TEXT){
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
    }*/
    cCanvas2D->Update();
    cCanvas2D->Write(variable + "VS" + variable2);
//////////////////////////////////////////


}//RPplots::PlotHistogram


