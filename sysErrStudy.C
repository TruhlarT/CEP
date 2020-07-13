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

const double speedOfLight = 299792458; // m/s
const double pionMass = 0.13957; // GeV /c^2
const double protonMass = 0.93827; // GeV /c^2
const double convertToDegree = 57.2957795;

enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};

TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};
TString sideLabel[nSides] = { TString("East"), TString("West")};
TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};

TFile* fout;
TFile* TPCeff;
TFile* TOFeff;
TFile* MyEff;
TFile* data;

TTree* recTree;

TH1F* hInvMass[4];

TH3F* hTPCeff;
TH3F* hTOFeff;
TH3F* hMyEff3D;
TH1* hMyEffPt;

Int_t nTracks, totalCharge, nTofTrks; 
UInt_t runNumber;
Double_t VPDSumEast, VPDSumWest, VPDTimeDiff;


Double_t chiPair[nParticles]; 
Double_t invMass[nParticles];
Double_t missingPt, deltaTOF, mSquared;
Double_t deltaDeltaTOF[nParticles];
Double_t deltaTOFExpected[nParticles];

/////////////////////////////////

Bool_t elastic, fourPiState;
Bool_t trigger[17];

UInt_t BBCSmall[nSides]; // BBC truncated sum, small tiles
UInt_t BBCLarge[nSides]; // BBC truncated sum, large tiles 
Double_t xCorrelationsRp[nSides];
Double_t yCorrelationsRp[nSides];
Double_t thetaRp[nSides];
Double_t phiRp[nSides];
Double_t timeRp[nSides];
Double_t pRp[nSides];
Double_t ptRp[nSides];
Double_t etaRp[nSides];
Double_t rpX[nSides];
Double_t rpZ[nSides];
Double_t rpY[nSides];
Double_t t[nSides];
Double_t xi[nSides];

Double_t ADC[nRomanPots][2]; // RP_ID   0,    1,    2,   3,   4,   5,   6, 7
Double_t TAC[nRomanPots][2]; // RP_name E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D

Double_t dEdx[4];
Double_t momentum[4];
Double_t transMomentum[4];
Double_t TOFtime[4];
Double_t TOFlength[4];
Double_t charge[4];
Double_t nSigmaTPC[nParticles][4];
Double_t vertexId[4];
Double_t vertexesZ[4];
Double_t DcaXY[4];
Double_t DcaZ[4];
Double_t NhitsFit[4];
Double_t NhitsDEdx[4];
Double_t Eta[4];
Double_t Phi[4];
Double_t Chi2[4];


const double nFitLim[2] = { 26.0 , 24.0}; 
const double nDEdxLim[2] = { 16.0 , 14.0}; 
const double dcaLim[2] = { 1.3, 1.7}; 

void Init();
void ConnectInput();
void Make();

void sysErrStudy()
{
    TString rafalEff = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/etaPhiEfficiency_16_01_19_delta015_twoRuns.root";
    TString myEff = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effMy.root";
    TString output = "/home/truhlar/Desktop/STAR/CEP/Analysis/Outputs/sysStudy.root";
    TString TOFeffInput = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/effWithBinningForSystematics.root";
    TString input = "/home/truhlar/Desktop/STAR/CEP/Analysis/Data/ppRun17.root";

    TPCeff = TFile::Open(rafalEff, "read");
    if (!TPCeff)
    {
        cout<<"Error: cannot open "<<rafalEff<<endl;
        return;
    }

    TOFeff = TFile::Open(TOFeffInput, "read");
    if (!TOFeff)
    {
        cout<<"Error: cannot open "<<TOFeffInput<<endl;
        return;
    }

    MyEff = TFile::Open(myEff, "read");
    if (!MyEff)
    {
        cout<<"Error: cannot open "<<myEff<<endl;
        return;
    }

    hTPCeff = (TH3F*)TPCeff -> Get(Form("hTPCEffiCD%i120",3));
    hTOFeff = (TH3F*)TOFeff -> Get(Form("hTOFEffiCD%i12",3));
//    hMyEff3D = (TH3F*)MyEff -> Get("effRafal");
    hMyEff3D = (TH3F*)MyEff -> Get("effMy_0");
    hMyEffPt = hMyEff3D -> Project3D("y");
    double nBinsXZ = 20*24;
    hMyEffPt->Scale(1/nBinsXZ);

    data = TFile::Open(input, "read");
    if (!data)
    {
        cout<<"Error: cannot open "<<input<<endl;
        return;
    }


    fout = new TFile(output,"RECREATE");
    hMyEffPt->Write("effPt");
    Init(); // Preparing histograms 
    ConnectInput(); // Connecting input


    Long64_t nev = recTree->GetEntries();
    cout<<"Proccesing "<<nev<<" events"<<endl;

    for(Long64_t iev=0; iev<nev; ++iev) 
    { //get the event
        recTree->GetEntry(iev); 
        Make();
    } 

    TCanvas* newCanvas = new TCanvas("newCanvas","newCanvas",800,700);
    gPad->SetMargin(0.11,0.02,0.105,0.02); // (Float_t left, Float_t right, Float_t bottom, Float_t top)
    gPad->SetTickx();
    gPad->SetTicky(); 
    gPad->SetLogy(0);
    gStyle->SetOptStat("11111");
    gStyle->SetOptStat("");
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);      //axis line
    gStyle->SetFrameLineWidth(2); //frame line
    TString stateLabel[] = {"#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}"}; 

    gPad->SetMargin(0.9,0.02,0.3,0.02);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"Y");
    gStyle->SetTickLength(0.03,"X");
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.295, 1, 1.0);
    pad1->SetTopMargin(0.04);
    pad1->SetRightMargin(0.02);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->Draw();             // Draw the upper pad: pad1
             
    newCanvas->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.2985);
    pad2->SetTopMargin(0);
    pad2->SetRightMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
/*
    Double_t scaleFactor; 
    for (int i = 0; i < 4; ++i)
    {
        scaleFactor =   1 /hInvMass[i]->Integral();
        hInvMass[i]->Scale(scaleFactor);
    }
*/

    TH1F *h4 = (TH1F*)hInvMass[0]->Clone("h4");
    pad1->cd();
    hInvMass[1]->SetMarkerColor(2);
    hInvMass[1]->SetMarkerSize(1);
    hInvMass[1]->SetMarkerStyle(24);
    hInvMass[1]->SetStats(false);

    hInvMass[2]->SetMarkerColor(1);
    hInvMass[2]->SetMarkerSize(1);
    hInvMass[2]->SetMarkerStyle(25);
    hInvMass[2]->SetStats(false);

    hInvMass[0]->SetStats(false);
    hInvMass[0]->SetTitle(" ; m(#pi^{+}#pi^{-}) [GeV/c^{2}]; counts");
    hInvMass[0]->GetXaxis()->SetTitleFont(42);
    hInvMass[0]->GetYaxis()->SetTitleFont(42);
    hInvMass[0]->GetXaxis()->SetLabelFont(42);
    hInvMass[0]->GetYaxis()->SetLabelFont(42);
    hInvMass[0]->GetXaxis()->SetTitleSize(0.045);
    hInvMass[0]->GetYaxis()->SetTitleSize(0.045);
    hInvMass[0]->GetXaxis()->SetTitleOffset(0.9);
    hInvMass[0]->GetYaxis()->SetTitleOffset(1.3);
    hInvMass[0]->GetYaxis()->SetRangeUser(10.0,hInvMass[2]->GetMaximum()*1.1);
    hInvMass[0]->SetMarkerColor(4);
    hInvMass[0]->SetMarkerSize(1);
    hInvMass[0]->SetMarkerStyle(26);
    hInvMass[0]->SetLineColor(4);
    hInvMass[0]->SetLineStyle(1);
    hInvMass[0]->SetLineWidth(1);
    hInvMass[0]->Draw("E");
    hInvMass[1]->Draw("ESAME");
    hInvMass[2]->Draw("ESAME");


    TPaveText *textPub = new TPaveText(0.7,0.75,0.92,0.88,"brNDC");
    textPub -> SetTextSize(0.04);
    textPub -> SetTextAlign(22);
    textPub -> SetFillColor(0);
    textPub -> SetTextFont(42);
    textPub -> AddText("p + p #rightarrow p + #pi^{+}#pi^{-} + p");
    textPub -> AddText("#sqrt{s} = 510 GeV");
    textPub -> Draw("same");

    TPaveText *textSTAR;
    textSTAR = new TPaveText(0.75,0.89,0.9,0.95,"brNDC");
    textSTAR -> SetTextSize(0.04);
    textSTAR -> SetFillColor(0);
    textSTAR -> SetTextFont(62);
    textSTAR->AddText("STAR Internal");
    textSTAR -> Draw("same");

    TLegend* leg1 = new TLegend(0.6, 0.50, 0.78, 0.74);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetTextFont(42);
    leg1->AddEntry(hInvMass[0],"Uncorr. standard cuts","p");
    leg1->AddEntry(hInvMass[1],"Uncorr. strict cuts","p");
    leg1->AddEntry(hInvMass[2],"Uncorr. loose cuts","p");
    leg1->Draw("same");

    pad2->cd();
    // Define the ratio plot
    TH1F *h3 = (TH1F*)hInvMass[0]->Clone("h3");
    //h3->GetYaxis()->SetLabelSize(0.);
    h3->SetMarkerColor(2);
    h3->SetLineColor(2);
    h3->SetMinimum(0.7);  // Define Y ..
    h3->SetMaximum(1.3); // .. range
    h3->Divide(hInvMass[1]);
    h3->SetMarkerStyle(20);
    h3->Draw("HIST");       // Draw the ratio plot

    
    //h4->GetYaxis()->SetLabelSize(0.);
    h4->SetMarkerColor(6);
    h4->SetLineColor(6);
    h4->Divide(hInvMass[2]);
    h4->SetMarkerStyle(20);
    h4->Draw("HISTSAME");       // Draw the ratio plot

    // Ratio plot (h3) settings
    h3->SetTitle(""); // Remove the ratio title

    // Y axis ratio plot settings
    h3->GetYaxis()->SetTitle("Ratio");
    h3->GetYaxis()->SetNdivisions(505);
    h3->GetYaxis()->SetTitleSize(25);
    h3->GetYaxis()->SetTitleFont(43);
    h3->GetYaxis()->SetTitleOffset(1.5);
    h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h3->GetYaxis()->SetLabelSize(0.0);

    // X axis ratio plot settings
    h3->GetXaxis()->SetTitleSize(25);
    h3->GetXaxis()->SetTitleFont(43);
    h3->GetXaxis()->SetTitleOffset(3.2);
    h3->GetXaxis()->SetLabelOffset(0.025);
    h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h3->GetXaxis()->SetLabelSize(20);

    TGaxis *axis = new TGaxis(h3->GetXaxis()->GetXmin(), 0.7, h3->GetXaxis()->GetXmin(), 1.3, 0.7, 1.3, 4, "");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(20);
    axis->Draw("same");

    TLine *unity = new TLine(h3->GetXaxis()->GetXmin(), 1., h3->GetXaxis()->GetXmax() , 1.);
    unity->SetLineColor(kBlack);
    unity->SetLineWidth(2);
    unity->Draw();


    gStyle->SetOptStat("");
    newCanvas->Update();
    newCanvas->Write("massComaprisson");
    newCanvas->Close();
    fout->Write();
    fout->Close();
    
}




void Make()
{
    double effTotal, effTPC, effTOF;
    double effMy, effVar;
    double effPt, effPtVar;
   // cout<< vertexesZ[0] <<" "<< NhitsFit[0]<<" "<<NhitsFit[1] <<" "<< NhitsDEdx[0]<<" "<<NhitsDEdx[1] <<" "<<DcaZ[0] <<" "<<DcaZ[1] <<" "<<DcaXY[0] <<" "<<DcaXY[1] <<" "<<Eta[0] <<" "<<Eta[1] <<" "<< !fourPiState<<endl; 
    if(vertexesZ[0] < 80 && vertexesZ[0] > -80 && DcaZ[0] < 1 && DcaZ[0] > -1 && DcaZ[1] < 1 && DcaZ[1] > -1 && Eta[0] > -0.7 && Eta[0] < 0.7 && Eta[1] > -0.7 && Eta[1] < 0.7 &&  t[0] < -0.12 && t[1] < -0.12 && t[0] > -1.0  && t[1] > -1.0 && !fourPiState)
    {
        effTotal = 1;
        effMy = 1;
        effPt = 1;
        if(chiPair[Pion] > 9 && chiPair[Kaon] > 9 && chiPair[Proton] < 9 && mSquared > 0.6) // it is... proton!
        {


        }
        else if(chiPair[Pion] > 9 && chiPair[Kaon] < 9 && chiPair[Proton] > 9 && mSquared > 0.15) // it is... kaon!
        {

        }
        else if( chiPair[Pion] < 12) // it is... pion!
        {
            if(transMomentum[0] > 0.2 && transMomentum[1] > 0.2)
            {
                if( NhitsFit[0] >= 25 && NhitsFit[1] >= 25 && NhitsDEdx[0] >= 15 && NhitsDEdx[1] >= 15 
                    && DcaXY[0] < 1.5 && DcaXY[1] < 1.5 )
                    hInvMass[0]->Fill(invMass[Pion]);
                for (int ind = 0; ind < 2; ++ind)
                {
                    if( NhitsFit[0] >= nFitLim[ind] && NhitsFit[1] >= nFitLim[ind] && NhitsDEdx[0] >= nDEdxLim[ind] && NhitsDEdx[1] >= nDEdxLim[ind] 
                        && DcaXY[0] < dcaLim[ind] && DcaXY[1] < dcaLim[ind] )
                        hInvMass[ind + 1]->Fill(invMass[Pion]);

                }

            }

        }

    }

}


void Init(){

    hInvMass[0] = new TH1F("invMass", "uncorrected standard inv. mass", 64, 0.3, 3.5);
    hInvMass[1] = new TH1F("invMassLow", "uncorrected inv. mass with strict cuts", 64, 0.3, 3.5);
    hInvMass[2] = new TH1F("invMassUp", "uncorrected inv. mass with loose cuts ", 64, 0.3, 3.5);
}

void ConnectInput(){
    recTree = dynamic_cast<TTree*>( data->Get("recTree") );
    if (!recTree)
    {
        cout<<"Error: cannot get recTree"<<endl;
        return;
    }

// PID and some quality event info
    recTree->SetBranchAddress("missingPt", &missingPt);
    recTree->SetBranchAddress("deltaTOF", &deltaTOF);
    recTree->SetBranchAddress("mSquared", &mSquared); 
    recTree->SetBranchAddress("nSigTrk1Pion", &nSigmaTPC[Pion][0]);
    recTree->SetBranchAddress("nSigTrk2Pion", &nSigmaTPC[Pion][1]);
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        recTree->SetBranchAddress("invMass" + particleLables[iPart], &invMass[iPart]);
        recTree->SetBranchAddress("chiPair" + particleLables[iPart], &chiPair[iPart]);
        recTree->SetBranchAddress("deltaTOFExpected" + particleLables[iPart], &deltaTOFExpected[iPart]);
        recTree->SetBranchAddress("deltaDeltaTOF" + particleLables[iPart], &deltaDeltaTOF[iPart]);  
    }


// Vertex info
    recTree->SetBranchAddress("vertexZ", &vertexesZ[0]);

// Central track info
    for (int i = 0; i < 4; ++i)
    {
        recTree->SetBranchAddress(Form("dEdx%i",i), &dEdx[i]);
        recTree->SetBranchAddress(Form("momentum%i",i), &momentum[i]);
        recTree->SetBranchAddress(Form("transMomentum%i",i), &transMomentum[i]);
        recTree->SetBranchAddress(Form("charge%i",i), &charge[i]);
        recTree->SetBranchAddress(Form("TOFtime%i",i), &TOFtime[i]);
        recTree->SetBranchAddress(Form("TOFlength%i",i), &TOFlength[i]);
        recTree->SetBranchAddress(Form("DcaXY%i",i), &DcaXY[i]);
        recTree->SetBranchAddress(Form("DcaZ%i",i), &DcaZ[i]);
        recTree->SetBranchAddress(Form("NhitsFit%i",i), &NhitsFit[i]);
        recTree->SetBranchAddress(Form("NhitsDEdx%i",i), &NhitsDEdx[i]);
        recTree->SetBranchAddress(Form("Eta%i",i), &Eta[i]);
        recTree->SetBranchAddress(Form("Phi%i",i), &Phi[i]);
        recTree->SetBranchAddress(Form("Chi2%i",i), &Chi2[i]);
    }
// RP track info  
    for (int i = 0; i < nSides; ++i)
    {
        recTree->SetBranchAddress("rpX" + sideLabel[i], &rpX[i]);
        recTree->SetBranchAddress("rpY" + sideLabel[i], &rpY[i]);
        recTree->SetBranchAddress("rpZ" + sideLabel[i], &rpZ[i]);
        recTree->SetBranchAddress("thetaRp" + sideLabel[i], &thetaRp[i]);
        recTree->SetBranchAddress("phiRp" + sideLabel[i], &phiRp[i]);
        recTree->SetBranchAddress("timeRp" + sideLabel[i], &timeRp[i]);
        recTree->SetBranchAddress("pRp" + sideLabel[i], &pRp[i]);
        recTree->SetBranchAddress("ptRp" + sideLabel[i], &ptRp[i]);
        recTree->SetBranchAddress("etaRp" + sideLabel[i], &etaRp[i]);
        recTree->SetBranchAddress("xCorrelationsRp" + sideLabel[i], &xCorrelationsRp[i]);
        recTree->SetBranchAddress("yCorrelationsRp" + sideLabel[i], &yCorrelationsRp[i]);
        recTree->SetBranchAddress("t" + sideLabel[i], &t[i]);
        recTree->SetBranchAddress("xi" + sideLabel[i], &xi[i]);
    }
// RP event info
    for (int i = 0; i < nRomanPots; ++i)
    {
        recTree->SetBranchAddress("ADC_" + rpNames[i] + "V", &ADC[i][0]);
        recTree->SetBranchAddress("ADC_" + rpNames[i] + "H", &ADC[i][1]);
        recTree->SetBranchAddress("TAC_" + rpNames[i] + "V", &TAC[i][0]);
        recTree->SetBranchAddress("TAC_" + rpNames[i] + "H", &TAC[i][1]);
    }

// event info
    recTree->SetBranchAddress("elastic", &elastic);
    recTree->SetBranchAddress("fourPiState", &fourPiState);
    recTree->SetBranchAddress("runNumber", &runNumber);
    recTree->SetBranchAddress("VPDTimeDiff", &VPDTimeDiff);
    recTree->SetBranchAddress("VPDSumWest", &VPDSumWest);
    recTree->SetBranchAddress("VPDSumEast", &VPDSumEast);
    for (int i = 0; i < nSides; ++i)
    {
        recTree->SetBranchAddress("BBCSmall" + sideLabel[i], &BBCSmall[i]);
        recTree->SetBranchAddress("BBCLarge" + sideLabel[i], &BBCLarge[i]); 
    }
    recTree->SetBranchAddress("RP_CPT2_570701", &trigger[3]);
    recTree->SetBranchAddress("RP_CPT2noBBCL_570705", &trigger[7]);
    recTree->SetBranchAddress("RP_CPT2_570711", &trigger[9]);
    recTree->SetBranchAddress("RP_CPT2_590701", &trigger[12]);
    recTree->SetBranchAddress("RP_CPT2noBBCL_590705", &trigger[14]);
    recTree->SetBranchAddress("RP_CPTnoBBCL_590708", &trigger[15]);
}