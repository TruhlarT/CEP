// Run by: ./Analysis input index
// e.g. ./Analysis /gpfs01/star/pwg/truhlar/Final/CPtrig/merge_files/StUPCRP_production.list -1
// or you can open just root file
// ./Analysis /star/data01/pwg_tasks/upc02/Part9/18143045/18143045.root


// c++ headers
#include <iostream>
#include <string>     // std::string, std::to_string
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

// picoDst headers
#include "StRPEvent.h"
#include "StUPCRpsTrack.h"
#include "StUPCRpsTrackPoint.h"

#include "StUPCEvent.h"
#include "StUPCTrack.h"
#include "StUPCBemcCluster.h"
#include "StUPCVertex.h"

using namespace std;

enum {kAll = 1,kCPtrig, kEl, kInEl, kTOF2t, kSameVrtx, kTotCH0, kMissPt, kMaxCount};
enum SIDE {E = 0, East = 0, W = 1, West = 1, nSides};
enum PARTICLES {Pion = 0, Kaon = 1, Proton = 2, nParticles};
enum XY_COORDINATE {X, Y, nCoordinates};
enum RP_ID {E1U, E1D, E2U, E2D, W1U, W1D, W2U, W2D, nRomanPots};
enum BRANCH_ID { EU, ED, WU, WD, nBranches };
enum ARM_ID { EU_WU, ED_WD, nArms };
enum STATION_ID { E1, E2, W1, W2, nStations };
enum RP_CONFIGURATION {EUD = 0, EDU = 1, IUU = 2, IDD = 3, nConfiguration}; 

const double particleMass[nParticles] = { 0.13957, 0.49367, 0.93827}; // GeV /c^2
const double speedOfLight = 299792458; // m/s
const double beamMomentum = 254.9; // GeV
const int nTriggers = 17;
const int triggerID[] = { 570209, 570219, 570229, 570701, 570702, 570703, 570704, 570705, 
                  570709, 570711, 570712, 570719, 590701, 590703, 590705, 590708, 590709};
// 570702 RP_UPC // 570712 RP_UPC // 570703 RP_SDT // 570709 RP_ET // 570719 RP_ET // 570701 RP_CPT2 // 570711 RP_CPT2 // 570705 RP_CPT2noBBCL // 570704 RP_Zerobias // 590703 RP_SDT // 590709 RP_ET // 590701 RP_CPT2 // 590705 RP_CPT2noBBCL // 590708 RP_CPTnoBBCL // 570209 JPsi*HTTP // 570219 JPsi*HTTP // 570229 JPsi*HTTP

Bool_t trigger[nTriggers];


TString branchName[nBranches] = { TString("EU"), TString("ED"), TString("WU"), TString("WD") };
TString summaryLabels[10] = { TString("All"), TString("CP Trigger"), TString("Elastic"), TString("Inelastic"), TString("2 TPC-TOF tracks"), 
                              TString("Same vertex"), TString("TotCharge 0"), TString("MissingPt < 0.1 GeV"), 
                              TString(""), TString("")};
TString summaryLabels2[10] = { TString("All"), TString("CP Trigger"), TString("Elastic"), TString("Inelastic"), TString("4 TPC-TOF tracks"), 
                              TString("Same vertex"), TString("TotCharge 0"), TString("MissingPt < 0.1 GeV"), 
                              TString(""), TString("")};
TString systemID[3] = { TString("Combi"), TString("InEl"), TString("El")};
TString systemState[4] = { TString("TPC2t"), TString("TOF2trk"), TString("Q0"), TString("Excl")};
TString configLabels[4] = { TString("EUD"), TString("EDU"), TString("IUU"), TString("IDD")};
TString systemLabels[2] = { TString("2part"), TString("4part")};
TString particleLables[nParticles] = { TString("Pion"), TString("Kaon"), TString("Proton")};
TString sideLabel[nSides] = { TString("East"), TString("West")};
TString rpNames[nRomanPots] = { TString("E1U"), TString("E1D"), TString("E2U"), TString("E2D"), TString("W1U"), TString("W1D"), TString("W2U"), TString("W2D")};

TH1I* hAnalysisFlow; // control plot
TH1I* hAnalysisFlow2; // control plot for 4 pion state
TH1F* hTriggerBits; // control plot
TH1F* hConfiguration;
TH1F* hNumberTrackPerBranch[nBranches]; // number of track per each brunch (west up, west down, east up, east down)
TH1F* hNumberTrack;
TH1F* hTrackPlanesUsed;
TH1F* hTrackPointPlanesUsed;
TH1F* hNumberTOFmatchTrack;
TH1F* hInvalidRPtrack; 

TH1F* hEtaTest[2];
TH2D* hEtaPhi[3];

TH1F* hMissingPt[12][2]; 
TH2D* hdEdxVsqP[12][2];

TFile *infile;
TFile *outfile;
TChain *upcChain;
TChain *rpChain;
StUPCEvent *upcEvt;
StRPEvent *rpEvt;

TTree *upcTree;
TTree *rpTree;

TTree *recTree, *bcgTree;
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


void Init();
void Make();
void FillPlots(int state, int configuration, int isFourPiState );
TFile *CreateOutputTree(const string& out);
bool ConnectInput(const string& in, int fileId);
void Clear();
void CalculatePID();

//_____________________________________________________________________________
int main(int argc, char** argv) {
    //open output file
    outfile = CreateOutputTree("AnalysisOutput.root"); 
    if(!outfile) 
    {
        cout << "Can not open output file." << endl; 
        return -1;
    }

    Init();
    int fileId = -1;
    const string& input = argv[1];
    if( input.find(".root") == string::npos )
        fileId = atoi(argv[2]);

    if(!ConnectInput(argv[1], fileId))
    {
        cout << "No input." << endl; 
        return 1;
    }
    //ask for number of events
    Long64_t nev = upcTree->GetEntries();
    cout<<"Proccesing "<<nev<<" events"<<endl;
    //return 0;
    //event loop
    //nev = 1000;
    for(Long64_t iev=0; iev<nev; ++iev) 
    { //get the event
        upcTree->GetEntry(iev); 
        Make();
    } 

    //close the outputs
    outfile->Write();
    outfile->Close();
    cout<<"Ending Analysis... GOOD BYE!"<<endl;
    return 0;
}//main

void Init(){
    hEtaTest[0] = new TH1F("EtaTestTrg", "EtaTestTrg", 100, -1.5, 1.5);
    hEtaTest[1] = new TH1F("EtaTestAss", "EtaTestAss", 100, -1.5, 1.5);
    hEtaPhi[0] = new TH2D("EtaPhiAll","EtaPhiAll",200,-2,2,100,-3.5,3.5);
    hEtaPhi[1] = new TH2D("EtaPhiMatched","EtaPhiMatched",200,-2,2,100,-3.5,3.5);
    hEtaPhi[2] = new TH2D("EtaPhiNotMtch","EtaPhiNotMtch",200,-2,2,100,-3.5,3.5);

    hAnalysisFlow = new TH1I("AnalysisFlow", "CutsFlow", kMaxCount-1, 1, kMaxCount);
    hAnalysisFlow2 = new TH1I("AnalysisFlow2", "CutsFlow", kMaxCount-1, 1, kMaxCount); 
    hInvalidRPtrack = new TH1F("invalidRPtrack","invalidRPtrack", 6,0,4);
    //hAnalysisFlow = new TH1F("AnalysisFlow", "AnalysisFlow", 10, -0.5, 9.5);
    for(int tb=1; tb<kMaxCount; ++tb) 
        hAnalysisFlow->GetXaxis()->SetBinLabel(tb, summaryLabels[tb-1]);

    for(int tb=1; tb<kMaxCount; ++tb) 
        hAnalysisFlow2->GetXaxis()->SetBinLabel(tb, summaryLabels2[tb-1]);

    hTriggerBits = new TH1F("TriggerBits", "TriggerBits", nTriggers, -0.5, 16.5);
    for(int tb=0; tb<nTriggers; ++tb){
        TString label; label.Form("%d",triggerID[tb]);
        hTriggerBits->GetXaxis()->SetBinLabel(tb+1, label);
    }

    hConfiguration = new TH1F("Configuration", "Track Configuration", 4, -0.5, 3.5);
    for(int tb=0; tb<4; ++tb) 
        hConfiguration->GetXaxis()->SetBinLabel(tb+1, configLabels[tb]);

    for(int i=0; i<nBranches; ++i)
        hNumberTrackPerBranch[i] = new TH1F("NumberTracksPerBranch,"+branchName[i],"Number of tracks in branch "+branchName[i], 8, -0.5, 7.5);

    hNumberTrack = new TH1F("NumberTracks", "Number of Tracks in RP", 40, -0.5, 39.5);
    hNumberTOFmatchTrack = new TH1F("hNumberTOFmatchTrack", "Number of TOF matched tracks", 40, -0.5, 39.5);
    hTrackPlanesUsed = new TH1F("hTrackPlanesUsed", "Number of planes used for track reconstruction", 40, -0.5, 39.5);
    hTrackPointPlanesUsed = new TH1F("hTrackPointPlanesUsed", "Number of planes used for track point reconstruction", 40, -0.5, 39.5); 

    for(int i=0; i<12;++i){
        if(i==0){
            outfile->mkdir("All")->cd();
        }else if(i==4){
            outfile->cd();
            outfile->mkdir("Inelastic")->cd();
        }else if(i==8){
            outfile->cd();
            outfile->mkdir("Elastic")->cd();
        }

        for( int state = 0; state < 2; ++state)
        {
            hdEdxVsqP[i][state] = new TH2D("dEdxVsqP_"+systemState[i%4]+"_"+systemID[i/4] + systemLabels[state],"dE/dx Vs #frac{q}{e} P",200,-2,2,100,0,20);
            hMissingPt[i][state] = new TH1F( "MissingPt_"+systemState[i%4]+"_"+systemID[i/4] + systemLabels[state], "p^{miss}_{T} [GeV/c]", 200, 0, 2 );
        }
    }


    outfile->cd();
}

void Make(){
    hAnalysisFlow->Fill(kAll);
    hAnalysisFlow2->Fill(kAll);

    Clear();
    bool CPTtrigger = false;
    for(int var = 0; var < nTriggers; ++var)
    {
        if(upcEvt->isTrigger(triggerID[var]))
        {
            hTriggerBits->Fill(var);
            trigger[var] = true;
            //if(var==3 || var==7 || var==9 || var==12 || var==14 || var==15)
                // RP_CPT2, RP_CPT2noBBCL, RP_CPT2, RP_CPT2, RP_CPT2noBBCL, RP_CPTnoBBCL
            if(var==7 || var==14 )
                // RP_CPT2noBBCL, RP_CPT2noBBCL
                CPTtrigger=true;
        }
    }

    if(!CPTtrigger)
        return;

    hAnalysisFlow->Fill(kCPtrig);
    hAnalysisFlow2->Fill(kCPtrig);
    nTracks = (Int_t) rpEvt->getNumberOfTracks();

    // Vector below will be filled with indices of good-quality tracks
    vector<int> rpTrackIdVec_perBranch[nBranches];
    vector<int> rpTrackIdVec_perSide[nSides];
    int numberOfTracks = 0;
    int numberOfTracksPerBranch[nBranches] = {0, 0, 0, 0};

    // Loop over all tracks reconstructed in Roman Pots  
    for(unsigned int k = 0; k < rpEvt->getNumberOfTracks(); ++k)
    {
        // Get pointer to k-th track in Roman Pot data collection
        StUPCRpsTrack *trk = rpEvt->getTrack(k);
        trk->setEvent(rpEvt);
        // Get ID of a branch in which this k-th track was reconstructed
        int j = trk->branch();
        int side = j<2 ? E : W;

        ++numberOfTracks;
        ++numberOfTracksPerBranch[j];
        hTrackPlanesUsed->Fill(trk->planesUsed());
        if(trk->getTrackPoint(0))
            hTrackPointPlanesUsed->Fill(trk->getTrackPoint(0)->planesUsed());

        if(trk->getTrackPoint(1))
            hTrackPointPlanesUsed->Fill(trk->getTrackPoint(1)->planesUsed());

        // If track is global (made of track-points in 1st and 2nd station)
        // and all 8 planes were used in reconstruction - store its ID
        if( trk->type()==StUPCRpsTrack::rpsGlobal && trk->planesUsed()==8)
        {
            rpTrackIdVec_perBranch[j].push_back( k );
            rpTrackIdVec_perSide[side].push_back( k );
        } 
        // a bit looser selection
        /* if( (trk->getTrackPoint(0) ? trk->getTrackPoint(0)->planesUsed()>=3 : false) &&
            (trk->getTrackPoint(1) ? trk->getTrackPoint(1)->planesUsed()>=3 : false) && 
            trk->type()==StUPCRpsTrack::rpsGlobal && trk->planesUsed()<=8)
        {
            rpTrackIdVec_perBranch[j].push_back( k );
            rpTrackIdVec_perSide[side].push_back( k );
        } */
    }

    hNumberTrack->Fill(numberOfTracks);
    for(int i=0; i<nBranches; ++i) 
        hNumberTrackPerBranch[i]->Fill(numberOfTracksPerBranch[i]);

// Loop over arms - check if have good-quality tracks, selecting branch combination
    for(int i=0; i<nConfiguration; ++i)
    { 
    // Define IDs of branches based on known ID of the arm
        int branch[nSides];
        switch(i) {
            case EUD :
            {
                branch[E] = EU;
                branch[W] = WD;
            }
            break;
            case EDU :
            {
                branch[E] = ED;
                branch[W] = WU;
            }
            break; 
            case IUU :
            {
                branch[E] = EU;
                branch[W] = WU;
            }
            break; 
            case IDD :
            {
                branch[E] = ED;
                branch[W] = WD;
            }
            break;        
        }

        //If exactly one good-quality track in each branch in the arm and there is no track in the other RP do some staff
        if( rpTrackIdVec_perBranch[ branch[E] ].size()==1 
        && rpTrackIdVec_perBranch[ branch[W] ].size()==1
        && rpTrackIdVec_perSide[E].size()==1 && rpTrackIdVec_perSide[W].size()==1)
        {
            // Get pointers to good-quality tracks
            if(i==EUD || i==EDU)
            {
                hAnalysisFlow->Fill(kEl);
                hAnalysisFlow2->Fill(kEl);
                elastic = true;
            }else
            {
                hAnalysisFlow->Fill(kInEl);
                hAnalysisFlow2->Fill(kInEl);
                elastic = false;
            }
            hConfiguration->Fill(i);
            TLorentzVector centralTracksTotalFourMomentum[nParticles];

            // loop over all TPC tracks
            for(int j=0; j<upcEvt->getNumberOfTracks(); ++j)
            {
            // get TPC track object
                const StUPCTrack* trk = upcEvt->getTrack(j);

                TLorentzVector trkLVector[nParticles];
                for (int iPart = 0; iPart < nParticles; ++iPart)
                    trk->getLorentzVector(trkLVector[iPart], particleMass[iPart]);

                if( !trk->getFlag(StUPCTrack::kPrimary)) 
                    continue;

                hEtaPhi[0]->Fill(trk->getEta(),trk->getPhi());

                if(!trk->getFlag(StUPCTrack::kTof) )
                {
                    hEtaPhi[2]->Fill(trk->getEta(),trk->getPhi());
                    continue;
                } 

                hEtaPhi[1]->Fill(trk->getEta(),trk->getPhi());
                if(nTofTrks>3)
                {
                    nTofTrks++;
                    continue;
                }

                // read basic information about the track     
                dEdx[nTofTrks] = trk->getDEdxSignal()*1000000;
                momentum[nTofTrks] = trkLVector[0].P();
                transMomentum[nTofTrks] = trk->getPt();
                charge[nTofTrks] = trk->getCharge();
                TOFtime[nTofTrks] = trk->getTofTime();
                TOFlength[nTofTrks] = trk->getTofPathLength();
                DcaXY[nTofTrks] = trk->getDcaXY();
                DcaZ[nTofTrks] = trk->getDcaZ();
                NhitsFit[nTofTrks] = trk->getNhitsFit();
                NhitsDEdx[nTofTrks] = trk->getNhitsDEdx();
                Eta[nTofTrks] = trk->getEta();
                Phi[nTofTrks] = trk->getPhi();
                Chi2[nTofTrks] = trk->getChi2();
                vertexId[nTofTrks] = trk->getVertexId();
                vertexesZ[nTofTrks] = trk->getVertex()->getPosZ();
                totalCharge += static_cast<int>( trk->getCharge() );
                nSigmaTPC[Pion][nTofTrks] = trk->getNSigmasTPCPion();
                nSigmaTPC[Kaon][nTofTrks] = trk->getNSigmasTPCKaon();
                nSigmaTPC[Proton][nTofTrks] = trk->getNSigmasTPCProton();

                for (int iPart = 0; iPart < nParticles; ++iPart)
                    centralTracksTotalFourMomentum[iPart] += trkLVector[iPart];
                
                ++nTofTrks;
            }

            hNumberTOFmatchTrack->Fill(nTofTrks);
            missingPt = (centralTracksTotalFourMomentum[0].Vect() + rpEvt->getTrack( rpTrackIdVec_perSide[E][0] )->pVec() + rpEvt->getTrack( rpTrackIdVec_perSide[W][0] )->pVec() ).Pt();
            
            for (int iPart = 0; iPart < nParticles; ++iPart)
                invMass[iPart] = centralTracksTotalFourMomentum[iPart].M();
      
            FillPlots(0,i,0);
            FillPlots(0,i,1);
            fourPiState = false;
            if( nTofTrks==4)
            {
                fourPiState = true;
                hAnalysisFlow2->Fill(kTOF2t);
            }else if (nTofTrks!=2)
            {
                return;
            }

            if(Eta[0] < -0.7 || Eta[0] > 0.7 || Eta[1] < -0.7 || Eta[1] > 0.7)
            {
                int tmp1 = 0;
                int tmp2 = 1;
                if(Eta[1] < -0.7 || Eta[1] > 0.7)
                {
                    tmp1 = 1;
                    tmp2 = 0;
                }    

                hEtaTest[0]->Fill(Eta[tmp1]);
                hEtaTest[1]->Fill(Eta[tmp2]);
            }
            if(!fourPiState)
                hAnalysisFlow->Fill(kTOF2t);

            if( vertexId[0] != vertexId[1])
                return;

            if(fourPiState)
            {
                if(vertexId[2] != vertexId[3] || vertexId[2] != vertexId[1])
                    return;
                hAnalysisFlow2->Fill(kSameVrtx);
                FillPlots(1,i,1);
            }else
            {
                hAnalysisFlow->Fill(kSameVrtx);
                FillPlots(1,i,0);
            }
      
            for (int iSide = 0; iSide < nSides; ++iSide)
            {
                StUPCRpsTrack* trackRP = rpEvt->getTrack(rpTrackIdVec_perSide[iSide][0]);
                if(trackRP)
                {   
                    StUPCRpsTrackPoint* trackPoint = trackRP->getTrackPoint(0);
                    if(trackPoint)
                    {
                        rpX[iSide] = trackPoint->x(); 
                        rpY[iSide] = trackPoint->y(); 
                        rpZ[iSide] = trackPoint->z();
                    }else
                    {
                        hInvalidRPtrack->Fill(3);
                        return;
                    }
                    xCorrelationsRp[iSide] = trackRP->pVec().X();
                    yCorrelationsRp[iSide] = trackRP->pVec().Y();
                    thetaRp[iSide] = trackRP->thetaRp();
                    phiRp[iSide] = trackRP->phiRp(); 
                    timeRp[iSide] = trackRP->time();
                    pRp[iSide] = trackRP->p();
                    ptRp[iSide] = trackRP->pt();
                    etaRp[iSide] = trackRP->eta();
                    t[iSide] = trackRP->t(beamMomentum);
                    xi[iSide] = trackRP->xi(beamMomentum);
                }else
                {
                    hInvalidRPtrack->Fill(1);
                    return;
                }
            }

            for(int iPart = 0; iPart < nParticles; ++iPart)
                chiPair[iPart] = nSigmaTPC[iPart][0]*nSigmaTPC[iPart][0] + nSigmaTPC[iPart][1]*nSigmaTPC[iPart][1];

            for (int iRP = 0; iRP < nRomanPots; ++iRP)
            {
                ADC[iRP][0] = rpEvt->adc(iRP,0);
                ADC[iRP][1] = rpEvt->adc(iRP,1);
                TAC[iRP][0] = rpEvt->tac(iRP,0);
                TAC[iRP][1] = rpEvt->tac(iRP,1);
            }

            runNumber = upcEvt->getRunNumber();
            BBCLarge[East] = upcEvt->getBBCLargeEast();
            BBCSmall[East] = upcEvt->getBBCSmallEast();
            BBCLarge[West] = upcEvt->getBBCLargeWest();
            BBCSmall[West] = upcEvt->getBBCSmallWest();
            VPDTimeDiff = upcEvt->getVPDTimeDiff();
            VPDSumWest = upcEvt->getVPDSumWest();
            VPDSumEast = upcEvt->getVPDSumEast();

            if(!fourPiState)
               CalculatePID();

            if(totalCharge != 0){
                if( missingPt < 0.1 )
                    bcgTree->Fill();
                return;
            }

            if(fourPiState)
            {
                hAnalysisFlow2->Fill(kTotCH0);
                FillPlots(2,i,1);
            } else
            {
                hAnalysisFlow->Fill(kTotCH0);
                FillPlots(2,i,0);      
            }

            if( missingPt > 0.1 )
                return;

            if(fourPiState)
            {
                hAnalysisFlow2->Fill(kMissPt);
                FillPlots(3,i,1);
            } else
            {
                hAnalysisFlow->Fill(kMissPt);
                FillPlots(3,i,0);      
            }

            recTree->Fill();
        }
    } // end of loop over arms
}


//_____________________________________________________________________________
TFile *CreateOutputTree(const string& out) {

    TFile *outputFile = TFile::Open(out.c_str(), "recreate");
    if(!outputFile) 
        return 0x0;

    //standard reconstructed tree
    recTree = new TTree("recTree", "recTree");

// PID and some quality event info
    recTree->Branch("missingPt", &missingPt, "missingPt/D");
    recTree->Branch("deltaTOF", &deltaTOF, "deltaTOF/D");
    recTree->Branch("mSquared", &mSquared, "mSquared/D"); 
    recTree->Branch("nSigTrk1Pion", &nSigmaTPC[Pion][0], "nSigTrk1Pion/D");
    recTree->Branch("nSigTrk2Pion", &nSigmaTPC[Pion][1], "nSigTrk2Pion/D");
    for (int iPart = 0; iPart < nParticles; ++iPart)
    {
        recTree->Branch("invMass" + particleLables[iPart], &invMass[iPart], "invMass"  + particleLables[iPart] + "/D");
        recTree->Branch("chiPair" + particleLables[iPart], &chiPair[iPart], "chiPair"  + particleLables[iPart] + "/D");
        recTree->Branch("deltaTOFExpected" + particleLables[iPart], &deltaTOFExpected[iPart], "deltaTOFExpected"  + particleLables[iPart] + "/D");
        recTree->Branch("deltaDeltaTOF" + particleLables[iPart], &deltaDeltaTOF[iPart], "deltaDeltaTOF"  + particleLables[iPart] + "/D");  
    }


// Vertex info
    recTree->Branch("vertexZ", &vertexesZ[0], "vertexZ/D");

// Central track info
    for (int i = 0; i < 4; ++i)
    {
        recTree->Branch(Form("dEdx%i",i), &dEdx[i], Form("dEdx%i/D",i));
        recTree->Branch(Form("momentum%i",i), &momentum[i], Form("momentum%i/D",i));
        recTree->Branch(Form("transMomentum%i",i), &transMomentum[i], Form("transMomentum%i/D",i));
        recTree->Branch(Form("charge%i",i), &charge[i], Form("charge%i/D",i));
        recTree->Branch(Form("TOFtime%i",i), &TOFtime[i], Form("TOFtime%i/D",i));
        recTree->Branch(Form("TOFlength%i",i), &TOFlength[i], Form("TOFlength%i/D",i));
        recTree->Branch(Form("DcaXY%i",i), &DcaXY[i], Form("DcaXY%i/D",i));
        recTree->Branch(Form("DcaZ%i",i), &DcaZ[i], Form("DcaZ%i/D",i));
        recTree->Branch(Form("NhitsFit%i",i), &NhitsFit[i], Form("NhitsFit%i/D",i));
        recTree->Branch(Form("NhitsDEdx%i",i), &NhitsDEdx[i], Form("NhitsDEdx%i/D",i));
        recTree->Branch(Form("Eta%i",i), &Eta[i], Form("Eta%i/D",i));
        recTree->Branch(Form("Phi%i",i), &Phi[i], Form("Phi%i/D",i));
        recTree->Branch(Form("Chi2%i",i), &Chi2[i], Form("Chi2%i/D",i));
    }
// RP track info  
    for (int i = 0; i < nSides; ++i)
    {
        recTree->Branch("rpX" + sideLabel[i], &rpX[i], "rpX" + sideLabel[i] + "/D");
        recTree->Branch("rpY" + sideLabel[i], &rpY[i], "rpY" + sideLabel[i] + "/D");
        recTree->Branch("rpZ" + sideLabel[i], &rpZ[i], "rpZ" + sideLabel[i] + "/D");
        recTree->Branch("thetaRp" + sideLabel[i], &thetaRp[i], "thetaRp" + sideLabel[i] + "/D");
        recTree->Branch("phiRp" + sideLabel[i], &phiRp[i], "phiRp" + sideLabel[i] + "/D");
        recTree->Branch("timeRp" + sideLabel[i], &timeRp[i], "time" + sideLabel[i] + "/D");
        recTree->Branch("pRp" + sideLabel[i], &pRp[i], "p" + sideLabel[i] + "/D");
        recTree->Branch("ptRp" + sideLabel[i], &ptRp[i], "pt" + sideLabel[i] + "/D");
        recTree->Branch("etaRp" + sideLabel[i], &etaRp[i], "eta" + sideLabel[i] + "/D");
        recTree->Branch("xCorrelationsRp" + sideLabel[i], &xCorrelationsRp[i], "xCorrelations" + sideLabel[i] + "/D");
        recTree->Branch("yCorrelationsRp" + sideLabel[i], &yCorrelationsRp[i], "yCorrelations" + sideLabel[i] + "/D");
        recTree->Branch("t" + sideLabel[i], &t[i], "t" + sideLabel[i] + "/D");
        recTree->Branch("xi" + sideLabel[i], &xi[i], "xi" + sideLabel[i] + "/D");
    }
// RP event info
    for (int i = 0; i < nRomanPots; ++i)
    {
        recTree->Branch("ADC_" + rpNames[i] + "V", &ADC[i][0], "ADC_" + rpNames[i] + "V/D");
        recTree->Branch("ADC_" + rpNames[i] + "H", &ADC[i][1], "ADC_" + rpNames[i] + "H/D");
        recTree->Branch("TAC_" + rpNames[i] + "V", &TAC[i][0], "TAC_" + rpNames[i] + "V/D");
        recTree->Branch("TAC_" + rpNames[i] + "H", &TAC[i][1], "TAC_" + rpNames[i] + "H/D");
    }

// event info
    recTree->Branch("elastic", &elastic, "elastic/O");
    recTree->Branch("fourPiState", &fourPiState, "fourPiState/O");
    recTree->Branch("runNumber", &runNumber, "runNumber/i");
    recTree->Branch("VPDTimeDiff", &VPDTimeDiff, "VPDTimeDiff/D");
    recTree->Branch("VPDSumWest", &VPDSumWest, "VPDSumWest/D");
    recTree->Branch("VPDSumEast", &VPDSumEast, "VPDSumEast/D");
    for (int i = 0; i < nSides; ++i)
    {
        recTree->Branch("BBCSmall" + sideLabel[i], &BBCSmall[i], "BBCSmall" + sideLabel[i] + "/i");
        recTree->Branch("BBCLarge" + sideLabel[i], &BBCLarge[i], "BBCLarge" + sideLabel[i] + "/i"); 
    }
    recTree->Branch("RP_CPT2_570701", &trigger[3], "RP_CPT2_570701/O");
    recTree->Branch("RP_CPT2noBBCL_570705", &trigger[7], "RP_CPT2noBBCL_570705/O");
    recTree->Branch("RP_CPT2_570711", &trigger[9], "RP_CPT2_570711/O");
    recTree->Branch("RP_CPT2_590701", &trigger[12], "RP_CPT2_590701/O");
    recTree->Branch("RP_CPT2noBBCL_590705", &trigger[14], "RP_CPT2noBBCL_590705/O");
    recTree->Branch("RP_CPTnoBBCL_590708", &trigger[15], "RP_CPTnoBBCL_590708/O");

// Setting background Tree
    bcgTree = recTree->CloneTree(0);
    bcgTree->SetName("Background");

    return outputFile;

}//CreateOutputTree

bool ConnectInput(const string& in, int fileId) {
    //input from file or chain
    upcTree = 0x0;
    if( in.find(".root") != string::npos ) 
    {
        cout << "Input from root file" << endl;
        infile = TFile::Open(in.c_str(), "read");
        if(!infile) 
            return false;
        upcTree = dynamic_cast<TTree*>( infile->Get("mUPCTree") );
    } else 
    {
        cout << "Input from chain" << endl;
        upcChain = new TChain("mUPCTree");
        ifstream instr(in.c_str());
        string line;
        int lineId=0;
        while(getline(instr, line)) 
        {
            if(fileId==lineId || fileId== -1)
                upcChain->AddFile(line.c_str());
            lineId++;
        }
        instr.close();
        upcTree = dynamic_cast<TTree*>( upcChain );
    }

    if(!upcTree) 
        return false;

    upcTree->SetBranchAddress("mUPCEvent", &upcEvt);
    upcTree->SetBranchAddress("mRPEvent", &rpEvt); 

    return true;

}//ConnectInput

void FillPlots(int state, int configuration, int isFourPiState){
    for(int tmp = 0; tmp < 3; ++tmp)
    {
        if(tmp == 1 && (configuration!=IUU && configuration!=IDD))
            continue;

        if(tmp == 2 && (configuration!=EUD && configuration!=EDU))
            continue;

        for(int trackId = 0; trackId < nTofTrks; ++trackId)
            hdEdxVsqP[4*tmp+state][isFourPiState]->Fill(momentum[trackId]*charge[trackId],dEdx[trackId]);

        hMissingPt[4*tmp+state][isFourPiState]->Fill(missingPt); 
    }
}

void Clear(){
    totalCharge = 0;
    nTofTrks = 0;
    runNumber = 0;

    for (int i = 0; i < nSides; ++i)
    {
        xCorrelationsRp[i] = -999;
        yCorrelationsRp[i] = -999;
        thetaRp[i] = -999;
        phiRp[i] = -999;
        timeRp[i] = -999;
        pRp[i] = -999;
        ptRp[i] = -999;
        etaRp[i] = -999;
        rpX[i] = -999;
        rpY[i] = -999;
        rpZ[i] = -999;
    }

    for (int i = 0; i < 4; ++i)
    {
        dEdx[i] = -999;
        momentum[i] = -999;
        transMomentum[i] = -999;
        TOFtime[i] = -999;
        TOFlength[i] = -999;
        charge[i] = -999;
        nSigmaTPC[Pion][i] = -999;
        nSigmaTPC[Kaon][i] = -999;
        nSigmaTPC[Proton][i] = -999;
        vertexId[i] = -999;
        vertexesZ[i] = -999;
        DcaXY[i] = -999;
        DcaZ[i] = -999;
        NhitsFit[i] = -999;
        NhitsDEdx[i] = -999;
        Eta[i] = -999;
        Phi[i] = -999;
        Chi2[i] = -999;    
    }

    for (int i = 0; i < nTriggers; ++i)
    {
        trigger[i] = false;
    }
}



void CalculatePID(){

    if(TOFtime[0] < 0 || TOFtime[1] < 0 || TOFlength[0] < 0 || TOFlength[1] < 0)
    {
        mSquared = -999.0;
        deltaTOF = -999;   
        for (int iPart = 0; iPart < nParticles; ++iPart)
        {
            deltaTOFExpected[iPart] = -999;
            deltaDeltaTOF[iPart] = -999;
        }

        return;
    }
    deltaTOF = TOFtime[1] - TOFtime[0];

    double speedOfLight2 = speedOfLight*speedOfLight;
    double speedOfLight4 = speedOfLight2*speedOfLight2;
    double length1Squared = TOFlength[0]*TOFlength[0]/(100*100); // convert TOFlength from cm to m
    double length2Squared = TOFlength[1]*TOFlength[1]/(100*100); // convert TOFlength from cm to m
    double deltaTime2 = (deltaTOF*deltaTOF)/(pow(10.0,18.0)); // convert TOFtime from ns to s
    double deltaTime4 = deltaTime2*deltaTime2;
    double oneOverMomentum1sq = 1/(momentum[0]*momentum[0]);
    double oneOverMomentum2sq = 1/(momentum[1]*momentum[1]);
    double cEq = -2*length1Squared*length2Squared + speedOfLight4*deltaTime4 + length2Squared*length2Squared + length1Squared*length1Squared -2*speedOfLight2*deltaTime2*(length2Squared + length1Squared);
    double bEq = -2*length1Squared*length2Squared*(oneOverMomentum1sq + oneOverMomentum2sq) + 2*length1Squared*length1Squared*oneOverMomentum1sq + 2*length2Squared*length2Squared*oneOverMomentum2sq -2*speedOfLight2*deltaTime2*(length1Squared*oneOverMomentum1sq + length2Squared*oneOverMomentum2sq);
    double aEq = -2*length1Squared*length2Squared*oneOverMomentum1sq*oneOverMomentum2sq + length1Squared*length1Squared*oneOverMomentum1sq*oneOverMomentum1sq + length2Squared*length2Squared*oneOverMomentum2sq*oneOverMomentum2sq;
    mSquared = (-bEq + sqrt(bEq*bEq-4*aEq*cEq)) / (2*aEq);
    
//////////////////////
      
/////////////// Calculate deltaDeltaTOFPion
    double expectedTime1, expectedTime2;

    for (int i = 0; i < nParticles; ++i)
    {
        expectedTime1 = (TOFlength[0] / speedOfLight) * sqrt(1 + pow(particleMass[i]/ momentum[0], 2)) * pow(10,7); // in ns
        expectedTime2 = (TOFlength[1] / speedOfLight) * sqrt(1 + pow(particleMass[i]/ momentum[1],2)) * pow(10,7); // in ns             
        deltaTOFExpected[i] = expectedTime2 - expectedTime1;        
        deltaDeltaTOF[i] = deltaTOF - deltaTOFExpected[i];
    }

}
        

