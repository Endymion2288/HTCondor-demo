#define ReadTree_cxx
//#define SHOW_DEBUG
#ifdef SHOW_DEBUG
#define VERBOSE
#endif

//#define DRAW_RAW

//#define CUT_GLOBAL_VTX_PROB
#define CUT_FROM_3J_MUON

#define CUT_FROM_3J_JPSI

//#define CUT_MUON_ID_LOOSE

#define CUT_MUON_ID_SOFT

//#define CUT_UPS_TRY

#define ALLOW_OVERLAP

#define NO_PLOT

#include "/afs/cern.ch/user/c/chiw/cernbox/JpsiJpsiUps/preCut/ParticleCand/interface/ReadTree.h"
#include "/afs/cern.ch/user/c/chiw/cernbox/JpsiJpsiUps/preCut/ParticleCand/src/ParticleCand.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <memory>
#include <cstdio>
// Include the header file for the roofit.
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooChebychev.h"

void ReadTree::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    printf("Entries: %ld\n", nentries);

    const unsigned int nBin = 40;
    const unsigned int nBin_cut = 20;
    const unsigned int nCandsAllowed = 20;

    // Use Roofit to draw the plot with proper error bars.
    // Define mass histograms for Jpsi, Ups and Pri passing the cut. Using Roofit.
    RooRealVar Jpsi1_mass_cut_var("Jpsi1_mass_cut", "Jpsi1_mass_cut", 2.5, 3.5);
    RooRealVar Jpsi2_mass_cut_var("Jpsi2_mass_cut", "Jpsi2_mass_cut", 2.5, 3.5);
    RooRealVar Ups_mass_cut_var("Ups_mass_cut","Ups_mass_cut", 8.0, 12.0);
    RooRealVar Pri_mass_cut_var("Pri_mass_cut","Pri_mass_cut", 0.0, 100.0);

    // Define dataset for Jpsi, Ups and Pri passing the cut. Using Roofit.
    RooDataSet Jpsi1_mass_cut_data("Jpsi1_mass_cut_data", "Jpsi1_mass_cut_data", RooArgList(Jpsi1_mass_cut_var));
    RooDataSet Jpsi2_mass_cut_data("Jpsi2_mass_cut_data", "Jpsi2_mass_cut_data", RooArgList(Jpsi2_mass_cut_var));
    RooDataSet Ups_mass_cut_data("Ups_mass_cut_data", "Ups_mass_cut_data", RooArgList(Ups_mass_cut_var));
    RooDataSet Pri_mass_cut_data("Pri_mass_cut_data", "Pri_mass_cut_data", RooArgList(Pri_mass_cut_var));


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        // Initialize the total set of candidates.
        std::vector<std::shared_ptr<ParticleCand> > CandList;
        ParticleCand tempCand;
        ParticleCand::PartIdxList_t tempList;
        double temp_massChi2;

        const double vtxProb_Thr = 0.05;

        // Marker
        #ifdef SHOW_DEBUG
        printf("\n>>>>> Begin new event %lld <<<<<\n", jentry);
        #endif

        // Display event number once 500 events.
        if(jentry % 500 == 0){
            printf("Processing event %lld\n", jentry);
        }
        // Loop over all candidates.
        size_t nCands = Jpsi_1_mass->size();

        #ifdef VERBOSE
        printf("number of candidate = %lld\n", nCands);
        #endif

        if(nCands == 0) continue;

        // Calculate muon pT, eta, phi from px, py and pz.
        std::vector<double> mu_pT;
        std::vector<double> mu_eta;
        std::vector<double> mu_phi;
        
        for(Long64_t iMuon=0; iMuon < nMu; iMuon++){
            // Calculate pT, eta, phi from px, py and pz.
            double pT  = sqrt(muPx->at(iMuon) * muPx->at(iMuon) + muPy->at(iMuon) * muPy->at(iMuon));
            double eta = 0.5 * log((pT + muPz->at(iMuon)) / (pT - muPz->at(iMuon)));
            double phi = atan2(muPy->at(iMuon), muPx->at(iMuon));
            mu_pT.push_back(pT);
            mu_eta.push_back(eta);
            mu_phi.push_back(phi);
        }

        for (Long64_t iCand=0; iCand < nCands; iCand++) {
            // Exclude failed fitting.
            if(Pri_mass->at(iCand) <= 0.0) continue;

            #ifdef CUT_GLOBAL_VTX_PROB
            // Cut from vertex probability
            if(Pri_VtxProb->at(iCand) < vtxProb_Thr) continue;
            if(Jpsi_1_VtxProb->at(iCand) < vtxProb_Thr) continue;
            if(Jpsi_2_VtxProb->at(iCand) < vtxProb_Thr) continue;
            if(Ups_VtxProb->at(iCand) < vtxProb_Thr) continue;
            #endif

            // Further cuts are applied here.
            bool passCut = true;

            // Prevent underflow or overflow of masses.
            if(Jpsi_1_mass->at(iCand) < 2.5 || Jpsi_1_mass->at(iCand) > 3.5){
                continue;
            }
            if(Jpsi_2_mass->at(iCand) < 2.5 || Jpsi_2_mass->at(iCand) > 3.5){
                continue;
            }
            if(Ups_mass->at(iCand) < 8.0 || Ups_mass->at(iCand) > 12.0){
                continue;
            }

            // Register the candidate final state muons.
            tempList.push_back(Jpsi_1_mu_1_Idx->at(iCand));
            tempList.push_back(Jpsi_1_mu_2_Idx->at(iCand));
            tempList.push_back(Jpsi_2_mu_1_Idx->at(iCand));
            tempList.push_back(Jpsi_2_mu_2_Idx->at(iCand));
            tempList.push_back(Ups_mu_1_Idx->at(iCand));
            tempList.push_back(Ups_mu_2_Idx->at(iCand));
            tempCand.AddParticle(ParticleCand::PartType::Muon, tempList);

            #ifdef CUT_MUON_ID_LOOSE
            // Check by muIsPatLooseMuon
            for(auto idx : tempList){
                if(!muIsPatLooseMuon->at(idx)){
                    passCut = false;
                    break;
                }
            }

            #endif

            #ifdef CUT_MUON_ID_SOFT
            // Check by muIsPatSoftMuon
            for(auto idx : tempList){
                if(!muIsPatSoftMuon->at(idx)){
                    passCut = false;
                    break;
                }
            }

            #endif

            if (!passCut){
                continue;
            }

            #ifdef CUT_FROM_3J_MUON
            // Cut from 3J
            // Muons that decay from Jpsi and Ups:
            // require pT > 3.5GeV/c for abs(eta) < 1.2 and pT > 2.5GeV/c for 1.2 < abs(eta) < 2.4
            for(auto idx : tempList){
                if(mu_pT[idx] < 3.5 && abs(mu_eta[idx]) < 1.2){
                    passCut = false;
                    break;
                }
                if(mu_pT[idx] < 2.5 && abs(mu_eta[idx]) > 1.2 && abs(mu_eta[idx]) < 2.4){
                    passCut = false;
                    break;
                }
                if(abs(mu_eta[idx]) > 2.4){
                    passCut = false;
                    break;
                }
            }

            if(!passCut){
                continue;
            }
            #endif

            #ifdef CUT_FROM_3J_JPSI

            // For Jpsi: require pT > 6GeV/c and abs(eta) < 2.4
            // Trying as precut: lower the pT cut to 2GeV/c
            if(Jpsi_1_pt->at(iCand) <= 6.0 || abs(Jpsi_1_eta->at(iCand)) > 2.4){
                continue ;
            }
            if(Jpsi_2_pt->at(iCand) <= 6.0 || abs(Jpsi_2_eta->at(iCand)) > 2.4){
                continue ;
            }


            #endif

            #ifdef CUT_UPS_TRY
            // For Upsilon: try a loose cut. pT > 6GeV/c and abs(eta) < 2.4
            if(Ups_pt->at(iCand) <= 6.0 || abs(Ups_eta->at(iCand)) > 2.4){
                passCut = false;
            }
            // Additional cut for muons from Ups: pT > 4 GeV/c
            if(mu_pT[Ups_mu_1_Idx->at(iCand)] <= 4.0 || mu_pT[Ups_mu_2_Idx->at(iCand)] <= 4.0){
                passCut = false;
            }
            if(!muIsPatMediumMuon->at(Ups_mu_1_Idx->at(iCand))){
                passCut = false;
            }
            if(!muIsPatMediumMuon->at(Ups_mu_2_Idx->at(iCand))){
                passCut = false;
	        }
            #endif

            if(!passCut){
                continue;
            }


            // Calculate Chi2 from massDiff and massErr of Jpsi and Ups.
            temp_massChi2 =   (Jpsi_1_massDiff->at(iCand) / Jpsi_1_massErr->at(iCand))
                            * (Jpsi_1_massDiff->at(iCand) / Jpsi_1_massErr->at(iCand))
                            + (Jpsi_2_massDiff->at(iCand) / Jpsi_2_massErr->at(iCand))
                            * (Jpsi_2_massDiff->at(iCand) / Jpsi_2_massErr->at(iCand))
                            + (Ups_massDiff->at(iCand) / Ups_massErr->at(iCand))
                            * (Ups_massDiff->at(iCand) / Ups_massErr->at(iCand));
            tempCand.SetScore(temp_massChi2);
            tempCand.SetId(iCand);
            CandList.push_back(std::make_shared<ParticleCand>(tempCand));
            tempList.clear();
            tempCand.Clear();
        }
        // Sort all candidates for later use.
        std::sort(CandList.begin(), CandList.end(), [](const std::shared_ptr<ParticleCand>& a, 
                                                       const std::shared_ptr<ParticleCand>& b){
            return a->GetScore() < b->GetScore();
        });

        // Print out the sorted candidates.
        #ifdef VERBOSE
        for(auto& cand : CandList){
            printf("[ %d ] %.4f ",cand->GetId() ,cand->GetScore());
            puts(cand->ToString().c_str());
        }
        #endif

        #ifdef ALLOW_OVERLAP
        std::vector<std::shared_ptr<ParticleCand> > SelectedCands_cut;
        for(auto& cand : CandList){
            if(SelectedCands_cut.size() < nCandsAllowed){
                SelectedCands_cut.push_back(cand);
            }
            else{
                break;
            }
        }

        #else
        std::vector<std::shared_ptr<ParticleCand> > SelectedCands_cut;
        for(auto& cand : CandList){
            if(SelectedCands_cut.empty()){
                SelectedCands_cut.push_back(cand);
            }
            else{
                bool isOverlap = false;
                for(auto& selCand : SelectedCands_cut){
                    if(cand->Overlap(*selCand)){
                        isOverlap = true;
                        break;
                    }
                }
                if(!isOverlap){
                    SelectedCands_cut.push_back(cand);
                }
            }
        }
        #endif

        // Print out the selected candidates with the highest score.
        #ifdef VERBOSE
        printf("\n Selected candidates: %lld\n", SelectedCands_raw.size());
        for(auto& cand : SelectedCands_raw){
            printf("[ %d ] %.4f ",cand->GetId() ,cand->GetScore());
            puts(cand->ToString().c_str());
        }
        #endif

        // Draw the selected candidates passing the cut.
        for(auto& cand : SelectedCands_cut){
            // Set the mass.
            Jpsi1_mass_cut_var.setVal(Jpsi_1_mass->at(cand->GetId()));
            Jpsi2_mass_cut_var.setVal(Jpsi_2_mass->at(cand->GetId()));
            Ups_mass_cut_var.setVal(Ups_mass->at(cand->GetId()));
            Pri_mass_cut_var.setVal(Pri_mass->at(cand->GetId()));
            // Set the mass error.
            Jpsi1_mass_cut_var.setError(Jpsi_1_massErr->at(cand->GetId()));
            Jpsi2_mass_cut_var.setError(Jpsi_2_massErr->at(cand->GetId()));
            Ups_mass_cut_var.setError(Ups_massErr->at(cand->GetId()));
            Pri_mass_cut_var.setError(Pri_massErr->at(cand->GetId()));
            // Add the data to the dataset.
            Jpsi1_mass_cut_data.add(RooArgSet(Jpsi1_mass_cut_var));
            Jpsi2_mass_cut_data.add(RooArgSet(Jpsi2_mass_cut_var));
            Ups_mass_cut_data.add(RooArgSet(Ups_mass_cut_var));
            Pri_mass_cut_data.add(RooArgSet(Pri_mass_cut_var));
        }

        // Use the pre-defined outputTree to store the filtered data.
        for(auto& cand : SelectedCands_cut){
            filtered_Jpsi_1_mass->push_back(Jpsi_1_mass->at(cand->GetId()));
            filtered_Jpsi_1_massErr->push_back(Jpsi_1_massErr->at(cand->GetId()));
            filtered_Jpsi_1_massDiff->push_back(Jpsi_1_massDiff->at(cand->GetId()));
            filtered_Jpsi_1_ctau->push_back(Jpsi_1_ctau->at(cand->GetId()));
            filtered_Jpsi_1_ctauErr->push_back(Jpsi_1_ctauErr->at(cand->GetId()));
            filtered_Jpsi_1_Chi2->push_back(Jpsi_1_Chi2->at(cand->GetId()));
            filtered_Jpsi_1_ndof->push_back(Jpsi_1_ndof->at(cand->GetId()));
            filtered_Jpsi_1_VtxProb->push_back(Jpsi_1_VtxProb->at(cand->GetId()));
            filtered_Jpsi_1_px->push_back(Jpsi_1_px->at(cand->GetId()));
            filtered_Jpsi_1_py->push_back(Jpsi_1_py->at(cand->GetId()));
            filtered_Jpsi_1_pz->push_back(Jpsi_1_pz->at(cand->GetId()));
            filtered_Jpsi_1_phi->push_back(Jpsi_1_phi->at(cand->GetId()));
            filtered_Jpsi_1_eta->push_back(Jpsi_1_eta->at(cand->GetId()));
            filtered_Jpsi_1_pt->push_back(Jpsi_1_pt->at(cand->GetId()));
            filtered_Jpsi_1_mu_1_Idx->push_back(Jpsi_1_mu_1_Idx->at(cand->GetId()));
            filtered_Jpsi_1_mu_2_Idx->push_back(Jpsi_1_mu_2_Idx->at(cand->GetId()));

            filtered_Jpsi_2_mass->push_back(Jpsi_2_mass->at(cand->GetId()));
            filtered_Jpsi_2_massErr->push_back(Jpsi_2_massErr->at(cand->GetId()));
            filtered_Jpsi_2_massDiff->push_back(Jpsi_2_massDiff->at(cand->GetId()));
            filtered_Jpsi_2_ctau->push_back(Jpsi_2_ctau->at(cand->GetId()));
            filtered_Jpsi_2_ctauErr->push_back(Jpsi_2_ctauErr->at(cand->GetId()));
            filtered_Jpsi_2_Chi2->push_back(Jpsi_2_Chi2->at(cand->GetId()));
            filtered_Jpsi_2_ndof->push_back(Jpsi_2_ndof->at(cand->GetId()));
            filtered_Jpsi_2_VtxProb->push_back(Jpsi_2_VtxProb->at(cand->GetId()));
            filtered_Jpsi_2_px->push_back(Jpsi_2_px->at(cand->GetId()));
            filtered_Jpsi_2_py->push_back(Jpsi_2_py->at(cand->GetId()));
            filtered_Jpsi_2_pz->push_back(Jpsi_2_pz->at(cand->GetId()));
            filtered_Jpsi_2_phi->push_back(Jpsi_2_phi->at(cand->GetId()));
            filtered_Jpsi_2_eta->push_back(Jpsi_2_eta->at(cand->GetId()));
            filtered_Jpsi_2_pt->push_back(Jpsi_2_pt->at(cand->GetId()));
            filtered_Jpsi_2_mu_1_Idx->push_back(Jpsi_2_mu_1_Idx->at(cand->GetId()));
            filtered_Jpsi_2_mu_2_Idx->push_back(Jpsi_2_mu_2_Idx->at(cand->GetId()));

            filtered_Pri_mass->push_back(Pri_mass->at(cand->GetId()));
            filtered_Pri_massErr->push_back(Pri_massErr->at(cand->GetId()));
            filtered_Pri_ctau->push_back(Pri_ctau->at(cand->GetId()));
            filtered_Pri_ctauErr->push_back(Pri_ctauErr->at(cand->GetId()));
            filtered_Pri_Chi2->push_back(Pri_Chi2->at(cand->GetId()));
            filtered_Pri_ndof->push_back(Pri_ndof->at(cand->GetId()));
            filtered_Pri_VtxProb->push_back(Pri_VtxProb->at(cand->GetId()));
            filtered_Pri_px->push_back(Pri_px->at(cand->GetId()));
            filtered_Pri_py->push_back(Pri_py->at(cand->GetId()));
            filtered_Pri_pz->push_back(Pri_pz->at(cand->GetId()));
            filtered_Pri_phi->push_back(Pri_phi->at(cand->GetId()));
            filtered_Pri_eta->push_back(Pri_eta->at(cand->GetId()));
            filtered_Pri_pt->push_back(Pri_pt->at(cand->GetId()));

            filtered_Ups_mass->push_back(Ups_mass->at(cand->GetId()));
            filtered_Ups_massErr->push_back(Ups_massErr->at(cand->GetId()));
            filtered_Ups_massDiff->push_back(Ups_massDiff->at(cand->GetId()));
            filtered_Ups_Chi2->push_back(Ups_Chi2->at(cand->GetId()));
            filtered_Ups_ndof->push_back(Ups_ndof->at(cand->GetId()));
            filtered_Ups_VtxProb->push_back(Ups_VtxProb->at(cand->GetId()));
            filtered_Ups_px->push_back(Ups_px->at(cand->GetId()));
            filtered_Ups_py->push_back(Ups_py->at(cand->GetId()));
            filtered_Ups_pz->push_back(Ups_pz->at(cand->GetId()));
            filtered_Ups_phi->push_back(Ups_phi->at(cand->GetId()));
            filtered_Ups_eta->push_back(Ups_eta->at(cand->GetId()));
            filtered_Ups_pt->push_back(Ups_pt->at(cand->GetId()));
            filtered_Ups_mu_1_Idx->push_back(Ups_mu_1_Idx->at(cand->GetId()));
            filtered_Ups_mu_2_Idx->push_back(Ups_mu_2_Idx->at(cand->GetId()));

            // For the muons
            filtered_Jpsi_1_mu_1_px->push_back(muPx->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_1_py->push_back(muPy->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_1_pz->push_back(muPz->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_1_eta->push_back(mu_eta[Jpsi_1_mu_1_Idx->at(cand->GetId())]);
            filtered_Jpsi_1_mu_1_phi->push_back(mu_phi[Jpsi_1_mu_1_Idx->at(cand->GetId())]);
            filtered_Jpsi_1_mu_1_pt->push_back(mu_pT[Jpsi_1_mu_1_Idx->at(cand->GetId())]);
            filtered_Jpsi_1_mu_1_isPatLooseMuon->push_back(muIsPatLooseMuon->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_1_isPatSoftMuon->push_back(muIsPatSoftMuon->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_1_isPatMediumMuon->push_back(muIsPatMediumMuon->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_1_isPatTightMuon->push_back(muIsPatTightMuon->at(Jpsi_1_mu_1_Idx->at(cand->GetId())));

            filtered_Jpsi_1_mu_2_px->push_back(muPx->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_2_py->push_back(muPy->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_2_pz->push_back(muPz->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_2_eta->push_back(mu_eta[Jpsi_1_mu_2_Idx->at(cand->GetId())]);
            filtered_Jpsi_1_mu_2_phi->push_back(mu_phi[Jpsi_1_mu_2_Idx->at(cand->GetId())]);
            filtered_Jpsi_1_mu_2_pt->push_back(mu_pT[Jpsi_1_mu_2_Idx->at(cand->GetId())]);
            filtered_Jpsi_1_mu_2_isPatLooseMuon->push_back(muIsPatLooseMuon->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_2_isPatSoftMuon->push_back(muIsPatSoftMuon->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_2_isPatMediumMuon->push_back(muIsPatMediumMuon->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_1_mu_2_isPatTightMuon->push_back(muIsPatTightMuon->at(Jpsi_1_mu_2_Idx->at(cand->GetId())));

            filtered_Jpsi_2_mu_1_px->push_back(muPx->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_1_py->push_back(muPy->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_1_pz->push_back(muPz->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_1_eta->push_back(mu_eta[Jpsi_2_mu_1_Idx->at(cand->GetId())]);
            filtered_Jpsi_2_mu_1_phi->push_back(mu_phi[Jpsi_2_mu_1_Idx->at(cand->GetId())]);
            filtered_Jpsi_2_mu_1_pt->push_back(mu_pT[Jpsi_2_mu_1_Idx->at(cand->GetId())]);
            filtered_Jpsi_2_mu_1_isPatLooseMuon->push_back(muIsPatLooseMuon->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_1_isPatSoftMuon->push_back(muIsPatSoftMuon->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_1_isPatMediumMuon->push_back(muIsPatMediumMuon->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_1_isPatTightMuon->push_back(muIsPatTightMuon->at(Jpsi_2_mu_1_Idx->at(cand->GetId())));

            filtered_Jpsi_2_mu_2_px->push_back(muPx->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_2_py->push_back(muPy->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_2_pz->push_back(muPz->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_2_eta->push_back(mu_eta[Jpsi_2_mu_2_Idx->at(cand->GetId())]);
            filtered_Jpsi_2_mu_2_phi->push_back(mu_phi[Jpsi_2_mu_2_Idx->at(cand->GetId())]);
            filtered_Jpsi_2_mu_2_pt->push_back(mu_pT[Jpsi_2_mu_2_Idx->at(cand->GetId())]);
            filtered_Jpsi_2_mu_2_isPatLooseMuon->push_back(muIsPatLooseMuon->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_2_isPatSoftMuon->push_back(muIsPatSoftMuon->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_2_isPatMediumMuon->push_back(muIsPatMediumMuon->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));
            filtered_Jpsi_2_mu_2_isPatTightMuon->push_back(muIsPatTightMuon->at(Jpsi_2_mu_2_Idx->at(cand->GetId())));

            filtered_Ups_mu_1_px->push_back(muPx->at(Ups_mu_1_Idx->at(cand->GetId())));
            filtered_Ups_mu_1_py->push_back(muPy->at(Ups_mu_1_Idx->at(cand->GetId())));
            filtered_Ups_mu_1_pz->push_back(muPz->at(Ups_mu_1_Idx->at(cand->GetId())));
            filtered_Ups_mu_1_eta->push_back(mu_eta[Ups_mu_1_Idx->at(cand->GetId())]);
            filtered_Ups_mu_1_phi->push_back(mu_phi[Ups_mu_1_Idx->at(cand->GetId())]);
            filtered_Ups_mu_1_pt->push_back(mu_pT[Ups_mu_1_Idx->at(cand->GetId())]);
            filtered_Ups_mu_1_isPatLooseMuon->push_back(muIsPatLooseMuon->at(Ups_mu_1_Idx->at(cand->GetId())));
            filtered_Ups_mu_1_isPatSoftMuon->push_back(muIsPatSoftMuon->at(Ups_mu_1_Idx->at(cand->GetId())));
            filtered_Ups_mu_1_isPatMediumMuon->push_back(muIsPatMediumMuon->at(Ups_mu_1_Idx->at(cand->GetId())));
            filtered_Ups_mu_1_isPatTightMuon->push_back(muIsPatTightMuon->at(Ups_mu_1_Idx->at(cand->GetId())));

            filtered_Ups_mu_2_px->push_back(muPx->at(Ups_mu_2_Idx->at(cand->GetId())));
            filtered_Ups_mu_2_py->push_back(muPy->at(Ups_mu_2_Idx->at(cand->GetId())));
            filtered_Ups_mu_2_pz->push_back(muPz->at(Ups_mu_2_Idx->at(cand->GetId())));
            filtered_Ups_mu_2_eta->push_back(mu_eta[Ups_mu_2_Idx->at(cand->GetId())]);
            filtered_Ups_mu_2_phi->push_back(mu_phi[Ups_mu_2_Idx->at(cand->GetId())]);
            filtered_Ups_mu_2_pt->push_back(mu_pT[Ups_mu_2_Idx->at(cand->GetId())]);
            filtered_Ups_mu_2_isPatLooseMuon->push_back(muIsPatLooseMuon->at(Ups_mu_2_Idx->at(cand->GetId())));
            filtered_Ups_mu_2_isPatSoftMuon->push_back(muIsPatSoftMuon->at(Ups_mu_2_Idx->at(cand->GetId())));
            filtered_Ups_mu_2_isPatMediumMuon->push_back(muIsPatMediumMuon->at(Ups_mu_2_Idx->at(cand->GetId())));
            filtered_Ups_mu_2_isPatTightMuon->push_back(muIsPatTightMuon->at(Ups_mu_2_Idx->at(cand->GetId())));
        }
        // Fill into the output tree.
        if(SelectedCands_cut.size() > 0){
            outputTree->Fill();
            // Clear the vectors.
            filtered_Jpsi_1_mass->clear();
            filtered_Jpsi_1_massErr->clear();
            filtered_Jpsi_1_massDiff->clear();
            filtered_Jpsi_1_ctau->clear();
            filtered_Jpsi_1_ctauErr->clear();
            filtered_Jpsi_1_Chi2->clear();
            filtered_Jpsi_1_ndof->clear();
            filtered_Jpsi_1_VtxProb->clear();
            filtered_Jpsi_1_px->clear();
            filtered_Jpsi_1_py->clear();
            filtered_Jpsi_1_pz->clear();
            filtered_Jpsi_1_phi->clear();
            filtered_Jpsi_1_eta->clear();
            filtered_Jpsi_1_pt->clear();
            filtered_Jpsi_1_mu_1_Idx->clear();
            filtered_Jpsi_1_mu_2_Idx->clear();

            filtered_Jpsi_2_mass->clear();
            filtered_Jpsi_2_massErr->clear();
            filtered_Jpsi_2_massDiff->clear();
            filtered_Jpsi_2_ctau->clear();
            filtered_Jpsi_2_ctauErr->clear();
            filtered_Jpsi_2_Chi2->clear();
            filtered_Jpsi_2_ndof->clear();
            filtered_Jpsi_2_VtxProb->clear();
            filtered_Jpsi_2_px->clear();
            filtered_Jpsi_2_py->clear();
            filtered_Jpsi_2_pz->clear();
            filtered_Jpsi_2_phi->clear();
            filtered_Jpsi_2_eta->clear();
            filtered_Jpsi_2_pt->clear();
            filtered_Jpsi_2_mu_1_Idx->clear();
            filtered_Jpsi_2_mu_2_Idx->clear();

            filtered_Pri_mass->clear();
            filtered_Pri_massErr->clear();
            filtered_Pri_ctau->clear();
            filtered_Pri_ctauErr->clear();
            filtered_Pri_Chi2->clear();
            filtered_Pri_ndof->clear();
            filtered_Pri_VtxProb->clear();
            filtered_Pri_px->clear();
            filtered_Pri_py->clear();
            filtered_Pri_pz->clear();
            filtered_Pri_phi->clear();
            filtered_Pri_eta->clear();
            filtered_Pri_pt->clear();

            filtered_Ups_mass->clear();
            filtered_Ups_massErr->clear();
            filtered_Ups_massDiff->clear();
            filtered_Ups_Chi2->clear();
            filtered_Ups_ndof->clear();
            filtered_Ups_VtxProb->clear();
            filtered_Ups_px->clear();
            filtered_Ups_py->clear();
            filtered_Ups_pz->clear();
            filtered_Ups_phi->clear();
            filtered_Ups_eta->clear();
            filtered_Ups_pt->clear();
            filtered_Ups_mu_1_Idx->clear();
            filtered_Ups_mu_2_Idx->clear();

            // For the muons
            filtered_Jpsi_1_mu_1_px->clear();
            filtered_Jpsi_1_mu_1_py->clear();
            filtered_Jpsi_1_mu_1_pz->clear();
            filtered_Jpsi_1_mu_1_eta->clear();
            filtered_Jpsi_1_mu_1_phi->clear();
            filtered_Jpsi_1_mu_1_pt->clear();
            filtered_Jpsi_1_mu_1_isPatLooseMuon->clear();
            filtered_Jpsi_1_mu_1_isPatSoftMuon->clear();
            filtered_Jpsi_1_mu_1_isPatMediumMuon->clear();
            filtered_Jpsi_1_mu_1_isPatTightMuon->clear();

            filtered_Jpsi_1_mu_2_px->clear();
            filtered_Jpsi_1_mu_2_py->clear();
            filtered_Jpsi_1_mu_2_pz->clear();
            filtered_Jpsi_1_mu_2_eta->clear();
            filtered_Jpsi_1_mu_2_phi->clear();
            filtered_Jpsi_1_mu_2_pt->clear();
            filtered_Jpsi_1_mu_2_isPatLooseMuon->clear();
            filtered_Jpsi_1_mu_2_isPatSoftMuon->clear();
            filtered_Jpsi_1_mu_2_isPatMediumMuon->clear();
            filtered_Jpsi_1_mu_2_isPatTightMuon->clear();

            filtered_Jpsi_2_mu_1_px->clear();
            filtered_Jpsi_2_mu_1_py->clear();
            filtered_Jpsi_2_mu_1_pz->clear();
            filtered_Jpsi_2_mu_1_eta->clear();
            filtered_Jpsi_2_mu_1_phi->clear();
            filtered_Jpsi_2_mu_1_pt->clear();
            filtered_Jpsi_2_mu_1_isPatLooseMuon->clear();
            filtered_Jpsi_2_mu_1_isPatSoftMuon->clear();
            filtered_Jpsi_2_mu_1_isPatMediumMuon->clear();
            filtered_Jpsi_2_mu_1_isPatTightMuon->clear();

            filtered_Jpsi_2_mu_2_px->clear();
            filtered_Jpsi_2_mu_2_py->clear();
            filtered_Jpsi_2_mu_2_pz->clear();
            filtered_Jpsi_2_mu_2_eta->clear();
            filtered_Jpsi_2_mu_2_phi->clear();
            filtered_Jpsi_2_mu_2_pt->clear();
            filtered_Jpsi_2_mu_2_isPatLooseMuon->clear();
            filtered_Jpsi_2_mu_2_isPatSoftMuon->clear();
            filtered_Jpsi_2_mu_2_isPatMediumMuon->clear();
            filtered_Jpsi_2_mu_2_isPatTightMuon->clear();

            filtered_Ups_mu_1_px->clear();
            filtered_Ups_mu_1_py->clear();
            filtered_Ups_mu_1_pz->clear();
            filtered_Ups_mu_1_eta->clear();
            filtered_Ups_mu_1_phi->clear();
            filtered_Ups_mu_1_pt->clear();
            filtered_Ups_mu_1_isPatLooseMuon->clear();
            filtered_Ups_mu_1_isPatSoftMuon->clear();
            filtered_Ups_mu_1_isPatMediumMuon->clear();
            filtered_Ups_mu_1_isPatTightMuon->clear();

            filtered_Ups_mu_2_px->clear();
            filtered_Ups_mu_2_py->clear();
            filtered_Ups_mu_2_pz->clear();
            filtered_Ups_mu_2_eta->clear();
            filtered_Ups_mu_2_phi->clear();
            filtered_Ups_mu_2_pt->clear();
            filtered_Ups_mu_2_isPatLooseMuon->clear();
            filtered_Ups_mu_2_isPatSoftMuon->clear();
            filtered_Ups_mu_2_isPatMediumMuon->clear();
            filtered_Ups_mu_2_isPatTightMuon->clear();
        }
        // Marker
        #ifdef SHOW_DEBUG
        puts(">>>>> End of event <<<<<");
        #endif
    }
    printf("Number of entries passing the cut: %d\n", Jpsi1_mass_cut_data.numEntries());
    // Display the histograms.
    // Mass histograms passing the cut are drawn using Roofit.
    #ifndef NO_PLOT
    TCanvas* c2 = new TCanvas("c2", "c2", 1600, 1200);
    c2->Divide(2,2);
    RooPlot* frame5 = Jpsi1_mass_cut_var.frame();
    RooPlot* frame6 = Jpsi2_mass_cut_var.frame();
    RooPlot* frame7 = Ups_mass_cut_var.frame();
    RooPlot* frame8 = Pri_mass_cut_var.frame();
    Jpsi1_mass_cut_data.plotOn(frame5);
    Jpsi2_mass_cut_data.plotOn(frame6);
    Ups_mass_cut_data.plotOn(frame7);
    Pri_mass_cut_data.plotOn(frame8);
    // Draw the histograms.
    c2->cd(1); frame5->Draw();
    c2->cd(2); frame6->Draw();
    c2->cd(3); frame7->Draw();
    c2->cd(4); frame8->Draw();
    // Save png file.
    c2->SaveAs("mass_cut.png");
    #endif

    // Storing the filtered candidates.
    // Create a new ROOT file to save the filtered data.
    TFile* outputFile = new TFile("filtered_data_preCut.root", "RECREATE");
    // Write the output file.
    outputTree->Write();
    outputFile->Close();
    
}

