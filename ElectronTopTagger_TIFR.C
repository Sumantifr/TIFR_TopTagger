#include "TSystem.h"
#include "TH1D.h"
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"
#include "TROOT.h"
#include "TNamed.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include "TClonesArray.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector.h"

#include <math.h>
#include <fstream>

#include "/home/suman/Delphes-3.4.1/modules/Delphes.h"   
#include "/home/suman/Delphes-3.4.1/classes/DelphesClasses.h"
#include "/home/suman/Delphes-3.4.1/classes/DelphesFactory.h"

//  "/home/suman/Delphes-3.4.1" <= directory of Delphes installation 
//  CHANGE IT for YOU! 

#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/PseudoJet.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;


double PhiInRange(const double& phi) {
  Float_t phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}


double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}

double eta_to_theta(double eta){
return(2*atan(exp(-1.*eta)));
}

double diff_func(double f1, double f2){
double ff = pow(f1-f2,2)*1./(pow(f1+f2,2));
return ff;
}

void ElectronTopTagger_TIFR()
{

const int nconstmax = 200;

Float_t pt_cut = 500;	// pt cut on input jet
Float_t jetrad = 0.8;   // jet radius

Float_t beta = 0;	// beta parameter of soft drop
Float_t z_cut = 0.1; ;	// z_cut parameter of soft drop

const int nkT_subjet_max = 2;  // no. of kT subjets


Float_t track_ptcut = 1.;

// PDG id //

Int_t lepPID = 11;//11;
Int_t nuPID = 12;//12;
Int_t bPID = 5;

// Particle masses //

float mass_e = 0.00051;
float mass_mu = 0.1057;
float w_mass = 80.4;
float mu_mass = 0.501;
float R_small = 0.2;
float t_pdg = 172.44;
float b_pdg = 4.2;

bool Muon_finder = false;	// if use muon track directly

const int njetmx = 2;	// how many jets you want to consider

int npfjet;
Float_t jetpt[njetmx];
Float_t jety[njetmx];
Float_t jetphi[njetmx];
Float_t jetener[njetmx];
Float_t jetmass[njetmx];

Float_t jetnhcalbye[njetmx];
Float_t jetneuembynhad[njetmx];
Float_t jetsdmass[njetmx];
Float_t jetchrad[njetmx];
Float_t subhaddiff[njetmx];
Float_t jettau2bytau1[njetmx];

Float_t subjet1_emfrac[njetmx];
Float_t subjet1_hadfrac[njetmx];
Float_t subjet2_emfrac[njetmx];
Float_t subjet2_hadfrac[njetmx];

Float_t lepener[njetmx];
Float_t R_new[njetmx];

Float_t lep_pt[njetmx];
Float_t lep_E[njetmx];
Float_t lep_eta[njetmx];
Float_t lep_phi[njetmx];

Float_t bjet_pt[njetmx];
Float_t bjet_E[njetmx];
Float_t bjet_eta[njetmx];
Float_t bjet_phi[njetmx];

Float_t Zb[njetmx];
Float_t Thetabl[njetmx];

const int nlepmax = 20;

TFile *fileout = new TFile("NTuple_ElectronicTopTagger.root","RECREATE");

TTree *T1 = new TTree("leptop","leptop");

T1->Branch("npfjet",&npfjet,"npfjet/I");
T1->Branch("jetpt",jetpt,"jetpt[npfjet]/F");
T1->Branch("jety",jety,"jety[npfjet]/F");
T1->Branch("jetphi",jetphi,"jetphi[npfjet]/F");
T1->Branch("jetener",jetener,"jetener[npfjet]/F");
T1->Branch("jetmass",jetmass, "jetmass[npfjet]/F");

T1->Branch("jetnhcalbye",jetnhcalbye,"jetnhcalbye[npfjet]/F");
T1->Branch("jetneuembynhad",jetneuembynhad,"jetneuembynhad[npfjet]/F");
T1->Branch("jetsdmass",jetsdmass,"jetsdmass[npfjet]/F");
T1->Branch("jetchrad",jetchrad,"jetchrad[npfjet]/F");
T1->Branch("subhaddiff",subhaddiff,"subhaddiff[npfjet]/F");
T1->Branch("jettau2bytau1",jettau2bytau1,"jettau2bytau1[npfjet]/F");

T1->Branch("lep_pt",lep_pt,"lep_pt[npfjet]/F");
T1->Branch("lep_E",lep_E,"lep_E[npfjet]/F");
T1->Branch("lep_eta",lep_eta,"lep_eta[npfjet]/F");
T1->Branch("lep_phi",lep_phi,"lep_phi[npfjet]/F");

T1->Branch("bjet_pt",bjet_pt,"bjet_pt[npfjet]/F");
T1->Branch("bjet_E",bjet_E,"bjet_E[npfjet]/F");
T1->Branch("bjet_eta",bjet_eta,"bjet_eta[npfjet]/F");
T1->Branch("bjet_phi",bjet_phi,"bjet_phi[npfjet]/F");

T1->Branch("Zb",Zb,"Zb[npfjet]/F");
T1->Branch("Thetabl",Thetabl,"Thetabl[npfjet]/F");

fastjet::Strategy               strategy = fastjet::Best;
fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;

fastjet::JetDefinition* jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, jetrad,recombScheme, strategy);
fastjet::JetDefinition* jetDefCA = new fastjet::JetDefinition(fastjet::cambridge_algorithm, jetrad,recombScheme, strategy);
fastjet::JetDefinition* jetDefkT = new fastjet::JetDefinition(fastjet::kt_algorithm, 10,recombScheme, strategy);
SoftDrop rsd(beta, z_cut);

vector<PseudoJet> inputList;

PseudoJet psjet;

TChain *chain = new TChain("Delphes");

// Write their names in a .log file as a list and use that file as an input while running the executable
//(it will ask for it in the beginning)

char datafile[100];
char rootfiles[100];

cout<<"Which files you wanna use? \n";
cin>>rootfiles;

TFileCollection fc("fileCol","",rootfiles);
chain->AddFileInfoList(fc.GetList());
std::cout << "Files found : " << fc.GetNFiles()<< std::endl;
fc.Print();


ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
ExRootResult *result = new ExRootResult();

Long64_t allEntries = treeReader->GetEntries();
cout<<"allEntries "<<allEntries<<endl;

Track *tracks;

TClonesArray *branch_pftrk = treeReader->UseBranch("EFlowTrack");
TClonesArray *branch_pfpho = treeReader->UseBranch("EFlowPhoton");
TClonesArray *branch_pfneuhad = treeReader->UseBranch("EFlowNeutralHadron");

Track *pftrack;
Tower *pfphoton;
Tower *pfneuhad;

TLorentzVector jet4;
TLorentzVector sdjet4;
TLorentzVector b_vec4;
TLorentzVector lep4;
TLorentzVector new_nu_vec4;
bool lep_finder = false;

for(unsigned int entry=0; entry < allEntries; ++entry){

if(entry%100==0) {cout<<"entry "<<entry<<endl ; }

 treeReader->ReadEntry(entry);

  npfjet = 0;

  inputList.clear();

	if((branch_pftrk->GetEntriesFast() > 0) || (branch_pfpho->GetEntriesFast() > 0) || (branch_pfneuhad->GetEntriesFast() > 0)){
		
		 jet4.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
         
         for(int ij=0; ij< branch_pftrk->GetEntriesFast(); ij++){
                pftrack = (Track*) branch_pftrk->At(ij); 
                psjet = PseudoJet((pftrack->P4()).Px(),(pftrack->P4()).Py(),(pftrack->P4()).Pz(),(pftrack->P4()).E());
                psjet.set_user_index(ij);
                inputList.push_back(psjet);
             }

		for(int ij=0; ij< branch_pfpho->GetEntriesFast(); ij++){
                pfphoton = (Tower*) branch_pfpho->At(ij); 
                psjet = PseudoJet((pfphoton->P4()).Px(),(pfphoton->P4()).Py(),(pfphoton->P4()).Pz(),(pfphoton->P4()).E());
                psjet.set_user_index(ij+(branch_pftrk->GetEntriesFast()));
                inputList.push_back(psjet);
            }
		
		for(int ij=0; ij< branch_pfneuhad->GetEntriesFast(); ij++){
                pfneuhad = (Tower*) branch_pfneuhad->At(ij); 
                psjet = PseudoJet((pfneuhad->P4()).Px(),(pfneuhad->P4()).Py(),(pfneuhad->P4()).Pz(),(pfneuhad->P4()).E());
                psjet.set_user_index(ij+(branch_pftrk->GetEntriesFast())+(branch_pfpho->GetEntriesFast()));
                inputList.push_back(psjet); 
              }

   vector <PseudoJet> sortedJets, CA_sortedJets;
   ClusterSequence clustSeq(inputList, *jetDef);

   sortedJets    = sorted_by_pt(clustSeq.inclusive_jets());
   if(sortedJets.size()<1) {continue;}
   
   int ijet = 0;
   if(sortedJets[ijet].pt()<pt_cut) { continue; }				// apply pt cut
   if(fabs(sortedJets[ijet].eta())>2.5)  { continue; }		// apply \eta cut
  
   vector <fastjet::PseudoJet> const_j = sortedJets[ijet].constituents(); 
 
   // SD//
   
   ClusterSequence clustSeq_CA(const_j, *jetDefCA);
   CA_sortedJets = sorted_by_pt(clustSeq_CA.inclusive_jets()); 
   fastjet::PseudoJet rsd_jet = rsd(CA_sortedJets[0]);
  
   sortedJets[ijet] = rsd_jet;
   
   //now jet is groomed 
   
   vector <fastjet::PseudoJet> const_0 = sortedJets[ijet].constituents(); 
   
   jet4.SetPxPyPzE(sortedJets[ijet].px(),sortedJets[ijet].py(),sortedJets[ijet].pz(),sortedJets[ijet].e());
  
   // basic kinematic variables from jet
  
   jetpt[ijet] = sortedJets[ijet].pt();
   jety[ijet] = sortedJets[ijet].rapidity();
   jetphi[ijet] = sortedJets[ijet].phi();
   jetener[ijet] = sortedJets[ijet].e();
   jetmass[ijet] = sortedJets[ijet].m(); 
   
   // calculating different energy components, and shape variables
   
   float em_e = 0, had_e = 0, pho_e = 0, neuhad_e = 0;
   int leadleptrack = -1;
   float leadleptrkpt = -100;
   
   if(const_0.size() > 0){
	
	float sumpt = 0;
	
	jetchrad[ijet] = 0;
	 
    for(unsigned icons=0 ; icons <const_0.size(); icons++){
					
		int id = const_0[icons].user_index();
		
		if(id < branch_pftrk->GetEntriesFast()){
			
			tracks = (Track*) branch_pftrk->At(id);
			
			if(((tracks->P4()).Pt() > leadleptrkpt) && (abs(tracks->PID) == lepPID)){
				leadleptrkpt = (tracks->P4()).Pt();
				leadleptrack = id;
				}
			
			if(abs(tracks->PID) == 11) { em_e += (tracks->P4()).E(); }										// EM energy of jet
			if((abs(tracks->PID) != 11) && (abs(tracks->PID) != 13)) { had_e += (tracks->P4()).E(); }		// Hadronic energy of jet
			
			jetchrad[ijet] += ((tracks->Charge) * (tracks->P4()).Pt() * delta2R((tracks->P4()).Eta(),(tracks->P4()).Phi(),jety[ijet],jetphi[ijet])); // charge radius
		
		}
			
		if((id >= branch_pftrk->GetEntriesFast()) && (id < (branch_pftrk->GetEntriesFast() + branch_pfpho->GetEntriesFast()))){
			
			id -= branch_pftrk->GetEntriesFast();
			
			em_e += (((Tower*) branch_pfpho->At(id))->P4()).E();  // EM energy of jet
			pho_e += (((Tower*) branch_pfpho->At(id))->P4()).E(); // Photon energy of jet
			
			}
			
	   	if((id >= (branch_pftrk->GetEntriesFast()+branch_pfpho->GetEntriesFast())) && (id < (branch_pftrk->GetEntriesFast() + branch_pfpho->GetEntriesFast() + branch_pfneuhad->GetEntriesFast()))){
			
			id -= (branch_pftrk->GetEntriesFast()+branch_pfpho->GetEntriesFast());
			
			had_e += (((Tower*) branch_pfneuhad->At(id))->P4()).E(); // Hadronic energy of jet
			neuhad_e += (((Tower*) branch_pfneuhad->At(id))->P4()).E(); // Neutral hadronic energy of jet
			
			 }
			
			sumpt += sortedJets[ijet].constituents()[icons].perp();
			
			} // loop over constituents
			
			jetnhcalbye[ijet] = (1. - had_e*1./jetener[ijet]); 
			jetneuembynhad[ijet] = (jetnhcalbye[ijet] > 1.e-9)? (pho_e*1./jetener[ijet])*1./jetnhcalbye[ijet] : -1000;
			
			jetchrad[ijet] *= 1./sumpt;
			
			// Subjettiness //
			
			Nsubjettiness nsub1_beta1(1,OnePass_KT_Axes(), UnnormalizedMeasure(1.));
			float tau1 = nsub1_beta1(sortedJets[ijet]); 
			Nsubjettiness nsub2_beta1(2,OnePass_KT_Axes(), UnnormalizedMeasure(1.));
			float tau2 = nsub2_beta1(sortedJets[ijet]);
			
			if(tau1 > 1.e-9) { jettau2bytau1[ijet] =  (tau2*1./tau1)  ; }
			
			// Subjettiness ends //
			
			sdjet4.SetPxPyPzE(rsd_jet.px(),rsd_jet.py(),rsd_jet.pz(),rsd_jet.e());
			
			jetsdmass[ijet] = rsd_jet.m();
			
			// subjet energy fractions
			
			Float_t em_sub[2]={0,0};
			Float_t had_sub[2]={0,0};
			
			vector<PseudoJet> pieces;
			pieces = rsd_jet.pieces();
			
			// variables from subjets
			
			if(pieces.size() > 0){
				
				TLorentzVector tmp0, tmp1;
				tmp0.SetPxPyPzE(pieces[0].px(),pieces[0].py(),pieces[0].pz(),pieces[0].e());
				tmp1.SetPxPyPzE(pieces[1].px(),pieces[1].py(),pieces[1].pz(),pieces[1].e());
			
				for(unsigned ip=0; ip<pieces.size(); ip++){
				  for(unsigned icons=0 ; icons < pieces[ip].constituents().size(); icons++){
				
					int id = (pieces[ip].constituents())[icons].user_index();
					
					if(id < branch_pftrk->GetEntriesFast()){
			
						tracks = (Track*) branch_pftrk->At(id);
			
						if(abs(tracks->PID) == 11) { em_sub[ip] += (tracks->P4()).E(); }
						if((abs(tracks->PID) != 11) && (abs(tracks->PID) != 13)) { had_sub[ip] += (tracks->P4()).E(); }
			
					}
			
					if((id >= branch_pftrk->GetEntriesFast()) && (id < (branch_pftrk->GetEntriesFast() + branch_pfpho->GetEntriesFast()))){
			
						id -= branch_pftrk->GetEntriesFast();
			
						em_sub[ip] += (((Tower*) branch_pfpho->At(id))->P4()).E();
						
					}
			
					if((id >= (branch_pftrk->GetEntriesFast()+branch_pfpho->GetEntriesFast())) && (id < (branch_pftrk->GetEntriesFast() + branch_pfpho->GetEntriesFast() + branch_pfneuhad->GetEntriesFast()))){
			
						id -= (branch_pftrk->GetEntriesFast()+branch_pfpho->GetEntriesFast());
			
						had_sub[ip] += (((Tower*) branch_pfneuhad->At(id))->P4()).E();
						
						}
					
					} // loop over subjet constituents
					
					  
					   em_sub[ip] *= 1./pieces[ip].e();
					   had_sub[ip] *= 1./pieces[ip].e();
					
				} // loop over subjets
				
					   subhaddiff[ijet] = diff_func(had_sub[0],had_sub[1]);
					 
					   subjet1_hadfrac[ijet] = had_sub[0];	
					   subjet1_emfrac[ijet] = em_sub[0];
					   
					   subjet2_hadfrac[ijet] = had_sub[1];	
					   subjet2_emfrac[ijet] = em_sub[1];
					 
					
			} // if any subjets
			
	
	 if(pieces.size() > 0){
		
		 // Lepton finder from subjets starts//
		
		 int lepsubtag = -1;
		 int lepmatch = -1;
		 vector <fastjet::PseudoJet> kTl_sortedJets;
		 
		 lep4.SetPxPyPzE(0,0,0,0);
		
		 if(subjet1_hadfrac[ijet]<subjet2_hadfrac[ijet]){ lepsubtag = 0; }
			else{ lepsubtag = 1;  }
			
		
		 int leadtrackid = -1;
		 
		 if(lepsubtag >= 0){
			 
			vector <fastjet::PseudoJet> const_lsub = pieces[lepsubtag].constituents();     
			
			int ntracksub = 0;
			for(int ic=0; ic<const_lsub.size(); ic++){
				int id = const_lsub[ic].user_index();
				if(id < branch_pftrk->GetEntriesFast()){
					ntracksub++;
				}
			}
		
			if(ntracksub == 0) { 
				lepsubtag = 1 - lepsubtag ; 
				const_lsub = pieces[lepsubtag].constituents();  
				}	
		 
					
			ClusterSequence clustSeq_kT(const_lsub, *jetDefkT);
			kTl_sortedJets = sorted_by_pt(clustSeq_kT.inclusive_jets());
			
			// finding highest pT track within EM enriched subjet
			
			float ptmax = -0.1; 
			for(int ic=0; ic<const_lsub.size(); ic++){
				int id = const_lsub[ic].user_index();
				if(id < branch_pftrk->GetEntriesFast()){
					tracks = (Track*) branch_pftrk->At(id);
					if(tracks->PT > ptmax){
						leadtrackid = id;
						ptmax = tracks->PT;
					}
				}
			}
		
		    if(kTl_sortedJets.size() > 0 && leadtrackid>=0){
				
				int ilms = -1;
				vector <fastjet::PseudoJet> lkT_subjets =  sorted_by_pt(kTl_sortedJets[0].exclusive_subjets_up_to(nkT_subjet_max));
				
				for(int kk=0; kk<lkT_subjets.size(); kk++){
					vector <fastjet::PseudoJet> const_kTsub = lkT_subjets[kk].constituents();     
						for(int icons=0; icons<const_kTsub.size(); icons++){
							int id = const_kTsub[icons].user_index();
							if(id < branch_pftrk->GetEntriesFast() && id==leadtrackid){
								lepmatch = 1;
								ilms = kk;
								break;
								}
							}
						}
					
				
					if(ilms>=0){
						lep4.SetPxPyPzE(lkT_subjets[ilms].px(),lkT_subjets[ilms].py(),lkT_subjets[ilms].pz(),lkT_subjets[ilms].e());
						  }
				
					lkT_subjets.clear();
				
					tracks = (Track*) branch_pftrk->At(leadtrackid);
					lep4.SetTheta(eta_to_theta((tracks->P4()).Eta()));
					lep4.SetPhi((tracks->P4()).Phi());
					
				}
		 
				const_lsub.clear();
				kTl_sortedJets.clear();	
			
			}
		
		 // Lepton Finder from subjets ends //		
		 
	if(Muon_finder){
		if(leadleptrack>=0){
			lepmatch = 1;
			tracks = (Track*) branch_pftrk->At(leadleptrack);
			lep4.SetPtEtaPhiE((tracks->P4()).Pt(),(tracks->P4()).Eta(),(tracks->P4()).Phi(),(tracks->P4()).E());
		}
	}	
		
	if(lepmatch>=0){
		lep_finder = true;
		b_vec4 = sdjet4 - lep4;
		b_vec4.SetPxPyPzE(pieces[1-lepsubtag].px(),pieces[1-lepsubtag].py(),pieces[1-lepsubtag].pz(),pieces[1-lepsubtag].e());
		lep_pt[ijet] = lep4.Pt();
		lep_eta[ijet] = lep4.Eta();
		lep_phi[ijet] = lep4.Phi();
		lep_E[ijet] = lep4.E();
		
		bjet_pt[ijet] = b_vec4.Pt();
		bjet_eta[ijet] = b_vec4.Eta();
		bjet_phi[ijet] = b_vec4.Phi();
		bjet_E[ijet] = b_vec4.E();
		
		lepener[ijet] = lep4.E();
		
		float nx = lep4.Px()/lep4.P();
		float ny = lep4.Py()/lep4.P();
		float nz = lep4.Pz()/lep4.P();

		TLorentzVector new_nu_vec4_perp;
			
		float bx = b_vec4.Px();
		float by = b_vec4.Py();
		float bz = b_vec4.Pz();
		
		b_pdg = b_vec4.M();
	   
		Thetabl[ijet] = (lep4.P()*1./b_vec4.P())*(1./(w_mass*w_mass - mass_mu*mass_mu))*(t_pdg*t_pdg - w_mass*w_mass - b_pdg*b_pdg - 2*b_vec4.Dot(lep4));
		
		float pll_dot = b_vec4.Vect() * lep4.Vect().Unit();
		
		R_new[ijet] = (2 * w_mass * w_mass  * 1./lep4.E()) * (b_vec4.E() - pll_dot) * 1./(t_pdg*t_pdg - w_mass*w_mass - b_pdg*b_pdg - 2*b_vec4.Dot(lep4));
		R_new[ijet] = sqrt(R_new[ijet]);

		float costheB = (3*b_vec4.E()-pll_dot)*R_new[ijet]/(4*b_vec4.P());

		float nupll = 0.5 * (w_mass*w_mass)*1./(lep4.E()*(sqrt(1+R_new[ijet]*R_new[ijet]) )- lep4.P());
		float nue = nupll * sqrt(1+R_new[ijet]*R_new[ijet]);
        
        Zb[ijet] = (b_vec4.E())*1./(nue+lep4.E());

			} //lepmatch
		} // using lepton from sd subjets
		 //lep_finder
			
	} // if any constituent in jet

    sortedJets.clear();
    CA_sortedJets.clear();
    
    if(isnan(jettau2bytau1[ijet])||(!lep_finder)) { jettau2bytau1[ijet] = -100; }
    if(isnan(jetnhcalbye[ijet])||(!lep_finder)) { jetnhcalbye[ijet] = -100; }
    if(isnan(jetneuembynhad[ijet])||(!lep_finder)) { jetneuembynhad[ijet] = -100; }
    if(isnan(jetchrad[ijet])||(!lep_finder)) { jetchrad[ijet] = -100; }
    if(isnan(subhaddiff[ijet])||(!lep_finder)) { subhaddiff[ijet] = -100; }
    if(isnan(jetsdmass[ijet])||(!lep_finder)) { jetsdmass[ijet] = -100; }
    
    if(isnan(Zb[ijet])||(!lep_finder)) { Zb[ijet] = -100; }
    if(isnan(Thetabl[ijet])||(!lep_finder)) {	Thetabl[ijet] = -100; }
    
    npfjet++;
    
    if(npfjet>njetmx) {continue;} 
  
      }// if any track, photon, neutral hadron
      
   T1->Fill();
  
    }	


fileout->cd();
fileout->Write();
fileout->Close();

}

int main()
{
ElectronTopTagger_TIFR() ;
}
