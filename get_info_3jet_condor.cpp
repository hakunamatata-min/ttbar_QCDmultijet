// get useful information from NanoAOD files for ttbar analysis
// each event contains the information for top quark pairs and their decay
// products written by Renqi Pan in 16th June, 2021.
//////////////////////////////////////////////////////////////////////////////
// declare global vraibles
Float_t nu_px, nu_py, neutrino_pz;
TLorentzVector mom_nu, mom_lep, mom_jets[45];
Int_t btag_num; // btag_num: number of bjets among the leading four jets
Float_t LepCharge, Jet_btagDeepB[45],jet_btagDeepB[45]; // charge of the satisfied leading lepton
UInt_t jet_num;     // count the number of jets satisfy pt and eta criteria
int bjets_index[2]; // store the indexes of two bjets with the hightest btaged
                  // score.
int jets_index[45]; // store the indexes of light jets
int jet_index[45];// index for selected jets
int bjet_lep, bjet_had, min_j1, min_j2; // denote the minimum likelihood case
Float_t mass_wlep, mass_whad, mass_tlep, mass_thad;
Double_t pro_wlep, pro_tlep, pro_thad, pro_whad, pro_twlep;
int bindex; // for loop over b-jet index

float xi_thad = 18.0, x0_thad = 179, xi_wlep = 2.0, x0_wlep = 80, xi_tlep = 8.5,
      x0_tlep = 169, xi_whad = 14.0, x0_whad = 84;
float mw_lep = 80, mt_had = 173, mt_lep = 173, mw_had = 80.0,
	  sigmat_had = 30, sigmaw_lep = 20.0,sigmat_lep = 20, sigmaw_had = 24.0;
TLorentzVector mom_top, mom_antitop;
float rectop_mass, recantitop_mass, rectop_pt, mass_tt, rapidity_tt,
      rectop_eta,rectop_rapidity,rectop_costheta;
Double_t nupz_min = -1000.0, nupz_max = 1000.0;
Double_t minimum;//minimum likelihood
float rapidity_bb,mass_bbjjl,deltaR_bb,rapidity_bl;

////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//define likelihood function with two b-jets
Double_t likelihood(Double_t *pz,Double_t *pars){
    Double_t nu_pz=pz[0];
    int j1=pars[0];
    int j2=pars[1];
    Float_t nu_E=sqrt(nu_px*nu_px+nu_py*nu_py+nu_pz*nu_pz);
    mom_nu= TLorentzVector(nu_px,nu_py,nu_pz,nu_E);
    mass_wlep=(mom_nu+mom_lep).M();
    mass_whad=(mom_jets[jets_index[j1]]+mom_jets[jets_index[j2]]).M();
    mass_tlep=(mom_nu+mom_lep+mom_jets[bjets_index[bindex]]).M();
    if(bindex==0){
          mass_thad=(mom_jets[jets_index[j1]]+mom_jets[jets_index[j2]]+mom_jets[bjets_index[1]]).M();
          }
    else
          mass_thad=(mom_jets[jets_index[j1]]+mom_jets[jets_index[j2]]+mom_jets[bjets_index[0]]).M();
    
    pro_wlep=ROOT::Math::gaussian_pdf(mass_wlep,sigmaw_lep,mw_lep);
    pro_tlep=ROOT::Math::gaussian_pdf(mass_tlep,sigmat_lep,mt_lep);
    pro_whad=ROOT::Math::gaussian_pdf(mass_whad,sigmaw_had,mw_had);
    pro_thad=ROOT::Math::gaussian_pdf(mass_thad,sigmat_had,mt_had);
    Double_t log_nupz;
    log_nupz=-TMath::Log(pro_tlep)-TMath::Log(pro_wlep)-TMath::Log(pro_thad)-TMath::Log(pro_whad); 
     
    return log_nupz;
}
/////////////////////////////////////////////////////////////
//define likelihood function with two b-jets
Double_t likelihood_3jet(Double_t *pz,Double_t *pars){
    sigmaw_lep=10; sigmat_had=35; sigmat_lep=20;
    Double_t nu_pz=pz[0];
    Float_t nu_E=sqrt(nu_px*nu_px+nu_py*nu_py+nu_pz*nu_pz);
    mom_nu= TLorentzVector(nu_px,nu_py,nu_pz,nu_E);
    mass_wlep=(mom_nu+mom_lep).M();
    mass_tlep=(mom_nu+mom_lep+mom_jets[bjets_index[bindex]]).M();
    if(bindex==0){
          mass_thad=(mom_jets[jets_index[0]]+mom_jets[bjets_index[1]]).M();
          }
    else
          mass_thad=(mom_jets[jets_index[0]]+mom_jets[bjets_index[0]]).M();
    
    pro_wlep=ROOT::Math::gaussian_pdf(mass_wlep,sigmaw_lep,mw_lep);
    pro_tlep=ROOT::Math::gaussian_pdf(mass_tlep,sigmat_lep,mt_lep);
    pro_thad=ROOT::Math::gaussian_pdf(mass_thad,sigmat_had,mt_had);
    Double_t log_nupz;
    log_nupz=-TMath::Log(pro_tlep)-TMath::Log(pro_wlep)-TMath::Log(pro_thad); 
     
    return log_nupz;
}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// reconstruct ttbar pair using likelihood function
// accorind to the leading four jets and two bjets among them.
void recons_tt() {
    //if (btag_num >= 2) {
    //reorder bjet scores acrroding to their scores
    int index[jet_num];
    for (int i = 0; i < jet_num; i++)
      index[i] = jet_index[i];

    for (int kk = 0; kk < 2; kk++) {
      int max = kk;
      for (int tt = kk + 1; tt < jet_num; tt++) {
        if (Jet_btagDeepB[index[tt]] > Jet_btagDeepB[index[max]]) {
          max = tt;
        }
      }
      int tmp = index[max];
      index[max] = index[kk];
      index[kk] = tmp;
    }
    bjets_index[0] = index[0];
    bjets_index[1] = index[1];
    int light_jets = 0;
    for (int k = 0; k < jet_num; k++) {
      if (jet_index[k] != bjets_index[0] && jet_index[k] != bjets_index[1]) {
        jets_index[light_jets] = jet_index[k]; //get index for light jets
        light_jets++; // light flavour jets
      }
    }

    Double_t minimum_likelihood, nupz;
    bjet_lep = 0, bjet_had = 1, min_j1 = 0, min_j2 = 51;
    //for at least 4 jets  
    if(light_jets >=2){
	    for (bindex = 0; bindex < 2; bindex++) {
	      for (int j1 = 0; j1 < jet_num - 2; j1++) {
	        for (int j2 = j1 + 1; j2 < jet_num - 2; j2++) {

	          TF1 *likelihood_fun =
	              new TF1("likelihood_fun", likelihood, nupz_min, nupz_max, 2);
	          likelihood_fun->SetParameters(j1, j2); // pass the index of j1 and j2 as parameters to a
	          minimum_likelihood = likelihood_fun->GetMinimum(nupz_min, nupz_max);
	          nupz = likelihood_fun->GetMinimumX(nupz_min, nupz_max);
	          if (bindex == 0 && j1 == 0 && j2 == 1) {
	            minimum = minimum_likelihood;
	            neutrino_pz = nupz;
	          } else if (minimum_likelihood < minimum) {
	            minimum = minimum_likelihood;
	            neutrino_pz = nupz;
	            bjet_lep = bindex;
	            bjet_had = bindex == 0 ? 1 : 0;
	            min_j1 = j1;
	            min_j2 = j2;
	          }
	          
	        }
	      }
	    }

	    Float_t nu_E =
	        sqrt(nu_px * nu_px + nu_py * nu_py + neutrino_pz * neutrino_pz);
	    mom_nu = TLorentzVector(nu_px, nu_py, neutrino_pz, nu_E);
	    mass_wlep = (mom_nu + mom_lep).M();
	    mass_whad =
	        (mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]]).M();
	    mass_tlep = (mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]]).M();
	    mass_thad = (mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]] +
	                 mom_jets[bjets_index[bjet_had]]).M();
	    if (LepCharge > 0) {
	      mom_top = mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]];
	      mom_antitop = mom_jets[jets_index[min_j1]] +
	                    mom_jets[jets_index[min_j2]] +
	                    mom_jets[bjets_index[bjet_had]];
	    } else {
	      mom_top = mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]] +
	                mom_jets[bjets_index[bjet_had]];
	      mom_antitop = mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]];
	    }
	}
	//for exactly 3 jets
	else if (light_jets==1){
		for (bindex = 0; bindex < 2; bindex++) {
	          TF1 *likelihood_fun =
	              new TF1("likelihood_fun", likelihood_3jet, nupz_min, nupz_max);
	          minimum_likelihood = likelihood_fun->GetMinimum(nupz_min, nupz_max);
	          nupz = likelihood_fun->GetMinimumX(nupz_min, nupz_max);
	          if (bindex == 0 ) {
	            minimum = minimum_likelihood;
	            neutrino_pz = nupz;
	          } else if (minimum_likelihood < minimum) {
	            minimum = minimum_likelihood;
	            neutrino_pz = nupz;
	            bjet_lep = bindex;
	            bjet_had = bindex == 0 ? 1 : 0;
	            //min_j1 = 0;
	            //min_j2 = 0;
	          }          
	    }
      minimum=minimum+4.0;
	    Float_t nu_E =sqrt(nu_px * nu_px + nu_py * nu_py + neutrino_pz * neutrino_pz);
	    mom_nu = TLorentzVector(nu_px, nu_py, neutrino_pz, nu_E);
	    mass_wlep = (mom_nu + mom_lep).M();
	    mass_whad =(mom_jets[jets_index[0]]).M();
	    mass_tlep = (mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]]).M();
	    mass_thad = (mom_jets[jets_index[0]] + mom_jets[bjets_index[bjet_had]]).M();
	    if (LepCharge > 0) {
	      mom_top = mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]];
	      mom_antitop = mom_jets[jets_index[0]] + mom_jets[bjets_index[bjet_had]];
	    } else {
	      mom_top = mom_jets[jets_index[0]] + mom_jets[bjets_index[bjet_had]];
	      mom_antitop = mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]];
	    }

	}
    rectop_mass = mom_top.M();
    recantitop_mass = mom_antitop.M();
    rectop_pt = mom_top.Pt();
    rectop_rapidity=mom_top.Rapidity();
    rectop_eta=mom_top.Eta();
    rectop_costheta=mom_top.CosTheta();
    rapidity_tt = mom_top.Rapidity() - mom_antitop.Rapidity();
    mass_tt = (mom_antitop + mom_top).M();
  //  deltaR_bb=mom_jets[bjets_index[bjet_lep]].DeltaR(mom_jets[bjets_index[bjet_had]]);
   // rapidity_bl=mom_jets[bjets_index[bjet_lep]].Rapidity()-mom_lep.Rapidity();
   // rapidity_bb=mom_jets[bjets_index[bjet_lep]].Rapidity()-mom_jets[bjets_index[bjet_had]].Rapidity();
   // mass_bbjjl=(mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]] +
    //            mom_jets[bjets_index[bjet_had]]+mom_lep + mom_jets[bjets_index[bjet_lep]]).M();
  //}
}
///////////////////////////////////////////////////////////////////////
// select the semileptonic final states and reconstruct top quark pairs.
void get_info_3jet_condor(TString inputFile,TString dataset) {
  TChain chain("Events");
  TString indir="root://cms-xrd-global.cern.ch/";
  //TString inputFile ="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root";
  //TString indir="./"; 
  chain.Add(indir+inputFile);
  TString output = "new_"+dataset+".root";
  TFile *file = new TFile(output, "RECREATE");
  TTree *mytree = new TTree("mytree", " tree with branches");
  TTree *rawtree = new TTree("rawtree", "tree without selection");
  Int_t nevents = 0, nevents2 = 0,count_4jet=0, count_3jet=0; // count the number of events written in tree

  cout << inputFile << " is reading and processing" << endl;
  cout << "total number of events: " << chain.GetEntries() << endl;
  
  ////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
  Float_t LHEPart_eta[9], LHEPart_mass[9], LHEPart_phi[9], LHEPart_pt[9];
  Int_t LHEPart_pdgId[9], LHEPart_status[9];
  UInt_t nLHEPart;
  UInt_t tt_efficiency;
  UInt_t LHE_nlep = 0, LHE_nhad = 0, LHE_tao = 0;
  Float_t M_tt_gen, delta_rapidity_gen, lep_charge;
  Float_t top_pt, top_eta, top_mass, top_phi, antitop_pt, antitop_eta,
      antitop_phi, antitop_mass;
  Float_t b_pt, b_eta, b_mass, b_phi, antib_pt, antib_eta, antib_phi,
      antib_mass;
  Float_t lep_pt, lep_eta, lep_mass, lep_phi, nu_pt, nu_eta, nu_phi, nu_mass;
  Float_t up_pt, up_eta, up_mass, up_phi, down_pt, down_eta, down_phi, down_mass;
  // Int_t LHE_had[6];
  

  /////////////////////////////////////////////////////////
  // difine branch for final state at detector level
  Float_t MET_pt, MET_phi;
  Float_t Electron_eta[9], Electron_mass[9], Electron_pt[9], Electron_phi[9];
  Float_t Muon_mass[9], Muon_phi[9], Muon_pt[9], Muon_eta[9];
  Float_t lepton_mass[18], lepton_phi[18], lepton_eta[18], lepton_pt[18];
  UInt_t nMuon, nElectron, nJet, nlepton, Jet_btaged[45], nBtag;
  Int_t Jet_partonFlavour[45], Muon_charge[9], Electron_charge[9],
      lepton_charge[18];
  Float_t Jet_btagCSVV2[45], Jet_eta[45], Jet_mass[45], Jet_phi[45], Jet_pt[45];
  Float_t jet_btagCSVV2[45], jet_eta[45], jet_mass[45], jet_phi[45], jet_pt[45],jet_partonFlavour[45],jet_btaged[45];
  Float_t MtW;

  chain.SetBranchAddress("Electron_phi", Electron_phi);
  chain.SetBranchAddress("Electron_pt", Electron_pt);
  chain.SetBranchAddress("Electron_mass", Electron_mass);
  chain.SetBranchAddress("Electron_eta", Electron_eta);
  chain.SetBranchAddress("nElectron", &nElectron);
  chain.SetBranchAddress("Electron_charge", Electron_charge);
  chain.SetBranchAddress("nMuon", &nMuon);
  chain.SetBranchAddress("nJet", &nJet);
  chain.SetBranchAddress("Muon_eta", Muon_eta);
  chain.SetBranchAddress("Muon_pt", Muon_pt);
  chain.SetBranchAddress("Muon_phi", Muon_phi);
  chain.SetBranchAddress("Muon_mass", Muon_mass);
  chain.SetBranchAddress("Muon_charge", Muon_charge);
  chain.SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour);
  chain.SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2);
  chain.SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB);
  chain.SetBranchAddress("Jet_eta", Jet_eta);
  chain.SetBranchAddress("Jet_pt", Jet_pt);
  chain.SetBranchAddress("Jet_phi", Jet_phi);
  chain.SetBranchAddress("Jet_mass", Jet_mass);

  //for physics object and event seletions
  Float_t Muon_pfRelIso04_all[9];
  Int_t Electron_cutBased[9], Jet_jetId[45],PV_npvsGood;
  Bool_t Muon_tightId[9], Muon_looseId[9];
  Float_t Electron_deltaEtaSC[9],Electron_dxy[9], Electron_dz[9];
  chain.SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all);
  chain.SetBranchAddress("Muon_tightId",Muon_tightId);
  chain.SetBranchAddress("Muon_looseId",Muon_looseId);
  chain.SetBranchAddress("Electron_cutBased",Electron_cutBased);
  chain.SetBranchAddress("Jet_jetId",Jet_jetId);
  chain.SetBranchAddress("Electron_deltaEtaSC",Electron_deltaEtaSC);
  chain.SetBranchAddress("Electron_dz",Electron_dz);
  chain.SetBranchAddress("Electron_dxy",Electron_dxy);
  chain.SetBranchAddress("PV_npvsGood",&PV_npvsGood);
//for systematic uncertainties
  Float_t btagWeight_DeepCSVB,Generator_weight;
  chain.SetBranchAddress("Generator_weight",&Generator_weight);
  mytree->Branch("Generator_weight",&Generator_weight,"Generator_weight/F");
  mytree->Branch("jet_num", &jet_num,"jet_num/I"); // number of jets satisfy the  seletion criteria
  rawtree->Branch("nJet", &nJet, "nJet/I");
  rawtree->Branch("Generator_weight",&Generator_weight,"Generator_weight/F");
  ////////////////////////////////////////////////////////////////
  // add information at reconstruction level.
  mytree->Branch("rectop_pt", &rectop_pt, "rectop_pt/F");
  mytree->Branch("rectop_eta",&rectop_eta,"rectop_eta/F");
  mytree->Branch("rectop_rapidity",&rectop_rapidity,"rectop_rapidity/F");
  mytree->Branch("rectop_costheta",&rectop_costheta,"rectop_costheta/F");
  mytree->Branch("rectop_mass", &rectop_mass, "rectop_mass/F");
  mytree->Branch("recantitop_mass", &recantitop_mass, "recantitop_mass/F");
  mytree->Branch("rapidity_tt", &rapidity_tt, "rapidity_tt/F");
  mytree->Branch("mass_tt", &mass_tt, "mass_tt/F");
  mytree->Branch("mass_wlep", &mass_wlep, "mass_wlep/F");
  mytree->Branch("mass_whad", &mass_whad, "mass_whad/F");
  mytree->Branch("neutrino_pz", &neutrino_pz, "neutrino_pz/F");
  mytree->Branch("mass_thad", &mass_thad, "mass_thad/F");
  mytree->Branch("mass_tlep", &mass_tlep, "mass_tlep/F");
  mytree->Branch("MtW",&MtW,"MtW/F");
  //mytree->Branch("rapidity_bb",&rapidity_bb,"rapidity_bb/F");
  //mytree->Branch("deltaR_bb",&deltaR_bb,"deltaR_bb/F");
  //mytree->Branch("mass_bbjjl",&mass_bbjjl,"mass_bbjjl/F");
  //mytree->Branch("rapidity_bl",&rapidity_bl,"rapidity_bl/F");
  mytree->Branch("likelihood",&minimum,"minimum/D" );
  //////////////////////////////////////////////////////////////////////
  // loop over entry
  cout << "infomation is writing. Please wait for a while" << endl;
  cout << "infomation is writing. Please wait for a while" << endl;
  Int_t njet_need =3; // the number of at least jets of semileptonic final state
  for (Int_t entry =0; entry < chain.GetEntries(); entry++) {

    chain.GetEntry(entry);
    /////////////////////////////////////////////////////
    // get information for final state at detector level
    //select satisfied leptons(electron or muon)
    TLorentzVector p4_lepton[18];
    nlepton=nMuon+nElectron;
    bool lepton_flag = false; // if true pass the selction
    int num_select1=0, num_select2=0;
    for (int i = 0; i < nlepton; i++) {
      if (i < nElectron) {
        p4_lepton[i].SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        lepton_charge[i]=Electron_charge[i];
        if(Electron_cutBased[i]>=2&& abs(Electron_eta[i]) <2.4 && (abs(Electron_eta[i])<1.4442
        	||abs(Electron_eta[i])>1.5660)&&Electron_pt[i]>15)
        {
        	if((abs(Electron_deltaEtaSC[i]+Electron_eta[i])<1.479&&abs(Electron_dxy[i])<0.05&&abs(Electron_dz[i])<0.1)
			    ||(abs(Electron_deltaEtaSC[i]+Electron_eta[i])>=1.479&&abs(Electron_dxy[i])<0.1&&abs(Electron_dz[i])<0.2))
        	{
	          num_select1++;
	          if(Electron_cutBased[i]==4&& abs(Electron_eta[i]) <2.4 && (abs(Electron_eta[i])<1.4442
	        	||abs(Electron_eta[i])>1.5660)&&Electron_pt[i]>30)
	            {  num_select2++;
	               mom_lep = p4_lepton[i];       // the lepton momenton for reconstrut
	               LepCharge = lepton_charge[i];
	               lepton_flag=true;
	              }
	        }      
        }
    }       
     else {
        p4_lepton[i].SetPtEtaPhiM(Muon_pt[i-nElectron],Muon_eta[i-nElectron],Muon_phi[i-nElectron],Muon_mass[i-nElectron]);
        lepton_charge[i]=Muon_charge[i-nElectron];  
        if(Muon_looseId[i-nElectron]==1&&Muon_pfRelIso04_all[i-nElectron]<0.25&&Muon_pt[i-nElectron]>15&&abs(Muon_eta[i-nElectron])<2.4)
          {    num_select1++;
             if(Muon_tightId[i-nElectron]==1&&Muon_pfRelIso04_all[i-nElectron]<0.15&&Muon_pt[i-nElectron]>30&&abs(Muon_eta[i-nElectron])<2.4)
              { num_select2++;
                mom_lep = p4_lepton[i];       // the lepton momenton for reconstruct
                LepCharge = lepton_charge[i]; //the lepton charge for reconstruct
                lepton_flag=true;
                }
          }

        }
      if(num_select1 > 1) {
          lepton_flag=false;
          break;
      }
    }
    nevents++;
    rawtree->Fill();
////////////////////////////////////////////////////////////////////
// select satisfied jets
    nBtag = 0;   // count number of bjet among all the jets
    jet_num = 0; // count number fot jets satisfy the selection criteria
    bool jet_flag = false; // if true pass the selection
    btag_num=0;
    if(lepton_flag==true){
      for (int i = 0; i < nJet; i++) {
          mom_jets[i].SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i],Jet_mass[i]);     
          if (abs(Jet_eta[i]) < 2.4 && Jet_pt[i] > 30 && Jet_jetId[i]>=6 && 
          mom_jets[i].DeltaR(mom_lep)>0.4 ) {
            mom_jets[i].SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i],Jet_mass[i]);
            jet_index[jet_num] = i;
            jet_num = jet_num + 1;
  
            if (Jet_btagDeepB[i] > 0.45) {
              Jet_btaged[i] = 1;
              nBtag++;
            } 
            else
              Jet_btaged[i] = 0;
            btag_num = btag_num + Jet_btaged[i];

        }
      }
      
    }

    if (jet_num >= njet_need && nBtag == 0){
        jet_flag = true;
    
        for(int i=0;i<jet_num;i++){
          jet_btaged[i]=Jet_btaged[jet_index[i]];
          jet_eta[i]=Jet_eta[jet_index[i]];
          jet_pt[i]=Jet_pt[jet_index[i]];
          jet_mass[i]=Jet_mass[jet_index[i]];
          jet_phi[i]=Jet_phi[jet_index[i]];
          jet_btagCSVV2[i]=Jet_btagCSVV2[jet_index[i]];
          jet_btagDeepB[i]=Jet_btagDeepB[jet_index[i]];
          jet_partonFlavour[i]=Jet_partonFlavour[jet_index[i]];
        }
        for (int i = 0; i < nlepton; i++) {
          lepton_pt[i] = p4_lepton[i].Pt();
          lepton_eta[i] = p4_lepton[i].Eta();
          lepton_phi[i] = p4_lepton[i].Phi();
          lepton_mass[i] = p4_lepton[i].M();
        }

    }
    
   //////////////////////////////////////////////////////////////////
    // select ttbar semiletopnic final state
    // select satisfied events
    if (jet_flag == true && lepton_flag == true && PV_npvsGood >=1 ) {
     //if (jet_flag == true ) {

      nu_px = MET_pt * cos(MET_phi);
      nu_py = MET_pt * sin(MET_phi);
      MtW=sqrt(2*(mom_lep.Pt()*MET_pt-mom_lep.Px()*nu_px-mom_lep.Py()*nu_py));
      
        recons_tt();
        if( minimum < 170.0  ){ 
        if(jet_num >= 4)
            count_4jet++;
        else if(jet_num==3)
            count_3jet++; 
        
        mytree->Fill();
        nevents2++;
      } // end of cut of minimum
    } //end of satisfied lepton and jets criteria  

  } // end loop over entries
  file->cd();
  mytree->Write();
  rawtree->Write();
 // delete mytree;
 // delete rawtree;
  cout << inputFile << " has " << chain.GetEntries() << " events" << endl;
  cout << output << " is created" << endl;
  cout << nevents << " events are written into "<< "rawtree." << endl;
  cout << nevents2 << " events are written into "<< "mytree." << endl;
  cout<<count_4jet<<" events have at least 4 jets"<<endl;
  cout<<count_3jet<<" events have 3 jets"<<endl;
  file->Close();

    //loop for different systematics
  //upmytree->Write("");

}