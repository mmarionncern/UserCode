
AdvFilter = cms.EDFilter('AdvLeptonFilter',

  ElectronInput = cms.InputTag('gsfElectrons'),
  MuonInput = cms.InputTag('muons'),
  TauInput = cms.InputTag('selectedPatTaus'),                          
  PhotonInput = cms.InputTag('photons'),
  PFJetInput = cms.InputTag('ak5PFJets'), #   ak5PFJetsL1L2L3                 

  MetInput = cms.InputTag('pfMet'), #pfType1CorrectedMet

### selections follow the convention :
### e -> electron, m -> muon, t -> tau (hadronic), p -> photon, j -> jet, h -> MET
### l -> (electron OR muon OR tau), L -> (electron OR muon)
###
### available ID working point : very loose V, loose L, medium M and tight T
### for electrons and photons -> official cut based IDs including isolation
### for muons -> tight ID corresponds to the nominal VBTF ID
### for taus -> ID corresponds to isolation for very loose, loose and medium IDs,
###             tight ID corresponds to medium ID + loose electron/muon discriminators                        
### for jets -> no very loose Id, official ID provided by the JetMET group
### for MET -> no ID available                         
###
### Examples :
### to select one muon of 20 GeV, the selction should be written as 1m_20
### for one electron of 25 GeV with loose ID : 1e_20L
### to select two muons of 20 GeV : 2m_20_20
### to select two muons or electrons of 20 GeV with one medium ID : 2L_20M_20
### to select three leptons with pt threshold of 5, 10 and 15 GeV, with a tight ID for the 5 GeV lepton : 3l_05T_10_15
### to select two electrons of 20 GeV, one jet of 50 GeV with loose ID and MET of 50 : 2e_20_20_1j_50L_1h_50
###
### WARNING : the filter will run up to 5 kinds of objects, no more. e.g 1h_xx_Ne_[...]_Nt_[...]_Nm_[...]_Nj_[...]
###           but an infinite number of selections can be called                         
###                         
### notice that the pt threshold have to be written with two characters and order of requirements is not important
### be sure that for a given number of objects the selection contains the corresponding number of pT/ID conditions
### notice also that a 1 is needed before the MET tag (h) to respect the encoding format
                         

  Selections = cms.untracked.vstring("2L_20V_20V_1h_20")                       
#,"2l_05T_05V_1m_20","2m_20_20L","2L_20_20_1m_20","2l_10T_15V","3m_20_20_20"
)

