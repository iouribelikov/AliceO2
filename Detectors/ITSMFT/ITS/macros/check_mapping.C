#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#endif

void check_mapping(std::string clusfile = "o2clus_its.root", std::string dictfile = "ITSdictionary.bin")
{
  const int Layer=5;

  // Geometry
  o2::base::GeometryManager::loadGeometry("");
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot,
                                            o2::math_utils::TransformType::L2G)); // request cached transforms

  const int nStav=gman->getNumberOfStaves(Layer);
  const int nChip=gman->getNumberOfChipsPerStave(Layer);
    
  std::cout<<nStav<<' '<<nChip<<'\n';

  std::string tit="Layer" + std::to_string(Layer);
  TH1F *h[nStav];
  for (int i=0; i<nStav; i++) {
    std::string title = tit + "_Stave" + std::to_string(i);
    h[i] = new TH1F(title.data(),title.data(),nChip,-0.5,nChip-0.5);
  }

  //TH2F *hl = new TH2F("hl","x vs z",300,-1.5,1.5,160,-0.8,0.8);
  TH2F *htl = new TH2F("htl","Low y vs z", 49,-73,75.1,16*32,-0.1,2.68);
  TH2F *htu = new TH2F("htu","Up  y vs z", 49,-73,75.1,16*32,-2.68,0.1);
  //TH2F *htu = new TH2F("htu","Up  y vs z", 49,-73,75.1,8*32,-2.68,2.68);
  
  // Dictionary
  if (dictfile.empty()) {
    dictfile = o2::base::NameConf::getDictionaryFileName(o2::detectors::DetID::ITS, "", ".bin");
  }
  o2::itsmft::TopologyDictionary dict;
  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    LOG(INFO) << "Running with dictionary: " << dictfile.c_str();
    dict.readBinaryFile(dictfile);
  } else {
    LOG(INFO) << "Running without dictionary !";
  }

  // Clusters
  TFile fileC(clusfile.data());
  TTree* clusTree = (TTree*)fileC.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);

  // Pixel patterns
  std::vector<unsigned char>* patternsPtr = nullptr;
  auto pattBranch = clusTree->GetBranch("ITSClusterPatt");
  if (pattBranch) {
    pattBranch->SetAddress(&patternsPtr);
  }

  // ROFrecords
  std::vector<o2::itsmft::ROFRecord> *rofRecVecP = nullptr;
  clusTree->SetBranchAddress("ITSClustersROF", &rofRecVecP);

  int ne=clusTree->GetEntries();
  for (int e=0; e<ne; e++) {
    clusTree->GetEntry(e);
    auto &clusters = *clusArr;
    
    auto pattIt = patternsPtr->cbegin();
    for (const auto &c : clusters) {
      auto chipID = c.getSensorID();
      int lay=gman->getLayer(chipID);
      if (lay != Layer) continue;

      auto pattID = c.getPatternID();
      if (pattID == o2::itsmft::CompCluster::InvalidPatternID ) {
          continue;
      }
      o2::math_utils::Point3D<float> locC;
      if (dict.isGroup(pattID)) {
          o2::itsmft::ClusterPattern patt(pattIt);
	  //npix = patt.getRowSpan()*patt.getColumnSpan(); // Just the bbox...
          locC = o2::itsmft::TopologyDictionary::getClusterCoordinates(c, patt);
      } else { 
	  //npix = dict.getNpixels(pattID);
          locC = dict.getClusterCoordinates(c);
      }
      
      // Inverse transformation to the local --> tracking
      auto traC = gman->getMatrixT2L(chipID) ^ locC;

      //std::cout<<traC.X()<<' '<<traC.Y()<<' '<<traC.Z()<<'\n';

      if (traC.X()<34.5)
      htl->Fill(traC.Z(), traC.Y());
      else
      htu->Fill(traC.Z(), traC.Y());
      
      auto sta = gman->getStave(chipID);
      auto idx = gman->getChipIdInStave(chipID);
      h[sta]->Fill(idx);

    }
  }

  new TCanvas;
  htl->Draw("colz");
  
  new TCanvas;
  htu->Draw("colz");
  

  TFile f("h.root","new");
  for (int i=0; i<nStav; i++) {
     f.WriteTObject(h[i]);
  }
  f.Close();

}

