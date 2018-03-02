// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrivialVertexer.cxx
/// \brief Implementation of the ITS trivial vertex finder

#include <map>

#include "FairLogger.h"

#include "ITSBase/GeometryTGeo.h"
#include "MathUtils/Cartesian3D.h"
#include "ITSReconstruction/LinearVertex.h"
#include "ITSReconstruction/TrivialVertexer.h"
#include "DataFormatsITSMFT/Cluster.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::ITSMFT;
using namespace o2::ITS;

TrivialVertexer::TrivialVertexer() = default;

TrivialVertexer::~TrivialVertexer() = default;

void TrivialVertexer::process(const std::vector<Cluster>& clusters, std::vector<std::array<Double_t, 3>>& vertices)
{
  if (mClsLabels == nullptr) {
    LOG(INFO) << "TrivialVertexer::process() : "
              << "No cluster labels available ! Running with a default MC vertex..." << FairLogger::endl;
    vertices.emplace_back(std::array<Double_t, 3>{ 0., 0., 0. });
    return;
  }

  using layer = std::vector<Int_t>;
  std::map<Int_t, std::pair<layer, layer>> events;
  
  auto gman = o2::ITS::GeometryTGeo::Instance();

  // Separate clusters coming from different MC events
  for (Int_t i = 0; i < clusters.size(); ++i) {
    auto mclab = (mClsLabels->getLabels(i))[0];
    if (mclab.getTrackID() == -1)
      continue; // noise
    auto id = mclab.getEventID();
    auto &event = events[id];
    const auto &c = clusters[i];
    auto r = c.getX();
    if (TMath::Abs(r - 2.2) < 0.5) { event.first.push_back(i); continue; }
    if (TMath::Abs(r - 3.0) < 0.5) event.second.push_back(i);
  }

  for (const auto &event : events) {
    auto mcEv = event.first;
    LinearVertex vtx;
    const auto &layer0 = event.second.first;
    for (auto i0 : layer0) {
        const auto &c0 = clusters[i0];
        auto mclab0 = (mClsLabels->getLabels(i0))[0];
        auto lab0 = mclab0.getTrackID();
        const auto &layer1 = event.second.second;
        for (auto i1 : layer1) {
           const auto &c1 = clusters[i1];
           auto mclab1 = (mClsLabels->getLabels(i1))[0];
           auto lab1 = mclab1.getTrackID();
	   if (lab0 != lab1) continue;
	   const auto p0 = c0.getXYZGloRot(*gman); 
	   const auto p1 = c1.getXYZGloRot(*gman); 
           std::array<Double_t, 3> p{ p0.X(), p0.Y(), p0.Z() };
           std::array<Double_t, 3> v{ p1.X() - p0.X(), p1.Y() - p0.Y(), p1.Z() - p0.Z() };
           auto sy2 = (c1.getSigmaY2() + c0.getSigmaY2());
           auto sz2 = (c1.getSigmaZ2() + c0.getSigmaZ2());
	   vtx.update(p, v, sy2, sz2);
	}
    }
    auto vx = vtx.getX();
    auto vy = vtx.getY();
    auto vz = vtx.getZ();

    vertices.emplace_back(std::array<Double_t, 3>{ vx, vy, vz });
    LOG(INFO) << "TrivialVertexer::process() : "
              << "MC event #" << mcEv << " with vertex (" << vx << ',' << vy << ',' << vz << ')' << FairLogger::endl;
  }
}
