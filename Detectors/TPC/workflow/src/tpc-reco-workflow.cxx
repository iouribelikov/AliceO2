// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// @file   tpc-reco-workflow.cxx
/// @author Matthias Richter
/// @since  2018-03-15
/// @brief  Basic DPL workflow for TPC reconstruction starting from digits

#include "Framework/WorkflowSpec.h"
#include "Framework/DeviceSpec.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CompletionPolicy.h"
#include "Framework/DispatchPolicy.h"
#include "Framework/PartRef.h"
#include "TPCWorkflow/RecoWorkflow.h"
#include "DataFormatsTPC/TPCSectorHeader.h"
#include "Algorithm/RangeTokenizer.h"
#include "SimConfig/ConfigurableParam.h"

#include <string>
#include <stdexcept>
#include <unordered_map>

// we need a global variable to know how many parts are expected in the completion policy check
bool gDoMC = true;
bool gDispatchPrompt = true;

// add workflow options, note that customization needs to be declared before
// including Framework/runDataProcessing
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  using namespace o2::framework;

  std::vector<ConfigParamSpec> options{
    {"input-type", VariantType::String, "digits", {"digitizer, digits, raw, clusters"}},
    {"output-type", VariantType::String, "tracks", {"digits, raw, clusters, tracks"}},
    {"ca-clusterer", VariantType::Bool, false, {"Use clusterer of GPUCATracking"}},
    {"disable-mc", VariantType::Bool, false, {"disable sending of MC information"}},
    {"tpc-sectors", VariantType::String, "0-35", {"TPC sector range, e.g. 5-7,8,9"}},
    {"tpc-lanes", VariantType::Int, 1, {"number of parallel lanes up to the tracker"}},
    {"dispatching-mode", VariantType::String, "prompt", {"determines when to dispatch: prompt, complete"}},
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings (e.g.: 'TPCHwClusterer.peakChargeThreshold=4;...')"}},
    {"configFile", VariantType::String, "", {"configuration file for configurable parameters"}}};

  std::swap(workflowOptions, options);
}

// customize dispatch policy, dispatch immediately what is ready
void customize(std::vector<o2::framework::DispatchPolicy>& policies)
{
  // we customize all devices to dispatch data immediately
  policies.push_back({"prompt-for-all", [](auto const&) { return gDispatchPrompt; }, o2::framework::DispatchPolicy::DispatchOp::WhenReady});
}

// customize clusterers and cluster decoders to process immediately what comes in
void customize(std::vector<o2::framework::CompletionPolicy>& policies)
{
  // we customize the processors to consume data as it comes
  policies.push_back({"tpc-sector-processors",
                      [](o2::framework::DeviceSpec const& spec) {
                        // the decoder should process immediately
                        bool apply = spec.name.find("decoder") != std::string::npos || spec.name.find("clusterer") != std::string::npos;
                        if (apply) {
                          LOG(INFO) << "Applying completion policy 'consume' to device " << spec.name;
                        }
                        return apply;
                      },
                      [](gsl::span<o2::framework::PartRef const> const& inputs) {
                        o2::framework::CompletionPolicy::CompletionOp policy = o2::framework::CompletionPolicy::CompletionOp::Wait;
                        int nValidParts = 0;
                        bool eod = false;
                        if (!gDoMC) {
                          for (auto const& part : inputs) {
                            if (part.header == nullptr) {
                              continue;
                            }
                            nValidParts++;
                            auto const* header = o2::header::get<o2::header::DataHeader*>(part.header->GetData(), part.header->GetSize());
                            auto const* sectorHeader = o2::header::get<o2::tpc::TPCSectorHeader*>(part.header->GetData(), part.header->GetSize());
                            if (sectorHeader && sectorHeader->sector < 0) {
                              eod = true;
                            }
                          }
                          if (nValidParts == inputs.size()) {
                            policy = o2::framework::CompletionPolicy::CompletionOp::Consume;
                          } else if (eod == false) {
                            policy = o2::framework::CompletionPolicy::CompletionOp::Consume;
                          }
                          return policy;
                        }
                        using IndexType = o2::header::DataHeader::SubSpecificationType;

                        std::set<IndexType> matchIndices;
                        for (auto const& part : inputs) {
                          if (part.header == nullptr) {
                            continue;
                          }
                          nValidParts++;
                          auto const* header = o2::header::get<o2::header::DataHeader*>(part.header->GetData(), part.header->GetSize());
                          assert(header != nullptr);
                          auto haveAlready = matchIndices.find(header->subSpecification);
                          if (haveAlready != matchIndices.end()) {
                            // inputs should be data-mc pairs if the index is already in the list
                            // the pair is complete and we can remove the index
                            matchIndices.erase(haveAlready);
                          } else {
                            // store the index in order to check if there is a pair
                            matchIndices.emplace(header->subSpecification);
                          }
                          auto const* sectorHeader = o2::header::get<o2::tpc::TPCSectorHeader*>(part.header->GetData(), part.header->GetSize());
                          if (sectorHeader && sectorHeader->sector < 0) {
                            eod = true;
                          }
                        }
                        if (nValidParts == inputs.size()) {
                          policy = o2::framework::CompletionPolicy::CompletionOp::Consume;
                        } else if (matchIndices.size() == 0 && eod == false) {
                          policy = o2::framework::CompletionPolicy::CompletionOp::Consume;
                        }
                        return policy;
                      }});
}

#include "Framework/runDataProcessing.h" // the main driver

using namespace o2::framework;

/// The workflow executable for the stand alone TPC reconstruction workflow
/// The basic workflow for TPC reconstruction is defined in RecoWorkflow.cxx
/// and contains the following default processors
/// - digit reader
/// - clusterer
/// - cluster raw decoder
/// - CA tracker
///
/// The default workflow can be customized by specifying input and output types
/// e.g. digits, raw, tracks.
///
/// MC info is processed by default, disabled by using command line option `--disable-mc`
///
/// This function hooks up the the workflow specifications into the DPL driver.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto tpcSectors = o2::RangeTokenizer::tokenize<int>(cfgc.options().get<std::string>("tpc-sectors"));
  // the lane configuration defines the subspecification ids to be distributed among the lanes.
  std::vector<int> laneConfiguration;
  auto nLanes = cfgc.options().get<int>("tpc-lanes");
  auto inputType = cfgc.options().get<std::string>("input-type");
  if (inputType == "digitizer") {
    // the digitizer is using a different lane setup so we have to force this for the moment
    laneConfiguration.resize(nLanes);
    std::iota(laneConfiguration.begin(), laneConfiguration.end(), 0);
  } else {
    laneConfiguration = tpcSectors;
  }

  // set up configuration
  o2::conf::ConfigurableParam::updateFromFile(cfgc.options().get<std::string>("configFile"));
  o2::conf::ConfigurableParam::updateFromString(cfgc.options().get<std::string>("configKeyValues"));
  o2::conf::ConfigurableParam::writeINI("o2tpcrecoworkflow_configuration.ini");

  gDoMC = not cfgc.options().get<bool>("disable-mc");
  auto dispmode = cfgc.options().get<std::string>("dispatching-mode");
  gDispatchPrompt = !(dispmode == "single");
  return o2::tpc::reco_workflow::getWorkflow(tpcSectors,                                     // sector configuration
                                             laneConfiguration,                              // lane configuration
                                             gDoMC,                                          //
                                             nLanes,                                         //
                                             inputType,                                      //
                                             cfgc.options().get<std::string>("output-type"), //
                                             cfgc.options().get<bool>("ca-clusterer")        //
  );
}
