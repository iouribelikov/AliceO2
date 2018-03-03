// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file LinearVertex.h
/// \brief Definition of a linearly fitted vertex
/// \author iouri.belikov@cern.ch

#ifndef ALICEO2_ITS_LINEARVERTEX_H
#define ALICEO2_ITS_LINEARVERTEX_H

#include <array>

namespace o2
{
namespace ITS
{
class LinearVertex
{
 public:
  LinearVertex() = default;
  LinearVertex(const LinearVertex& t) = default;
  LinearVertex& operator=(const LinearVertex& tr) = default;
  ~LinearVertex() = default;

  auto getX() const { return mX; }
  auto getY() const { return mY; }
  auto getZ() const { return mZ; }
  const auto& getCovariance() const { return mCov; }
  auto getCovX2() const { return mCov[0]; }
  auto getCovXY() const { return mCov[1]; }
  auto getCovY2() const { return mCov[2]; }
  auto getCovXZ() const { return mCov[3]; }
  auto getCovYZ() const { return mCov[4]; }
  auto getCovZ2() const { return mCov[5]; }
  auto getNumberOfProngs() const { return mProngs; }
  auto getChi2() const { return mChi2; }
  
  Bool_t update(const std::array<Double_t, 3>& p, const std::array<Double_t, 3>& v, Double_t sy2, Double_t sz2);

 private:
  Double_t mX = 0;                                  ///< vertex position
  Double_t mY = 0;                                  ///< vertex position
  Double_t mZ = 0;                                  ///< vertex position
  std::array<Double_t, 6> mCov{ 9, 0, 9, 0, 0, 9 }; ///< vertex covariance matrix
  Int_t mProngs = 0;                                ///< number of prongs
  Double_t mChi2 = 0;                               ///< chi2
  static constexpr Double_t mMaxChi2 = 4;           ///< maximal accepted chi2 increment

  ClassDef(LinearVertex, 1)
};
} // namespace ITS
} // namespace o2
#endif /* ALICEO2_ITS_LINEARVERTEX_H */
