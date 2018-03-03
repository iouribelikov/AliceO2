// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file LinearVertex.cxx
/// \brief Implementation of a linearly fitted vertex
/// \author iouri.belikov@cern.ch

#include <TMath.h>
#include <TMatrixD.h>

#include "ITSReconstruction/LinearVertex.h"

ClassImp(o2::ITS::LinearVertex)

using namespace o2::ITS;

Bool_t LinearVertex::update(const std::array<Double_t, 3>& p, const std::array<Double_t, 3>& v, Double_t sy2,
                            Double_t sz2)
{
  //--------------------------------------------------------------------
  // Linear update of the vertex parameters
  //--------------------------------------------------------------------

  // Get the vertex weight matrix
  TMatrixD wv(3, 3);
  wv(0, 0) = mCov[0];
  wv(0, 1) = mCov[1];
  wv(0, 2) = mCov[3];
  wv(1, 0) = mCov[1];
  wv(1, 1) = mCov[2];
  wv(1, 2) = mCov[4];
  wv(2, 0) = mCov[3];
  wv(2, 1) = mCov[4];
  wv(2, 2) = mCov[5];
  wv.Invert();
  if (!wv.IsValid())
    return kFALSE;

  // Get the tracklet weight matrix
  TMatrixD wt(3, 3);
  wt(0, 0) = 0.;
  wt(0, 1) = 0.;
  wt(0, 2) = 0.;
  wt(1, 0) = 0.;
  wt(1, 1) = 1 / sy2;
  wt(1, 2) = 0.;
  wt(2, 0) = 0.;
  wt(2, 1) = 0.;
  wt(2, 2) = 1 / sz2; // FIXME : In "parallel" system

  auto l = v[0];
  auto m = v[1];
  auto n = v[2];
  auto phi = TMath::ATan2(m, l);
  auto sp = TMath::Sin(phi);
  auto cp = TMath::Cos(phi);
  auto tgl = n / TMath::Sqrt(l * l + m * m);
  auto cl = 1 / TMath::Sqrt(1. + tgl * tgl);
  auto sl = tgl * cl;

  TMatrixD p2g(3, 3); //"parallel" --> global transformation
  p2g(0, 0) = cp * cl;
  p2g(1, 0) = sp * cl;
  p2g(2, 0) = sl;
  p2g(0, 1) = -sp;
  p2g(1, 1) = cp;
  p2g(2, 1) = 0.;
  p2g(0, 2) = -sl * cp;
  p2g(1, 2) = -sl * sp;
  p2g(2, 2) = cl;
  wt = p2g * wt * TMatrixD(TMatrixD::kTransposed, p2g); // Now, in global system

  // Check the possible chi2 increment
  TMatrixD cv(wv);
  cv += wt;
  cv.Invert();
  if (!cv.IsValid())
    return kFALSE;

  wv = wv*cv*wt;

  auto lmn = TMath::Sqrt(l * l + m * m + n * n);
  auto cosx = l / lmn, cosy = m / lmn, cosz = n / lmn;
  auto pp = (l * (p[0] - mX) + m * (p[1] - mY) + n * (p[2] - mZ)) / lmn;
  auto tx = p[0] - pp * cosx;
  auto ty = p[1] - pp * cosy;
  auto tz = p[2] - pp * cosz;
  Double_t delta[3]{ tx - mX, ty - mY, tz - mZ };
  Double_t chi2 = 0;
  for (Int_t i = 0; i < 3; i++) {
    Double_t s = 0.;
    for (Int_t j = 0; j < 3; j++)
      s += wv(i, j) * delta[j];
    chi2 += s * delta[i];
  }
  if (chi2 > mMaxChi2)
     return kFALSE;

  // Update the vertex covariance
  mCov[0] = cv(0, 0);
  mCov[1] = cv(1, 0);
  mCov[2] = cv(1, 1);
  mCov[3] = cv(2, 0);
  mCov[4] = cv(2, 1);
  mCov[5] = cv(2, 2);

  // Update the vertex position, chi2, and the number of prongs
  TMatrixD cwt(cv, TMatrixD::kMult, wt);
  for (Int_t j = 0; j < 3; j++)
    mX += cwt(0, j) * delta[j];
  for (Int_t j = 0; j < 3; j++)
    mY += cwt(1, j) * delta[j];
  for (Int_t j = 0; j < 3; j++)
    mZ += cwt(2, j) * delta[j];
  mChi2 += chi2;
  mProngs++;

  return kTRUE;
}
