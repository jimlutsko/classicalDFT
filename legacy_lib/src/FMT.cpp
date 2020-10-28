#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "FMT.h"
#include "Enskog.h"


ostringstream Eta_Too_Large_Exception::cnvt;


FundamentalMeasures FMT::getWeightedDensities(long i, vector<Species*> &allSpecies)
{
  // the weighted densities for the lattice site under consideration (i).
  // The basic densities are the scalars eta and s, the vector v and the tensor T.
  // Note that s0,s1 etc only differ in constant prefactors if there is only one species present.
  FundamentalMeasures fm;
  
  // Collect the contributions to the various weighted densities at lattice position
  for(Species* &generic_species : allSpecies)
    {
      FMT_Species *species = dynamic_cast<FMT_Species*>(generic_species);
      if(species)
	{      
	  double hsd = species->getHSD();
  
	  fm.eta += species->getEta(i); //Eta(dd).r(i);

	  fm.s0 += species->getS(i)/(hsd*hsd) ;
	  fm.s1 += species->getS(i)/hsd;
	  fm.s2 += species->getS(i);

	  for(int j=0;j<3;j++)
	    {
	      fm.v1[j] += species->getV(j,i)/hsd;
	      fm.v2[j] += species->getV(j,i);
	      for(int k=0;k<3;k++)
		fm.T[j][k] += species->getT(j,k,i);
	    }
	}
    }
  /*
  // touch up the normalization of the tensor density to account for any numerical errors
  double ss = fm.T[0][0]+fm.T[1][1]+fm.T[2][2];

  fm.T[0][0] += (fm.s2-ss)/3;
  fm.T[1][1] += (fm.s2-ss)/3;
  fm.T[2][2] += (fm.s2-ss)/3;
  */
  fm.calculate_derived_quantities();
  
  return fm;
}


double FMT::Phi(const FundamentalMeasures& fm) const
{
  // These are the eta-dependendent cofactors that lie at the heart of FMT
  double eta = fm.eta;
  
  double f1 = f1_(eta);
  double f2 = f2_(eta);
  double f3 = f3_(eta);

  // Now, construct phi.
  double s0 = fm.s0;
  double s1 = fm.s1;
  double s2 = fm.s2;
  double v1_v2 = fm.v1_v2;
  
  double phi = 0;
  phi -= (1/M_PI)*s0*f1; 
  phi += (1/(2*M_PI))*(s1*s2-v1_v2)*f2;
  phi += Phi3(fm)*f3;

  return phi;
}

void FMT::calculate_dPhi_wrt_fundamental_measures(const FundamentalMeasures& fm, FundamentalMeasures& dPhi) const
{
  // These are the eta-dependendent cofactors that lie at the heart of FMT
  double eta = fm.eta;
  
  double f1 = f1_(eta);
  double f2 = f2_(eta);
  double f3 = f3_(eta);

  double f1p = f1p_(eta);
  double f2p = f2p_(eta);
  double f3p = f3p_(eta);
  
// Now, construct the derivatives
  double s0 = fm.s0;
  double s1 = fm.s1;
  double s2 = fm.s2;
  double v1_v2 = fm.v1_v2;
  
  dPhi.eta -= (1/M_PI)*s0*f1p; 
  dPhi.eta += (1/(2*M_PI))*(s1*s2-v1_v2)*f2p;
  dPhi.eta += Phi3(fm)*f3p;

  dPhi.s0 += -(1/M_PI)*f1;
  dPhi.s1 += (1/(2*M_PI))*s2*f2;
  dPhi.s2 += (1/(2*M_PI))*s1*f2;
  dPhi.s2 += dPhi3_dS2(fm)*f3;

  for(int k=0;k<3;k++)
    {
      dPhi.v1[k] += (1/(2*M_PI))*(-fm.v2[k])*f2;
      dPhi.v2[k] += (1/(2*M_PI))*(-fm.v1[k])*f2;
      dPhi.v2[k] += dPhi3_dV2(k, fm)*f3; 	  
    }
	
  if(needsTensor())	
    for(int k=0;k<3;k++)
      for(int j=0;j<3;j++)
	dPhi.T[j][k] = dPhi3_dT(j,k,fm)*f3; 		
}

void FMT::D2Phi(const FundamentalMeasures& fm, vector<vector<double>>& d2Phi) const
{
  // These are the eta-dependendent cofactors that lie at the heart of FMT
  double eta = fm.eta;
  
  double f1 = f1_(eta);
  double f2 = f2_(eta);
  double f3 = f3_(eta);

  double f1p = f1p_(eta);
  double f2p = f2p_(eta);
  double f3p = f3p_(eta);

  double f1pp = f1pp_(eta);
  double f2pp = f2pp_(eta);
  double f3pp = f3pp_(eta);  
  
// Now, construct the derivatives
  double s0 = fm.s0;
  double s1 = fm.s1;
  double s2 = fm.s2;
  double v1_v2 = fm.v1_v2;
  
  //eta-eta
  d2Phi[0][0] -= (1/M_PI)*s0*f1pp; 
  d2Phi[0][0] += (1/(2*M_PI))*(s1*s2-v1_v2)*f2pp;
  d2Phi[0][0] += Phi3(fm)*f3pp;

  //eta-s0
  d2Phi[0][1] -= (1/M_PI)*f1p;

  //eta-s1
  d2Phi[0][2] += (1/(2*M_PI))*(s2)*f2p;

  //eta-s2
  d2Phi[0][3] += (1/(2*M_PI))*(s1)*f2p;
  d2Phi[0][3] += dPhi3_dS2(fm)*f3p;

  //eta-v1
  d2Phi[0][4] += -(1/(2*M_PI))*fm.v2[0]*f2p;
  d2Phi[0][5] += -(1/(2*M_PI))*fm.v2[1]*f2p;
  d2Phi[0][6] += -(1/(2*M_PI))*fm.v2[2]*f2p;

  //eta-v2
  d2Phi[0][7] += -(1/(2*M_PI))*fm.v1[0]*f2p + dPhi3_dV2(0,fm)*f3p;
  d2Phi[0][8] += -(1/(2*M_PI))*fm.v1[1]*f2p + dPhi3_dV2(1,fm)*f3p;
  d2Phi[0][9] += -(1/(2*M_PI))*fm.v1[2]*f2p + dPhi3_dV2(2,fm)*f3p;

  // s0-s0, s0-s1, s1-s1
  d2Phi[1][1] = 0.0;
  d2Phi[1][2] = 0.0;
  d2Phi[2][2] = 0.0;

  // s1-s2
  d2Phi[2][3] = (1/(2*M_PI))*f2;

  //s2-s2
  d2Phi[3][3] = dPhi3_dS2_dS2(fm)*f3;

  // s2-v2
  d2Phi[3][7] += dPhi3_dV2_dS2(0, fm)*f3;
  d2Phi[3][8] += dPhi3_dV2_dS2(1, fm)*f3; 	  
  d2Phi[3][9] += dPhi3_dV2_dS2(2, fm)*f3;

  // v1-v2
  d2Phi[4][7] = -(1/(2*M_PI))*f2;
  d2Phi[5][8] = -(1/(2*M_PI))*f2;
  d2Phi[6][9] = -(1/(2*M_PI))*f2;

  // v2-v2
  d2Phi[7][7] = dPhi3_dV2_dV2(0,0,fm)*f3;
  d2Phi[7][8] = dPhi3_dV2_dV2(0,1,fm)*f3;
  d2Phi[7][9] = dPhi3_dV2_dV2(0,2,fm)*f3;

  d2Phi[8][8] = dPhi3_dV2_dV2(1,1,fm)*f3;
  d2Phi[8][9] = dPhi3_dV2_dV2(1,2,fm)*f3;

  d2Phi[9][9] = dPhi3_dV2_dV2(2,2,fm)*f3;  

 if(needsTensor())
   {
     throw std::runtime_error("D2PHI not implemented for tensor theories");
   }

  
  for(int i=0;i<10;i++)
    for(int j=0;j<i;j++)
      d2Phi[i][j] = d2Phi[j][i];
}

double FMT::dPHI(long i, vector<Species*> &allSpecies)
{
  FundamentalMeasures fm = getWeightedDensities(i, allSpecies);

  double phi = Phi(fm);

  if(etaMax_ > 1.0)
    if(fm.eta > 0.5 && 1-fm.eta < 0.0)
      throw Eta_Too_Large_Exception();
  
  // Also add in the contributions to the derivative of phi (at lattice site i) wrt the various weighted densities
  // (part of the chain-rule evaluation of dPhi/drho(j) = dPhi/deta(i) * deta(i)/drho(j) + ...)
  // Note that at the level of the fundamental measures, we only keep track of one of each class since the others are trivially related
  // by factors of hsd.

  // Here, we fill dPhi with dPhi/d n_{alpha}(i) so I re-use the FundamentalMeasures structure even though the meaining here is different.
  FundamentalMeasures dPhi;
  calculate_dPhi_wrt_fundamental_measures(fm,dPhi);
  
  for(Species* &generic_species : allSpecies)
      generic_species->set_fundamental_measure_derivatives(dPhi, i, needsTensor());

  return phi;
}


#pragma omp declare reduction(SummationPlus: Summation: omp_out += omp_in) 

double FMT::calculateFreeEnergy(vector<Species*> &allSpecies)
{
  // Do FFT of density and compute the fundamental measures by convolution
  for(auto s: allSpecies)
    {
      FMT_Species *f = dynamic_cast<FMT_Species*>(s);
      if(f)
	f->convoluteDensities(needsTensor());
    }
  
  // Now compute the free energy. Here we loop over all lattice sites and compute Phi(r_i) for each one. This presupposes that we did the convolution above. 
  long Ntot = allSpecies.front()->getLattice().Ntot();  
  
  // There is a problem throwing exceptions from an OMP loop - I think that only the affected thread stops and the others continue.
  // So we eat the exception, remember it, and rethrow it when all loops are finished.
  bool hadCatch = false;

  Summation F;
  int chunk = Ntot/20;
  long i;

  //    schedule(static,chunk)		\
  
#pragma omp parallel for			\
  shared( chunk, allSpecies )			\
  private(i)					\
  schedule(auto)				\
  reduction(SummationPlus:F)
  for(i=0;i<Ntot;i++)
    {
      try {
	F += dPHI(i,allSpecies);
      } catch( Eta_Too_Large_Exception &e) {
	hadCatch = true;
      }
    }

  // For the AO species, there is additional work to do for both the free energy and the forces. 
  // Do FFT of density and compute the fundamental measures by convolution
  double FAO = 0;
  for(auto s: allSpecies)
    {
      FMT_AO_Species *fao_species = dynamic_cast<FMT_AO_Species*>(s);
      if(fao_species)
	{
	  FAO += fao_species->free_energy_post_process(needsTensor());

	  // The weights now hold Upsilon for each measure. We now need Upsilon-bar
	  // PSI will be used to hold intermediate results

	  const int Nfmt = (needsTensor ? fao_species->size() : 5);  // number of fmt densities is <= this value. 	  
	  for(int a=0;a<Nfmt;a++)
	    {  
	      for(long pos=0; pos<Ntot;pos++)
		{
		  FundamentalMeasures fm;
		  fao_species->getExtendedWeightedDensity(pos, fm);

		  vector< vector<double> > d2PHI(Nfmt,vector<double>(Nfmt,0.0));
		  D2Phi(fm,d2PHI);
	  
		  double sum = d2PHI[a][0]*fao_species->getEta(pos);

		  // all fucked up: which numbering scheme, which hsd?
		  sum += d2PHI[a][1]*fao_species->getS0(pos);
		  sum += d2PHI[a][2]*fao_species->getS1(pos);
		  sum += d2PHI[a][3]*fao_species->getS2(pos);
		  for(int k=0; k<3; k++)
		    {
		      sum += d2PHI[a][4+k]*fao_species->getV1(pos,k);
		      sum += d2PHI[a][7+k]*fao_species->getV2(pos,k);
		      //		      for(int l=k; l<3; l++)	      
		      //			sum += d2PHI[a][10+tcount]*fao_species->getMeasure(tcount,pos);
		    }
		  fao_species->setPSI(pos,-sum*fao_species->getReservoirDensity());
		}
	      fao_species->computeForceContribution(a);
	    }
	  
	}
    }
  // rethrow exception if it occurred: this messiness is do to the parallel evaluation. 
  if(hadCatch) 
    throw Eta_Too_Large_Exception();
  return F;
}

// Calculate dF[i] = dPhi/drho(i)
//                 = sum_j dV * dPhi(j)/drho(i)
//                 = sum_j dV * sum_alpha dPhi(j)/dn_{alpha}(j)  dn_{alpha}(j)/drho(i)
//                 = sum_j dV * sum_alpha dPhi(j)/dn_{alpha}(j)  w_{alpha}(j,i)
// This is done with convolutions: FMT_Weighted_Density is an array with the index (alpha)
// and holds the weights, w_{alpha}(j,i) and their FFT's AND dPhi(j)/dn_{alpha}(j).
// It therefore FFT's both of these and adds them to dPhi_.Four.
// Once this is done of all alpha, dPhi_ FFT's back to real space and the result is put into dF (with a factor of dV thrown in).

double FMT::calculateFreeEnergyAndDerivatives(vector<Species*> &allSpecies)
{
  double F = 0;
  try {
    F = calculateFreeEnergy(allSpecies);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  // The  derivatives: for each species s  we need the deriv wrt the species' density at each lattice site: dF/d n_{s}(i) 
  double dV = allSpecies.front()->getLattice().dV();

  DFT_FFT dPhi_(allSpecies.front()->getLattice().Nx(),
		allSpecies.front()->getLattice().Ny(),
		allSpecies.front()->getLattice().Nz());
  
  for(auto &s: allSpecies)
    {
      FMT_Species *species = dynamic_cast<FMT_Species*>(s);
      if(species)
	{      
	  dPhi_.Four().zeros();
	  species->Accumulate_dPhi(dPhi_.Four(), needsTensor());
	  
	  dPhi_.do_fourier_2_real();

	  dPhi_.Real().MultBy(dV);
	  species->addToForce(dPhi_.cReal());
	}
    }
  return F*dV;
};


// For this we need sum_J sum_K dV (d2PHI(K)/dn_a(K) dn_b(K)) w_a(K-I) w_b(K-J)v(J)
// Define psi_b(K) = dV conv(w_b,v;K) = sum_J w_b(K-J)V(J) dV
// Then we need sum_K (d2PHI(K)/dn_a(K) dn_b(K)) w_a(K-I) psi_b(K)
// or conv((d2PHI(K)/dn_a(K) dn_b(K))psi_b(K), w; I)
//
void FMT::add_second_derivative(vector<DFT_FFT> &v, vector<DFT_Vec> &d2F, vector<Species*> &allSpecies)
{
  if(needsTensor()) throw std::runtime_error("Second derivatives not implemented for tensor densities");
  if(allSpecies.size() < 1)  throw std::runtime_error("No species for FMT::add_second_derivative to work with");

  const int Nfmt = 10;  // number of fmt densities (tensor density not included)

  // Construct psi

  const Density &density1 = allSpecies[0]->getDensity();

  int Nspecies = allSpecies.size();
  long Ntot = density1.Ntot();
  int Nx    = density1.Nx();
  int Ny    = density1.Ny();
  int Nz    = density1.Nz();
  double dV = density1.dV();

  // Some working space
  DFT_FFT result(Nx,Ny,Nz);

  // Get Psi
  vector<DFT_FFT> Psi(Nfmt);
  for(auto &p: Psi) p.initialize(Nx,Ny,Nz);

  for(int s = 0;s<Nspecies;s++)
    {
      FMT_Species *s2 = dynamic_cast<FMT_Species*>(allSpecies[s]);
      if(!s2) continue; // Not an FMT_Species
      
      double hsd = s2->getHSD();
      
      //the fmt weights have all necessary normalization factors built in, including dV
      bool bConjugate = false;
      s2->convolute_eta_weight_with(v[s], result, bConjugate);
      Psi[0].Real().IncrementBy(result.Real());

      // s0,s1,s2
      s2->convolute_s_weight_with(v[s], result, bConjugate);
      Psi[1].Real().IncrementBy_Scaled_Vector(result.Real(), 1.0/(hsd*hsd));
      Psi[2].Real().IncrementBy_Scaled_Vector(result.Real(), 1.0/hsd);
      Psi[3].Real().IncrementBy(result.Real());

      // v1,v2
      for(int i=0;i<3;i++)
	{
	  s2->convolute_v_weight_with(i, v[s], result, bConjugate);
	  Psi[4+i].Real().IncrementBy_Scaled_Vector(result.Real(), 1.0/hsd);
	  Psi[7+i].Real().IncrementBy(result.Real());
	}      
    }

  // Get Lambda
  vector<DFT_FFT> Lambda(Nfmt);
  for(auto &p: Lambda) p.initialize(Nx,Ny,Nz);

  for(long pos = 0; pos < Ntot; pos++)
    {
      FundamentalMeasures fm = getWeightedDensities(pos, allSpecies);
      vector<vector<double>> d2Phi(Nfmt,vector<double>(Nfmt,0.0));
      D2Phi(fm, d2Phi);

      for(int a=0;a<Nfmt;a++)
	for(int b=0;b<Nfmt;b++)	    
	  Lambda[a].Real().IncrementBy(pos,d2Phi[a][b]*Psi[b].cReal().get(pos));
    }

  // Do the last convolution to get d2F
  for(int a=0;a<Nfmt;a++)
    Lambda[a].do_real_2_fourier();
  
  for(int s = 0;s<Nspecies;s++)
    {
      FMT_Species *s1 = dynamic_cast<FMT_Species*>(allSpecies[s]);
      if(!s1) continue; // Not an FMT_Species

      double hsd = s1->getHSD();      
            
      //the fmt weights have all necessary normalization factors built in, including dV
      bool bConjugate = true;
      
      s1->convolute_eta_weight_with(Lambda[0], result, bConjugate);
      d2F[s].IncrementBy(result.Real());

      // s0
      s1->convolute_s_weight_with(Lambda[1], result, bConjugate);
      d2F[s].IncrementBy_Scaled_Vector(result.Real(), 1.0/(hsd*hsd));

      // s1
      s1->convolute_s_weight_with(Lambda[2], result, bConjugate);
      d2F[s].IncrementBy_Scaled_Vector(result.Real(), 1.0/hsd);

      //s2
      s1->convolute_s_weight_with(Lambda[3], result, bConjugate);
      d2F[s].IncrementBy(result.Real());

      // v1 and v2
      for(int i=0;i<3;i++)
	{
	  s1->convolute_v_weight_with(i, Lambda[4+i], result, bConjugate);
	  d2F[s].IncrementBy_Scaled_Vector(result.Real(), 1.0/hsd);

	  s1->convolute_v_weight_with(i, Lambda[7+i], result, bConjugate);
	  d2F[s].IncrementBy(result.Real());
	}
    }

  for(auto &f: d2F)
    f.MultBy(dV);
}

// Brute-force evaluation of second derivatives
double FMT::d2Phi_dn_dn(int I[3], int si, int J[3], int sj, vector<Species*> &allSpecies)
{
  FMT_Species *s1 = dynamic_cast<FMT_Species*>(allSpecies[si]);
  if(!s1) return 0; // Not an FMT_Species

  FMT_Species *s2 = dynamic_cast<FMT_Species*>(allSpecies[sj]);
  if(!s2) return 0; // Not an FMT_Species  

  const int Nfmt = 10;  // number of fmt densities (tensor density not included)

  // Construct psi

  const Density &density1 = allSpecies[0]->getDensity();

  int Nx    = density1.Nx();
  int Ny    = density1.Ny();
  int Nz    = density1.Nz();
  double dV = density1.dV();

  double f = 0;
  
  for(int ix = 0;ix<Nx;ix++)
    for(int iy = 0;iy<Ny;iy++)
      for(int iz = 0;iz<Nz;iz++)
	{
	  long K = density1.get_PBC_Pos(ix,iy,iz);
	  long KI = density1.get_PBC_Pos(ix-I[0],iy-I[1],iz-I[2]);
	  long KJ = density1.get_PBC_Pos(ix-J[0],iy-J[1],iz-J[2]);
	  
	  FundamentalMeasures fm = getWeightedDensities(K, allSpecies);
	  vector<vector<double>> d2Phi(Nfmt,vector<double>(Nfmt,0.0));
	  D2Phi(fm, d2Phi);	      	      
	  
	  for(int a=0;a<Nfmt;a++)
	    for(int b=0;b<Nfmt;b++)
	      f += d2Phi[a][b]*s1->getExtendedWeight(KI,a)*s2->getExtendedWeight(KJ,b);
	}
  return f*dV;
}
