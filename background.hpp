//////////////////////////
// background.hpp
//////////////////////////
// 
// code components related to background evolution
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: February 2016
//
//////////////////////////

#ifndef BACKGROUND_HEADER
#define BACKGROUND_HEADER

#include <gsl/gsl_integration.h>

double FermiDiracIntegrand(double q, void * w)
{
	return q * q * sqrt(q * q + *(double *)w) / (exp(q) + 1.0l);
}

double FermiDiracIntegral(double &w)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	
	f.function = &FermiDiracIntegrand;
	f.params = &w;
		
	gsl_integration_qng(&f, 0.0l, 24.0l, 5.0e-7, 1.0e-7, &result, &err, &n);
	
	return result;
}


//////////////////////////
// bg_ncdm
//////////////////////////
// Description:
//   computes the background model for ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//
// Note:
//   For optimization, the last value of a is stored in a static variable such that
//   multiple calls at the same value of a will not result in multiple integrations
//   being carried out. This assumes that the cosmological model should not change!
//
// Returns: value for the background model
// 
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo)
{
	double w;
	static double result = -1.0;
	static double a_prev = -1.0;
	
	if (a != a_prev)
	{
		result = 0.0;
		a_prev = a;
		
		for (int p = 0; p < cosmo.num_ncdm; p++)
		{
			w = a * cosmo.m_ncdm[p] / (pow(cosmo.Omega_g * cosmo.h * cosmo.h / 4.48147e-7, 0.25) * cosmo.T_ncdm[p] * 8.61733e-5);
			w *= w;
		
			result += FermiDiracIntegral(w) * cosmo.Omega_ncdm[p] * pow(cosmo.Omega_g * cosmo.h * cosmo.h / 4.48147e-7, 0.25) * cosmo.T_ncdm[p] * 4.77921357e-5 / cosmo.m_ncdm[p];
		}
	}
	
	return result / a;
}


//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
// 
//////////////////////////

double Hconf(const double a, const double fourpiG, const cosmology cosmo)
{	
	return sqrt((2. * fourpiG / 3.) * (((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a) + (cosmo.Omega_Lambda * a * a) + (cosmo.Omega_rad / a / a)));
}


double Omega_m(const double a, const cosmology cosmo) { return cosmo.Omega_m / (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) + cosmo.Omega_Lambda * a * a * a + cosmo.Omega_rad / a); }

double Omega_rad(const double a, const cosmology cosmo) { return (cosmo.Omega_rad + (bg_ncdm(a, cosmo) + cosmo.Omega_cdm + cosmo.Omega_b - cosmo.Omega_m) * a) / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * a + cosmo.Omega_Lambda * a * a * a * a + cosmo.Omega_rad); }

double Omega_Lambda(const double a, const cosmology cosmo) { return cosmo.Omega_Lambda / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a / a / a + cosmo.Omega_Lambda + cosmo.Omega_rad / a / a / a / a); }


//////////////////////////
// rungekutta4bg
//////////////////////////
// Description:
//   integrates the Friedmann equation for the background model using a fourth-order
//   Runge-Kutta method
// 
// Arguments:
//   a          scale factor (will be advanced by dtau)
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   dtau       time step by which the scale factor should be advanced
//
// Returns:
// 
//////////////////////////

void rungekutta4bg(double &a, const double fourpiG, const cosmology cosmo, const double dtau)
{
	double k1a, k2a, k3a, k4a;

	k1a = a * Hconf(a, fourpiG, cosmo); 
	k2a = (a + k1a * dtau / 2.) * Hconf(a + k1a * dtau / 2., fourpiG, cosmo);
	k3a = (a + k2a * dtau / 2.) * Hconf(a + k2a * dtau / 2., fourpiG, cosmo);
	k4a = (a + k3a * dtau) * Hconf(a + k3a * dtau, fourpiG, cosmo);

	a += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
}

#endif

