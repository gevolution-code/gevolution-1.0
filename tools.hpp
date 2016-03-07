//////////////////////////
// tools.hpp
//////////////////////////
// 
// Collection of analysis tools for gevolution
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: January 2016
//
//////////////////////////

#ifndef TOOLS_HEADER
#define TOOLS_HEADER

#ifndef Cplx
#define Cplx Imag
#endif

#define KTYPE_GRID      0
#define KTYPE_LINEAR    1

using namespace std;
using namespace LATfield2;


#ifdef FFT3D
//////////////////////////
// extractCrossSpectrum
//////////////////////////
// Description:
//   generates the cross spectrum for two Fourier images
// 
// Arguments:
//   fld1FT     reference to the first Fourier image for which the cross spectrum should be extracted
//   fld2FT     reference to the first Fourier image for which the cross spectrum should be extracted
//   kbin       allocated array that will contain the central k-value for the bins
//   power      allocated array that will contain the average power in each bin
//   kscatter   allocated array that will contain the k-scatter for each bin
//   pscatter   allocated array that will contain the scatter in power for each bin
//   occupation allocated array that will count the number of grid points contributing to each bin
//   numbin     number of bins (minimum size of all arrays)
//   ktype      flag indicating which definition of momentum to be used
//                  0: grid momentum
//                  1: linear (default)
//
// Returns:
// 
//////////////////////////

void extractCrossSpectrum(Field<Cplx> & fld1FT, Field<Cplx> & fld2FT, Real * kbin, Real * power, Real * kscatter, Real * pscatter, int * occupation, const int numbins, const bool deconvolve = true, const int ktype = KTYPE_LINEAR)
{
	int i, weight;
	const int linesize = fld1FT.lattice().size(1);
	Real * typek2;
	Real * sinc;
	Real k2max, k2, s;
	rKSite k(fld1FT.lattice());
	Cplx p;
	
	typek2 = (Real *) malloc(linesize * sizeof(Real));
	sinc = (Real *) malloc(linesize * sizeof(Real));
	
	if (ktype == KTYPE_GRID)
	{
		for (i = 0; i < linesize; i++)
		{
			typek2[i] = 2. * (Real) linesize * sin(M_PI * (Real) i / (Real) linesize);
			typek2[i] *= typek2[i];
		}
	}
	else
	{
		for (i = 0; i <= linesize/2; i++)
		{
			typek2[i] = 2. * M_PI * (Real) i;
			typek2[i] *= typek2[i];
		}
		for (; i < linesize; i++)
		{
			typek2[i] = 2. * M_PI * (Real) (linesize-i);
			typek2[i] *= typek2[i];
		}
	}
	
	sinc[0] = 1.;
	if (deconvolve)
	{
		for (i = 1; i <= linesize / 2; i++)
		{
			sinc[i] = sin(M_PI * (float) i / (float) linesize) * (float) linesize / (M_PI * (float) i);
		}
	}
	else
	{
		for (i = 1; i <= linesize / 2; i++)
		{
			sinc[i] = 1.;
		}
	}
	for (; i < linesize; i++)
	{
		sinc[i] = sinc[linesize-i];
	}
	
	k2max = 3. * typek2[linesize/2];
	
	for (i = 0; i < numbins; i++)
	{
		kbin[i] = 0.;
		power[i] = 0.;
		kscatter[i] = 0.;
		pscatter[i] = 0.;
		occupation[i] = 0;
	}
	
	for (k.first(); k.test(); k.next())
	{
		if (k.coord(0) == 0 && k.coord(1) == 0 && k.coord(2) == 0)
			continue;
		else if (k.coord(0) == 0)
			weight = 1;
		else if ((k.coord(0) == linesize/2) && (linesize % 2 == 0))
			weight = 1;
		else
			weight = 2;
			
		k2 = typek2[k.coord(0)] + typek2[k.coord(1)] + typek2[k.coord(2)];
		s = sinc[k.coord(0)] * sinc[k.coord(1)] * sinc[k.coord(2)];
		s *= s;
		
		if (fld1FT.symmetry() == LATfield2::symmetric)
		{
			p = fld1FT(k, 0, 1) * fld2FT(k, 0, 1).conj();
			p += fld1FT(k, 0, 2) * fld2FT(k, 0, 2).conj();
			p += fld1FT(k, 1, 2) * fld2FT(k, 1, 2).conj();
			p *= 2.;
			p += fld1FT(k, 0, 0) * fld2FT(k, 0, 0).conj();
			p += fld1FT(k, 1, 1) * fld2FT(k, 1, 1).conj();
			p += fld1FT(k, 2, 2) * fld2FT(k, 2, 2).conj();
		}
		else
		{
			p = Cplx(0., 0.);
			for (i = 0; i < fld1FT.components(); i++)
				p += fld1FT(k, i) * fld2FT(k, i).conj();
		}
		
		i = (int) floor((double) ((Real) numbins * sqrt(k2 / k2max)));
		if (i < numbins) 
		{
			kbin[i] += weight * sqrt(k2);
			kscatter[i] += weight * k2;
			power[i] += weight * p.real() * k2 * sqrt(k2) / s;
			pscatter[i] += weight * p.real() * p.real() * k2 * k2 * k2 / s / s;
			occupation[i] += weight;
		}
	}
	
	free(typek2);
	free(sinc);
	
	parallel.sum<Real>(kbin, numbins);
	parallel.sum<Real>(kscatter, numbins);
	parallel.sum<Real>(power, numbins);
	parallel.sum<Real>(pscatter, numbins);
	parallel.sum<int>(occupation, numbins);
	
	for (i = 0; i < numbins; i++)
	{
		if (occupation[i] > 0)
		{
			kscatter[i] = sqrt(kscatter[i] * occupation[i] - kbin[i] * kbin[i]) / occupation[i];
			if (!isfinite(kscatter[i])) kscatter[i] = 0.;
			kbin[i] = kbin[i] / occupation[i];
			power[i] /= occupation[i];
			pscatter[i] = sqrt(pscatter[i] / occupation[i] - power[i] * power[i]);
			if (!isfinite(pscatter[i])) pscatter[i] = 0.;
		}
	}
}


//////////////////////////
// extractPowerSpectrum
//////////////////////////
// Description:
//   generates the power spectrum for a Fourier image
// 
// Arguments:
//   fldFT      reference to the Fourier image for which the power spectrum should be extracted
//   kbin       allocated array that will contain the central k-value for the bins
//   power      allocated array that will contain the average power in each bin
//   kscatter   allocated array that will contain the k-scatter for each bin
//   pscatter   allocated array that will contain the scatter in power for each bin
//   occupation allocated array that will count the number of grid points contributing to each bin
//   numbin     number of bins (minimum size of all arrays)
//   ktype      flag indicating which definition of momentum to be used
//                  0: grid momentum
//                  1: linear (default)
//
// Returns:
// 
//////////////////////////

void extractPowerSpectrum(Field<Cplx> & fldFT, Real * kbin, Real * power, Real * kscatter, Real * pscatter, int * occupation, const int numbins, const bool deconvolve = true, const int ktype = KTYPE_LINEAR)
{
	extractCrossSpectrum(fldFT, fldFT, kbin, power, kscatter, pscatter, occupation, numbins, deconvolve, ktype);
}
#endif


//////////////////////////
// writePowerSpectrum
//////////////////////////
// Description:
//   writes power spectra as tabulated data into ASCII file
// 
// Arguments:
//   kbin           array containing the central values of k for each bin
//   power          array containing the central values of P(k) for each bin
//   kscatter       array containing the statistical error on k for each bin
//   pscatter       array containing the statistical error on P(k) for each bin
//   occupation     array containing the number of k-modes contributing to each bin
//   numbins        total number of bins (length of the arrays)
//   rescalek       unit conversion factor for k
//   rescalep       unit conversion factor for P(k)
//   filename       output file name
//   description    descriptive header
//   a              scale factor for this spectrum
//
// Returns:
// 
//////////////////////////

void writePowerSpectrum(Real * kbin, Real * power, Real * kscatter, Real * pscatter, int * occupation, const int numbins, const Real rescalek, const Real rescalep, const char * filename, const char * description, const double a)
{
	if (parallel.isRoot())
	{
		FILE * outfile = fopen(filename, "w");
		if (outfile == NULL)
		{
			cout << " error opening file for power spectrum output!" << endl;
		}
		else
		{
			fprintf(outfile, "# %s\n", description);
			fprintf(outfile, "# redshift z=%f\n", (1./a)-1.);
			fprintf(outfile, "# k              Pk             sigma(k)       sigma(Pk)      count\n");
			for (int i = 0; i < numbins; i++)
			{
				if (occupation[i] > 0)
					fprintf(outfile, "  %e   %e   %e   %e   %d\n", kbin[i]/rescalek, power[i]/rescalep, kscatter[i]/rescalek, pscatter[i]/rescalep/ sqrt(occupation[i]), occupation[i]);
			}
			fclose(outfile);
		}
	}
}


//////////////////////////
// computeVectorDiagnostics
//////////////////////////
// Description:
//   computes some diagnostics for the spin-1 perturbation
// 
// Arguments:
//   Bi         reference to the real-space vector field to analyze
//   mdivB      will contain the maximum value of the divergence of Bi
//   mcurlB     will contain the maximum value of the curl of Bi
//
// Returns:
// 
//////////////////////////

void computeVectorDiagnostics(Field<Real> & Bi, Real & mdivB, Real & mcurlB)
{
	Real b1, b2, b3, b4;
	const Real linesize = (Real) Bi.lattice().sizeLocal(0);
	Site x(Bi.lattice());
	
	mdivB = 0.;
	mcurlB = 0.;
	
	for (x.first(); x.test(); x.next())
	{
		b1 = fabs((Bi(x,0)-Bi(x-0,0)) + (Bi(x,1)-Bi(x-1,1)) + (Bi(x,2)-Bi(x-2,2))) * linesize;
		if (b1 > mdivB) mdivB = b1;
		b1 = 0.5 * (Bi(x,0) + Bi(x+0,1) - Bi(x+1,0) - Bi(x,1) + Bi(x+2,0) + Bi(x+0+2,1) - Bi(x+1+2,0) - Bi(x+2,1)) * linesize;
		b2 = 0.5 * (Bi(x,0) + Bi(x+0,2) - Bi(x+2,0) - Bi(x,2) + Bi(x+1,0) + Bi(x+0+1,2) - Bi(x+2+1,0) - Bi(x+1,2)) * linesize;
		b3 = 0.5 * (Bi(x,2) + Bi(x+2,1) - Bi(x+1,2) - Bi(x,1) + Bi(x+0,2) + Bi(x+2+0,1) - Bi(x+1+0,2) - Bi(x+0,1)) * linesize;
		b4 = sqrt(b1 * b1 + b2 * b2 + b3 * b3);
		if (b4 > mcurlB) mcurlB = b4;
	}
	
	parallel.max<Real>(mdivB);
	parallel.max<Real>(mcurlB);
}


//////////////////////////
// computeTensorDiagnostics
//////////////////////////
// Description:
//   computes some diagnostics for the spin-2 perturbation
// 
// Arguments:
//   hij        reference to the real-space tensor field to analyze
//   mdivh      will contain the maximum value of the divergence of hij
//   mtraceh    will contain the maximum value of the trace of hij
//   mnormh     will contain the maximum value of the norm of hij
//
// Returns:
// 
//////////////////////////

void computeTensorDiagnostics(Field<Real> & hij, Real & mdivh, Real & mtraceh, Real & mnormh)
{
	Real d1, d2, d3;
	const Real linesize = (Real) hij.lattice().sizeLocal(0);
	Site x(hij.lattice());
	
	mdivh = 0.;
	mtraceh = 0.;
	mnormh = 0.;
	
	for (x.first(); x.test(); x.next())
	{
		d1 = (hij(x+0, 0, 0) - hij(x, 0, 0) + hij(x, 0, 1) - hij(x-1, 0, 1) + hij(x, 0, 2) - hij(x-2, 0, 2)) * linesize;
		d2 = (hij(x+1, 1, 1) - hij(x, 1, 1) + hij(x, 0, 1) - hij(x-0, 0, 1) + hij(x, 1, 2) - hij(x-2, 1, 2)) * linesize;
		d3 = (hij(x+2, 2, 2) - hij(x, 2, 2) + hij(x, 0, 2) - hij(x-0, 0, 2) + hij(x, 1, 2) - hij(x-1, 1, 2)) * linesize;
		d1 = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
		if (d1 > mdivh) mdivh = d1;
		d1 = fabs(hij(x, 0, 0) + hij(x, 1, 1) + hij(x, 2, 2));
		if (d1 > mtraceh) mtraceh = d1;
		d1 = sqrt(hij(x, 0, 0) * hij(x, 0, 0) + 2. * hij(x, 0, 1) * hij(x, 0, 1) + 2. * hij(x, 0, 2)* hij(x, 0, 2) + hij(x, 1, 1) * hij(x, 1, 1) + 2. * hij(x, 1, 2) * hij(x, 1, 2) + hij(x, 2, 2) * hij(x, 2, 2));
		if (d1 > mnormh) mnormh = d1;
	}
	
	parallel.max<Real>(mdivh);
	parallel.max<Real>(mtraceh);
	parallel.max<Real>(mnormh);
}


//////////////////////////
// hourMinSec
//////////////////////////
// Description:
//   generates formatted output for cpu-time: hh..h:mm:ss.s
// 
// Arguments:
//   seconds    number of seconds
//
// Returns:
//   formatted string
// 
//////////////////////////

string hourMinSec(double seconds)
{
	string output;
	char ptr[20];
	int h, m, s, f;

	h = (int) floor(seconds / 3600.);
	seconds -= 3600. * h;
	m = (int) floor(seconds / 60.);
	seconds -= 60. * m;
	s = (int) floor(seconds);
	seconds -= s;
	f = (int) floor(10. * seconds);
	sprintf(ptr, "%d:%02d:%02d.%d", h, m, s, f);

	output.reserve(20);
	output.assign(ptr);

	return output;
}

#endif
