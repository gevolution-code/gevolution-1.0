//////////////////////////
// ic_basic.hpp
//////////////////////////
// 
// basic initial condition generator for gevolution
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: February 2016
//
//////////////////////////

#ifndef IC_BASIC_HEADER
#define IC_BASIC_HEADER

#include "prng_engine.hpp"
#include <gsl/gsl_spline.h>

#ifndef Cplx
#define Cplx Imag
#endif  

#define MAX_LINESIZE 256

using namespace std;
using namespace LATfield2;


// header structure for GADGET-2 files [V. Springel, N. Yoshida, and S.D. White, New Astron. 6 (2001) 79
// and V. Springel, Mon. Not. R. Astron. Soc. 364 (2005) 1105], used only for loading a homogeneous template

#ifndef GADGET2_HEADER
#define GADGET2_HEADER
struct gadget2_header
{
	int npart[6];
	double mass[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];   /* fills to 256 Bytes */
};
#endif


// should be larger than maximum Ngrid
#ifndef HUGE_SKIP
#define HUGE_SKIP   65536
#endif


// trivial move function

void move_particles_ic_basic(double dtau, double lat_resolution, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
	for (int l = 0; l < 3; l++) (*part).pos[l] += dtau*(*part).vel[l];
}


// initializes velocities proportional to gradient of phi

Real initialize_q_ic_basic(double dtau, double lat_resolution, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
#define phi (fields[0])
#define xP (sites[0])
	
	double gradPhi[3] = {0, 0, 0};
	double v2 = 0.;
	
	gradPhi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*phi)(xP+0) - (*phi)(xP));
	gradPhi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*phi)(xP+1) - (*phi)(xP));
	gradPhi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*phi)(xP+2) - (*phi)(xP));
	gradPhi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*phi)(xP+1+0) - (*phi)(xP+1));
	gradPhi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*phi)(xP+1+0) - (*phi)(xP+0));
	gradPhi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*phi)(xP+2+0) - (*phi)(xP+0));
	gradPhi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*phi)(xP+2+0) - (*phi)(xP+2));
	gradPhi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*phi)(xP+2+1) - (*phi)(xP+2));
	gradPhi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*phi)(xP+2+1) - (*phi)(xP+1));
	gradPhi[0] += ref_dist[1] * ref_dist[2] * ((*phi)(xP+2+1+0) - (*phi)(xP+2+1));
	gradPhi[1] += ref_dist[0] * ref_dist[2] * ((*phi)(xP+2+1+0) - (*phi)(xP+2+0));
	gradPhi[2] += ref_dist[0] * ref_dist[1] * ((*phi)(xP+2+1+0) - (*phi)(xP+1+0));
	
	gradPhi[0] /= lat_resolution;
	gradPhi[1] /= lat_resolution;
	gradPhi[2] /= lat_resolution;  
	
	for (int i = 0 ; i < 3; i++)
	{
		(*part).vel[i] = -gradPhi[i] * params[0];
		v2 += (*part).vel[i] * (*part).vel[i];
	}
	
	return v2;
	
#undef phi
#undef xP
}


//////////////////////////
// loadHomogeneousTemplate
//////////////////////////
// Description:
//   loads a homogeneous template from a GADGET-2 file
// 
// Arguments:
//   filename   string containing the path to the template file
//   numpart    will contain the number of particles of the template
//   partdata   will contain the particle positions (memory will be allocated)
//
// Returns:
// 
//////////////////////////

void loadHomogeneousTemplate(const char * filename, long & numpart, float * & partdata)
{
	int i;

	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * templatefile;
		int blocksize1, blocksize2, num_read;
		gadget2_header filehdr;
		
		templatefile = fopen(filename, "r");
		
		if (templatefile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unable to open template file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		fread(&blocksize1, sizeof(int), 1, templatefile);
		if (blocksize1 != sizeof(filehdr))
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unknown template file format - header not recognized." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}		
		fread(&filehdr, sizeof(filehdr), 1, templatefile);
		fread(&blocksize2, sizeof(int), 1, templatefile);
		if (blocksize1 != blocksize2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unknown template file format - block size mismatch while reading header." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		// analyze header for compatibility
		if (filehdr.num_files != 1)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Multiple input files (" << filehdr.num_files << ") currently not supported." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		if (filehdr.BoxSize <= 0.)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! BoxSize = " << filehdr.BoxSize << " not allowed." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		if (filehdr.npart[1] <= 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! No particles declared." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		partdata = (float *) malloc(3 * sizeof(float) * filehdr.npart[1]);
		if (partdata == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Memory error." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		fread(&blocksize1, sizeof(int), 1, templatefile);
		if (filehdr.npart[0] > 0)
		{
			if (fseek(templatefile, 3 * sizeof(float) * filehdr.npart[0], SEEK_CUR))
			{
				cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unsuccesful attempt to skip gas particles (" << filehdr.npart[0] << ")." << endl;
				fclose(templatefile);
				parallel.abortForce();
			}
		}
		num_read = fread(partdata, sizeof(float), 3 * filehdr.npart[1], templatefile);
		if (num_read != 3 * filehdr.npart[1])
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unable to read particle data." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		for (i = 2; i < 6; i++)
		{
			if (filehdr.npart[i] > 0)
			{
				if (fseek(templatefile, 3 * sizeof(float) * filehdr.npart[i], SEEK_CUR))
				{
					cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unsuccesful attempt to skip particles (" << filehdr.npart[i] << ")." << endl;
					parallel.abortForce();
				}
			}
		}
		fread(&blocksize2, sizeof(int), 1, templatefile);
		if (blocksize1 != blocksize2)
			{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unknown template file format - block size mismatch while reading particles." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		fclose(templatefile);
		
		// reformat and check particle data
		for (i = 0; i < 3 * filehdr.npart[1]; i++)
		{
			partdata[i] /= filehdr.BoxSize;
			if (partdata[i] < 0. || partdata[i] > 1.)
			{
				cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Particle data corrupted." << endl;
				parallel.abortForce();
			}
		}
		numpart = (long) filehdr.npart[1];
		
		parallel.broadcast_dim0<long>(numpart, 0);
	}
	else
	{
		parallel.broadcast_dim0<long>(numpart, 0);
		
		if (numpart <= 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Communication error." << endl;
			parallel.abortForce();
		}
		
		partdata = (float *) malloc(3 * sizeof(float) * numpart);
		if (partdata == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<float>(partdata, 3 * numpart, 0);
}


//////////////////////////
// loadPowerSpectrum
//////////////////////////
// Description:
//   loads a tabulated matter power spectrum from a file
// 
// Arguments:
//   filename   string containing the path to the template file
//   pkspline   will point to the gsl_spline which holds the tabulated
//              power spectrum (memory will be allocated)
//   boxsize    comoving box size (in the same units as used in the file)
//
// Returns:
// 
//////////////////////////

void loadPowerSpectrum(const char * filename, gsl_spline * & pkspline, const double boxsize)
{
	int i = 0, numpoints = 0;
	double * k;
	double * pk;

	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * pkfile;
		char line[MAX_LINESIZE];
		double dummy1, dummy2;
		
		line[MAX_LINESIZE-1] = 0;
		
		pkfile = fopen(filename, "r");
		
		if (pkfile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Unable to open file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		while (!feof(pkfile) && !ferror(pkfile))
		{
			fgets(line, MAX_LINESIZE, pkfile);
			if (line[MAX_LINESIZE-1] != 0)
			{
				cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Character limit (" << (MAX_LINESIZE-1) << "/line) exceeded in file " << filename << "." << endl;
				fclose(pkfile);
				parallel.abortForce();
			}
			
			if (sscanf(line, " %lf %lf", &dummy1, &dummy2) == 2 && !feof(pkfile) && !ferror(pkfile)) numpoints++;
		}
		
		if (numpoints < 2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! No valid data found in file " << filename << "." << endl;
			fclose(pkfile);
			parallel.abortForce();
		}
		
		k = (double *) malloc(sizeof(double) * numpoints);
		pk = (double *) malloc(sizeof(double) * numpoints);
		
		if (k == NULL || pk == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Memory error." << endl;
			fclose(pkfile);
			parallel.abortForce();
		}
		
		rewind(pkfile);
		
		while (!feof(pkfile) && !ferror(pkfile))
		{
			fgets(line, MAX_LINESIZE, pkfile);
			
			if (sscanf(line, " %lf %lf", &dummy1, &dummy2) == 2 && !feof(pkfile) && !ferror(pkfile))
			{
				if (dummy1 < 0. || dummy2 < 0.)
				{
					cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Negative entry encountered." << endl;
					free(k);
					free(pk);
					fclose(pkfile);
					parallel.abortForce();
				}
				
				if (i > 0)
				{
					if (k[i-1] >= dummy1 * boxsize)
					{
						cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! k-values are not strictly ordered." << endl;
						free(k);
						free(pk);
						fclose(pkfile);
						parallel.abortForce();
					}
				}
				
				k[i] = dummy1 * boxsize;
				pk[i] = sqrt(0.5 * dummy2 * boxsize);
				i++;
			}
		}
		
		fclose(pkfile);
		
		if (i != numpoints)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! File may have changed or file pointer corrupted." << endl;
			free(k);
			free(pk);
			parallel.abortForce();
		}
		
		parallel.broadcast_dim0<int>(numpoints, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numpoints, 0);
		
		if (numpoints < 2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Communication error." << endl;
			parallel.abortForce();
		}
		
		k = (double *) malloc(sizeof(double) * numpoints);
		pk = (double *) malloc(sizeof(double) * numpoints);
		
		if (k == NULL || pk == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<double>(k, numpoints, 0);
	parallel.broadcast_dim0<double>(pk, numpoints, 0);
	
	pkspline = gsl_spline_alloc(gsl_interp_cspline, numpoints);
	
	gsl_spline_init(pkspline, k, pk, numpoints);
	
	free(k);
	free(pk);
}


//////////////////////////
// generateCICKernel
//////////////////////////
// Description:
//   generates convolution kernel for CIC projection
// 
// Arguments:
//   ker        reference to allocated field that will contain the convolution kernel
//   numpcl     number of particles in the pcldata pointer
//   pcldata    raw particle data from which the kernel will be constructed;
//              a standard kernel will be provided in case no particles are specified
//   numtile    tiling factor used for the particle template
//
// Returns:
// 
//////////////////////////

void generateCICKernel(Field<Real> & ker, const long numpcl = 0, float * pcldata = NULL, const int numtile = 1)
{
	const long linesize = ker.lattice().sizeLocal(0);
	Real renorm = linesize * linesize;
	long i, oct, sx, sy, sz;
	float wx, wy, wz, q1, q2, q3, q4, ww;
	Site x(ker.lattice());
	
	for (x.first(); x.test(); x.next())
	{
		ker(x) = 0.;
	}
	
	if (numpcl == 0 || pcldata == NULL) // standard kernel
	{
		if (x.setCoord(0, 0, 0))
		{
			ker(x) = 6. * renorm;
			ker(x+0) = -renorm;
			x.setCoord(linesize-1, 0, 0);
			ker(x) = -renorm;
		}
	
		if (x.setCoord(0, 1, 0))
			ker(x) = -renorm;
	
		if (x.setCoord(0, 0, 1))
			ker(x) = -renorm;
	
		if (x.setCoord(0, linesize-1, 0))
			ker(x) = -renorm;
	
		if (x.setCoord(0, 0, linesize-1))
			ker(x) = -renorm;
		
		return;
	}
	
	// compute kernel explicitly
	
	renorm /= (Real) (numpcl * (long) numtile * (long) numtile * (long) numtile) / (Real) (linesize * linesize * linesize);
	
	for (i = 0; i < numpcl; i++)
	{
		for (oct = 0; oct < 8; oct++)
		{
			if (oct == 0)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
				
				// particle is in the first octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = 1;
				sy = 1;
				sz = 1;
			}
			else if (oct == 1)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
					
				// particle is in the second octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = -1;
				sy = 1;
				sz = 1;
			}
			else if (oct == 2)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
					
				// particle is in the third octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = -1;
				sy = -1;
				sz = 1;
			}
			else if (oct == 3)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
					
				// particle is in the fourth octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = 1;
				sy = -1;
				sz = 1;
			}
			else if (oct == 4)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
				
				// particle is in the fifth octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = 1;
				sy = 1;
				sz = -1;
			}
			else if (oct == 5)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
					
				// particle is in the sixth octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = -1;
				sy = 1;
				sz = -1;
			}
			else if (oct == 6)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
					
				// particle is in the seventh octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = -1;
				sy = -1;
				sz = -1;
			}
			else if (oct == 7)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
					
				// particle is in the eight-th octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = 1;
				sy = -1;
				sz = -1;
			}
			else
				continue;
			
			// 0-direction
			
			ww = wy*wz*renorm;		   
			q1 = (wx > 0.9) ? 2. : 1.;
			
			if (x.setCoord(0, 0, 0))
			{
				ker(x) += ww * wy * wz * q1;
				x.setCoord((linesize+sx)%linesize, 0, 0);
				ker(x) -= ww * wy * wz;
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, 0, 0);
					ker(x) -= ww * wy * wz;
				}
			}
			
			if (x.setCoord(0, 0, (linesize+sz)%linesize))
			{
				ker(x) += ww * wy * (1.-wz) * q1;
				x.setCoord((linesize+sx)%linesize, 0, (linesize+sz)%linesize);
				ker(x) -= ww * wy * (1.-wz);
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, 0, (linesize+sz)%linesize);
					ker(x) -= ww * wy * (1.-wz);
				}
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, 0))
			{
				ker(x) += ww * (1.-wy) * wz * q1;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, 0);
				ker(x) -= ww * (1.-wy) * wz;
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, (linesize+sy)%linesize, 0);
					ker(x) -= ww * (1.-wy) * wz;
				}
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, (linesize+sz)%linesize))
			{
				ker(x) += ww * (1.-wy) * (1.-wz) * q1;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wy) * (1.-wz);
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
					ker(x) -= ww * (1.-wy) * (1.-wz);
				}
			}
			
			// 1-direction
			
			ww = wx*wz*renorm;
			q1 = (wy > 0.9) ? 2. : 1.;
			
			if (x.setCoord(0, 0, 0))
			{
				ker(x) += ww * wx * wz * q1;
				x.setCoord((linesize+sx)%linesize, 0, 0);
				ker(x) += ww * (1.-wx) * wz * q1;
			}
			
			if (x.setCoord(0, 0, (linesize+sz)%linesize))
			{
				ker(x) += ww * wx * (1.-wz) * q1;
				x.setCoord((linesize+sx)%linesize, 0, (linesize+sz)%linesize);
				ker(x) += ww * (1.-wx) * (1.-wz) * q1;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, 0))
			{
				ker(x) -= ww * wx * wz;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, 0);
				ker(x) -= ww * (1.-wx) * wz;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, (linesize+sz)%linesize))
			{
				ker(x) -= ww * wx * (1.-wz);
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wx) * (1.-wz);
			}
			
			if (wy > 0.9)
			{
				if (x.setCoord(0, (linesize-sy)%linesize, 0))
				{
					ker(x) -= ww * wx * wz;
					x.setCoord((linesize+sx)%linesize, (linesize-sy)%linesize, 0);
					ker(x) -= ww * (1.-wx) * wz;
				}
				
				if (x.setCoord(0, (linesize-sy)%linesize, (linesize+sz)%linesize))
				{
					ker(x) -= ww * wx * (1.-wz);
					x.setCoord((linesize+sx)%linesize, (linesize-sy)%linesize, (linesize+sz)%linesize);
					ker(x) -= ww * (1.-wx) * (1.-wz);
				}
			}
						
			// 2-direction
			
			ww = wx*wy*renorm;
			q1 = (wz > 0.9) ? 2. : 1.;
			
			if (x.setCoord(0, 0, 0))
			{
				ker(x) += ww * wx * wy * q1;
				x.setCoord((linesize+sx)%linesize, 0, 0);
				ker(x) += ww * (1.-wx) * wy * q1;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, 0))
			{
				ker(x) += ww * wx * (1.-wy) * q1;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, 0);
				ker(x) += ww * (1.-wx) * (1.-wy) * q1;
			}
			
			if (x.setCoord(0, 0, (linesize+sz)%linesize))
			{
				ker(x) -= ww * wx * wy;
				x.setCoord((linesize+sx)%linesize, 0, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wx) * wy;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, (linesize+sz)%linesize))
			{
				ker(x) -= ww * wx * (1.-wy);
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wx) * (1.-wy);
			}
			
			if (wz > 0.9)
			{
				if (x.setCoord(0, 0, (linesize-sz)%linesize))
				{
					ker(x) -= ww * wx * wy;
					x.setCoord((linesize+sx)%linesize, 0, (linesize-sz)%linesize);
					ker(x) -= ww * (1.-wx) * wy;
				}
				
				if (x.setCoord(0, (linesize+sy)%linesize, (linesize-sz)%linesize))
				{
					ker(x) -= ww * wx * (1.-wy);
					x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize-sz)%linesize);
					ker(x) -= ww * (1.-wx) * (1.-wy);
				}
			}
		}
	}
}


#ifdef FFT3D

//////////////////////////
// generateDisplacementField
//////////////////////////
// Description:
//   generates particle displacement field
// 
// Arguments:
//   potFT      reference to allocated field that contains the convolution kernel relating the potential
//              (generating the displacement field) with the bare density perturbation; will contain the
//              Fourier image of the potential generating the displacement field
//   coeff      gauge correction coefficient "H_conformal^2"
//   pkspline   pointer to a gsl_spline which holds a tabulated power spectrum
//   seed       initial seed for random number generator
//   ksphere    flag to indicate that only a sphere in k-space should be initialized
//              (default = 0: full k-space cube is initialized)
//
// Returns:
// 
//////////////////////////

void generateDisplacementField(Field<Cplx> & potFT, const Real coeff, const gsl_spline * pkspline, const unsigned int seed, const int ksphere = 0)
{
	const int linesize = potFT.lattice().size(1);
	const int kmax = (linesize / 2) - 1;
	rKSite k(potFT.lattice());
	int kx, ky, kz, i, j;
	int kymin, kymax, kzmin, kzmax;
	long jumpy, jumpz;
	float r1, r2, k2, s;
	float * sinc;
	sitmo::prng_engine prng;
	uint64_t huge_skip = HUGE_SKIP;
	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	
	sinc = (float *) malloc(linesize * sizeof(float));
	
	sinc[0] = 1.;
	for (i = 1; i < linesize; i++)
	{
		sinc[i] = sin(M_PI * (float) i / (float) linesize) * (float) linesize / (M_PI * (float) i);
	}
	
	k.initialize(potFT.lattice(), potFT.lattice().siteLast());
	kymax = k.coord(1);
	kzmax = k.coord(2);
	k.initialize(potFT.lattice(), potFT.lattice().siteFirst());
	kymin = k.coord(1);
	kzmin = k.coord(2);
		
	if (kymin < (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
	{
		prng.seed(seed);
		   
		if (kymin == 0 && kzmin == 0)
		{
			k.setCoord(0, 0, 0);
			potFT(k) = Cplx(0.,0.);
			kx = 1;
		}
		else
		{
			kx = 0;
			prng.discard(((uint64_t) kzmin * huge_skip + (uint64_t) kymin) * huge_skip); 
		}
		
		for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				for (i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) (ky * ky) + (float) (kz * kz);
					
					if (kx >= kmax || ky >= kmax || kz >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[ky] * sinc[kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s * (1. + 7.5 * coeff / k2) / potFT(k);
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
				kx = 0;
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymax >= (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + (uint64_t) kzmin) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		
		for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				for (kx = 0, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
										
					k2 = (float) (kx * kx) + (float) ((linesize-ky) * (linesize-ky)) + (float) (kz * kz);
					
					if (kx >= kmax || (linesize-ky) >= kmax || kz >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[linesize-ky] * sinc[kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
						
						potFT(k) = Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s * (1. + 7.5 * coeff / k2) / potFT(k);
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymin < (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) kymin) * huge_skip);
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) (ky * ky) + (float) ((linesize-kz) * (linesize-kz));
					
					if (kx >= kmax || ky >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[ky] * sinc[linesize-kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s * (1. + 7.5 * coeff / k2) / potFT(k);
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
		
		prng.seed(seed);
		prng.discard(((uint64_t) (linesize - kzmax) * huge_skip + (uint64_t) kymin) * huge_skip);
		kx = 0;
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				k.setCoord(kx, ky, kz);
					
				k2 = (float) (ky * ky) + (float) ((linesize-kz) * (linesize-kz));
				i = 0;
				
				if (ky >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
				{
					potFT(k) = Cplx(0., 0.);
				}
				else
				{
					s = sinc[ky] * sinc[linesize-kz];
					k2 *= 4. * M_PI * M_PI;
					do
					{
						r1 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					}
					while (r1 == 0.);
					r2 = (float) prng() / (float) sitmo::prng_engine::max();
					i++;
				
					potFT(k) = Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s * (1. + 7.5 * coeff / k2) / potFT(k);
				}
								
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymax >= (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) ((linesize-ky) * (linesize-ky)) + (float) ((linesize-kz) * (linesize-kz));
					
					if (kx >= kmax || (linesize-ky) >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[linesize-ky] * sinc[linesize-kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s * (1. + 7.5 * coeff / k2) / potFT(k);
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
		
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		kx = 0;
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				k.setCoord(kx, ky, kz);
					
				k2 = (float) ((linesize-ky) * (linesize-ky)) + (float) ((linesize-kz) * (linesize-kz));
				i = 0;
				
				if ((linesize-ky) >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
				{
					potFT(k) = Cplx(0., 0.);
				}
				else
				{
					s = sinc[linesize-ky] * sinc[linesize-kz];
					k2 *= 4. * M_PI * M_PI;
					do
					{
						r1 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					}
					while (r1 == 0.);
					r2 = (float) prng() / (float) sitmo::prng_engine::max();
					i++;
				
					potFT(k) = Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s * (1. + 7.5 * coeff / k2) / potFT(k);
				}
				
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	gsl_interp_accel_free(acc);
	free(sinc);
}
#endif


//////////////////////////
// initializeParticlePositions
//////////////////////////
// Description:
//   initializes particle positions using a homogeneous template and a displacement field
// 
// Arguments:
//   pot        reference to the potential generating the displacement field
//   numpart    number of particles of the template
//   partdata   particle positions in the template
//   numtile    tiling factor for homogeneous template - total particle number will be
//              numpart * numtile^3
//   pcls       reference to (empty) particle object which will contain the new particle ensemble
//   coeff      scaling coefficient to be applied to the displacement field
//
// Returns:
//   maximum displacement
// 
//////////////////////////

Real initializeParticlePositions(Field<Real> & pot, const long numpart, const float * partdata, const int numtile, Particles<part_simple,part_simple_info,part_simple_dataType> & pcls, const float coeff)
{
	long xtile, ytile, ztile, i;
	Real x, dx, fx, x0, y, dy, fy, y0, z, dz, fz, z0;
	Real dmax = 0.;
	Site p(pot.lattice());
	
	part_simple part;
	
	for (ztile = (pot.lattice().coordSkip()[0] * numtile) / pot.lattice().size(2); ztile <= ((pot.lattice().coordSkip()[0] + pot.lattice().sizeLocal(2)) * numtile) / pot.lattice().size(2); ztile++)
	{
		if (ztile >= numtile) break;
		
		for (ytile = (pot.lattice().coordSkip()[1] * numtile) / pot.lattice().size(1); ytile <= ((pot.lattice().coordSkip()[1] + pot.lattice().sizeLocal(1)) * numtile) / pot.lattice().size(1); ytile++)
		{
			if (ytile >= numtile) break;
			
			for (xtile = 0; xtile < numtile; xtile++)
			{
				for (i = 0; i < numpart; i++)
				{
					x = ((Real) xtile + partdata[3*i]) / (Real) numtile;
					y = ((Real) ytile + partdata[3*i+1]) / (Real) numtile;
					z = ((Real) ztile + partdata[3*i+2]) / (Real) numtile;
					
					fx = modf(x * pot.lattice().size(0), &x0);
					fy = modf(y * pot.lattice().size(1), &y0);
					fz = modf(z * pot.lattice().size(2), &z0);
					
					if (!p.setCoord((int) x0, (int) y0, (int) z0)) continue;
					
					dx = (1.-fy) * (1.-fz) * (pot(p+0) - pot(p));
					dy = (1.-fx) * (1.-fz) * (pot(p+1) - pot(p));
					dz = (1.-fx) * (1.-fy) * (pot(p+2) - pot(p));
					dx += fy * (1.-fz) * (pot(p+1+0) - pot(p+1));
					dy += fx * (1.-fz) * (pot(p+1+0) - pot(p+0));
					dz += fx * (1.-fy) * (pot(p+2+0) - pot(p+0));
					dx += (1.-fy) * fz * (pot(p+2+0) - pot(p+2));
					dy += (1.-fx) * fz * (pot(p+2+1) - pot(p+2));
					dz += (1.-fx) * fy * (pot(p+2+1) - pot(p+1));
					dx += fy * fz * (pot(p+2+1+0) - pot(p+2+1));
					dy += fx * fz * (pot(p+2+1+0) - pot(p+2+0));
					dz += fx * fy * (pot(p+2+1+0) - pot(p+1+0));
					
					part.ID = i + numpart * (xtile + (long) numtile * (ytile + (long) numtile * ztile));
					part.pos[0] = x;
					part.pos[1] = y;
					part.pos[2] = z;
					part.vel[0] = dx * pot.lattice().size(0) / coeff;
					part.vel[1] = dy * pot.lattice().size(0) / coeff;
					part.vel[2] = dz * pot.lattice().size(0) / coeff;
					pcls.addParticle_global(part);
					
					if (dmax < dx * dx + dy * dy + dz * dz) dmax = dx * dx + dy * dy + dz * dz;
				}
			}
		}
	}
	
	pcls.moveParticles(move_particles_ic_basic, 1./coeff);
	
	return sqrt(dmax) * pot.lattice().size(0) / (coeff * coeff);
}


//////////////////////////
// applyMomentumDistribution
//////////////////////////
// Description:
//   adds a random momentum vector drawn from a Fermi-Dirac distribution to
//   each particle. The current implementation uses a simple rejection-sampling
//   method based on the ziggurat algorithm [G. Marsaglia and W.W. Tsang,
//   J. Stat. Softw. 5 (2000) 1]. The "ziggurat" is hard-coded and was precomputed
//   for the ultrarelativistic limit of a Fermi-Dirac distribution with zero
//   chemical potential. A method to construct the ziggurat on-the-fly for more
//   general distribution functions could be implemented in the future.
// 
// Arguments:
//   pcls   pointer to particle handler
//   seed   seed for random number generator
//   T_m    dimensionless parameter in the distribution function
//          in most cases the ratio of temperature and fundamental mass
//
// Returns: sum of momenta over all particles (for reduction)
// 
//////////////////////////

double applyMomentumDistribution(Particles<part_simple,part_simple_info,part_simple_dataType> * pcls, unsigned int seed, float T_m = 0.)
{	
	Site xPart(pcls->lattice());	
	typename std::list<part_simple>::iterator it;
	sitmo::prng_engine prng;
	float r1, r2, q;
	uint32_t i, r;
	double sum_q = 0.0;
	
	const float ql[] = {0.0f,       0.0453329523f, 0.0851601009f, 0.115766097f,
	                    0.142169202f, 0.166069623f, 0.188283033f, 0.209275639f,
	                    0.229344099f, 0.248691555f, 0.267464827f, 0.285774552f,
	                    0.303706948f, 0.321331123f, 0.338703809f, 0.355872580f,
	                    0.372878086f, 0.389755674f, 0.406536572f, 0.423248791f,
	                    0.439917819f, 0.456567164f, 0.473218795f, 0.489893502f,
	                    0.506611198f, 0.523391180f, 0.540252353f, 0.557213439f,
	                    0.574293166f, 0.591510448f, 0.608884565f, 0.626435339f,
	                    0.644183319f, 0.662149971f, 0.680357893f, 0.698831041f,
	                    0.717594982f, 0.736677192f, 0.756107381f, 0.775917886f,
	                    0.796144127f, 0.816825143f, 0.838004248f, 0.859729814f,
	                    0.882056241f, 0.905045149f, 0.928766878f, 0.953302387f,
	                    0.978745698f, 1.00520708f,  1.03281729f,  1.06173322f,
	                    1.09214584f,  1.12429121f,  1.15846661f,  1.19505456f,
	                    1.23456031f,  1.27767280f,  1.32536981f,  1.37911302f,
	                    1.44124650f,  1.51592808f,  1.61180199f,  1.75307820f, 29.0f};
	                    
	const float qr[] = {29.0f,       11.8477879f, 10.3339062f, 9.58550750f,
	                    9.08034038f, 8.69584518f, 8.38367575f, 8.11975908f,
	                    7.89035001f, 7.68685171f, 7.50352431f, 7.33634187f,
	                    7.18236952f, 7.03940031f, 6.90573157f, 6.78002122f,
	                    6.66119176f, 6.54836430f, 6.44081188f, 6.33792577f,
	                    6.23919052f, 6.14416529f, 6.05246957f, 5.96377208f,
	                    5.87778212f, 5.79424263f, 5.71292467f, 5.63362289f,
	                    5.55615183f, 5.48034280f, 5.40604138f, 5.33310512f,
	                    5.26140176f, 5.19080749f, 5.12120558f, 5.05248501f,
	                    4.98453932f, 4.91726547f, 4.85056276f, 4.78433177f,
	                    4.71847331f, 4.65288728f, 4.58747144f, 4.52212011f,
	                    4.45672261f, 4.39116154f, 4.32531058f, 4.25903193f,
	                    4.19217309f, 4.12456265f, 4.05600493f, 3.98627278f,
	                    3.91509779f, 3.84215661f, 3.76705129f, 3.68928014f,
	                    3.60819270f, 3.52291720f, 3.43223655f, 3.33436084f,
	                    3.22646839f, 3.10364876f, 2.95592669f, 2.75624893f, 0.0f};
	                    
	const float f[] = {0.0f,          0.0010042516f, 0.0034718142f, 0.00631345904f,
	                  0.00938886471f, 0.0126471707f, 0.0160614805f, 0.0196150986f,
	                   0.0232967063f, 0.0270982040f, 0.0310135938f, 0.0350383391f,
	                   0.0391689707f, 0.0434028310f, 0.0477379005f, 0.0521726764f,
	                   0.0567060859f, 0.0613374223f, 0.0660662983f, 0.0708926105f,
	                   0.0758165136f, 0.0808384013f, 0.0859588926f, 0.0911788222f,
	                   0.0964992355f, 0.101921386f,  0.107446735f,  0.113076958f,
	                   0.118813945f,  0.124659815f,  0.130616922f,  0.136687871f,
	                   0.142875536f,  0.149183078f,  0.155613967f,  0.162172016f,
	                   0.168861408f,  0.175686736f,  0.182653051f,  0.189765914f,
	                   0.197031455f,  0.204456455f,  0.212048433f,  0.219815749f,
	                   0.227767740f,  0.235914877f,  0.244268957f,  0.252843349f,
	                   0.261653294f,  0.270716296f,  0.280052614f,  0.289685921f,
	                   0.299644172f,  0.309960782f,  0.320676286f,  0.331840691f,
	                   0.343516979f,  0.355786485f,  0.368757588f,  0.382580625f,
	                   0.397475563f,  0.413789108f,  0.432131941f,  0.453799050f, 0.482830296f};
	                    
	for (xPart.first(); xPart.test(); xPart.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for (it=(pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				prng.seed(seed);
				prng.discard((uint64_t) (7l * (*it).ID));
				
				while (true)
				{
					r = prng();
					i = r % 64;
					r /= 64;
					
					q = ql[i] + 64.0f * (qr[i]-ql[i]) * ((float) r / (float) sitmo::prng_engine::max());
					
					if (q > ql[i+1] && q < qr[i+1]) break;
					
					if (f[i] + (f[i+1]-f[i]) * ((float) prng() / (float) sitmo::prng_engine::max()) < q * q / (exp(q) + 1.0f)) break;
				}
				
				r1 = acos(2. * ((float) prng() / (float) sitmo::prng_engine::max()) - 1.);
				r2 = 2 * M_PI * ((float) prng() / (float) sitmo::prng_engine::max());
				
				q *= T_m;
				
				(*it).vel[0] += cos(r2) * sin(r1) * q;
				(*it).vel[1] += sin(r2) * sin(r1) * q;
				(*it).vel[2] += cos(r1) * q;
				
				sum_q += q;
			}
		}
	}
	
	return sum_q;
}


#ifdef FFT3D

//////////////////////////
// generateIC_basic
//////////////////////////
// Description:
//   basic initial condition generator
// 
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   dtau           time step
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   max_vel_ncdm   array that will contain the maximum q/m/a (max. velocity)
//                  of non-cold species
//   phi            pointer to allocated field
//   psi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_psi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//
// Returns: maximum velocity of DM within domain
// 
//////////////////////////

double generateIC_basic(metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, const double dtau, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, double * max_vel_ncdm, Field<Real> * phi, Field<Real> * psi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_psi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij)
{
	int i, j, p;
	double a = 1. / (1. + sim.z_in);
	float * pcldata = NULL;
	float * ncdm_pcldata = NULL;
	gsl_spline * pkspline = NULL;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
	double max_displacement;
	double max_vel;
	double rescale;
	double mean_q;
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_ncdm_info[MAX_PCL_SPECIES];
	part_simple_dataType pcls_ncdm_dataType;
	Real boxSize[3] = {1.,1.,1.};
	
	loadHomogeneousTemplate(ic.pclfile, sim.numpcl[0], pcldata);
	
	if (pcldata == NULL)
	{
		COUT << " error: particle data was empty!" << endl;
		parallel.abortForce();
	}
	
	loadPowerSpectrum(ic.pkfile, pkspline, sim.boxsize);
	
	if (pkspline == NULL)
	{
		COUT << " error: power spectrum was empty!" << endl;
		parallel.abortForce();
	}
	
	if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
	{
		generateCICKernel(*psi, sim.numpcl[0], pcldata, ic.numtile);
	}
	else
		generateCICKernel(*psi);	
	
	plan_psi->execute(FFT_FORWARD);
	generateDisplacementField(*scalarFT, sim.gr_flag * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo), pkspline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);	
	plan_psi->execute(FFT_BACKWARD);
	
	for (x.first(); x.test(); x.next())
	{
		(*psi)(x) *= Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) / sim.boxsize / sim.boxsize;
	}
	
	psi->updateHalo();
	
	strcpy(pcls_cdm_info.type_name, "part_simple");
	pcls_cdm_info.mass = (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) (sim.numpcl[0]*(long)ic.numtile*(long)ic.numtile*(long)ic.numtile);
	pcls_cdm_info.relativistic = false;
	
	pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);
	
	max_displacement = initializeParticlePositions(*psi, sim.numpcl[0], pcldata, ic.numtile, *pcls_cdm, Hconf(a, fourpiG, cosmo));
	parallel.max(max_displacement);
	
	COUT << " " << sim.numpcl[0] * (long) ic.numtile * (long) ic.numtile * (long) ic.numtile << " cdm particles initialized: maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
	
	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		ncdm_pcldata = (float *) malloc(3 * sizeof(float) * sim.numpcl[0] * ic.psd_samples[p]);
	
		if (ncdm_pcldata == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in generateIC_simple! Memory error." << endl;
			parallel.abortForce();
		}
	
		for (i = sim.numpcl[0]-1; i >= 0; i--)
		{
			for (j = 0; j < ic.psd_samples[p]; j++)
			{
				ncdm_pcldata[3*(ic.psd_samples[p]*i+j)] = pcldata[3*i];
				ncdm_pcldata[3*(ic.psd_samples[p]*i+j)+1] = pcldata[3*i+1];
				ncdm_pcldata[3*(ic.psd_samples[p]*i+j)+2] = pcldata[3*i+2];
			}
		}
		sim.numpcl[p+1] = sim.numpcl[0] * ic.psd_samples[p];
		
		strcpy(pcls_ncdm_info[p].type_name, "part_simple");
		pcls_ncdm_info[p].mass = cosmo.Omega_ncdm[p] / (Real) (sim.numpcl[p+1]*(long)ic.numtile*(long)ic.numtile*(long)ic.numtile);
		pcls_ncdm_info[p].relativistic = true;
		
		pcls_ncdm[p].initialize(pcls_ncdm_info[p], pcls_ncdm_dataType, &(phi->lattice()), boxSize);
		
		max_displacement = initializeParticlePositions(*psi, sim.numpcl[p+1], ncdm_pcldata, ic.numtile, pcls_ncdm[p], Hconf(a, fourpiG, cosmo));
		parallel.max(max_displacement);
		
		sim.numpcl[p+1] *= (long) ic.numtile * (long) ic.numtile * (long) ic.numtile;
	
		COUT << " " << sim.numpcl[p+1] << " ncdm particles initialized for species " << p+1 << ": maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
		
		free(ncdm_pcldata);
	}
	
	free(pcldata);
	gsl_spline_free(pkspline);
	
	sim.numpcl[0] *= (long) ic.numtile * (long) ic.numtile * (long) ic.numtile;
	
	projection_init(source);
	scalarProjectionCIC_project(pcls_cdm, source);
	for (p = 0; p < cosmo.num_ncdm; p++)
		scalarProjectionCIC_project(pcls_ncdm+p, source);
	scalarProjectionCIC_comm(source);
	
	plan_source->execute(FFT_FORWARD);
	
	kFT.first();
	if (kFT.coord(0) == 0 && kFT.coord(1) == 0 && kFT.coord(2) == 0)
		(*scalarFT)(kFT) = Cplx(0.,0.);
				
	solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a, 3. * sim.gr_flag * (Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) + fourpiG * cosmo.Omega_m / a));
	
	plan_phi->execute(FFT_BACKWARD);

	for (x.first(); x.test(); x.next())
	{
		(*psi)(x) = (*phi)(x);
		(*chi)(x) = 0.;
	}
	
	psi->updateHalo();
	phi->updateHalo();
	chi->updateHalo();
	
	rungekutta4bg(a, fourpiG, cosmo, -0.5*dtau);
	rescale = a / Hconf(a, fourpiG, cosmo) / (1.5 * Omega_m(a, cosmo) + 2. * Omega_rad(a, cosmo));
	max_vel = pcls_cdm->updateVel(initialize_q_ic_basic, 0., &psi, 1, &rescale) / a;
			
	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		rescale = a / Hconf(a, fourpiG, cosmo) / (1.5 * Omega_m(a, cosmo) + 2. * Omega_rad(a, cosmo));
		pcls_ncdm[p].updateVel(initialize_q_ic_basic, 0., &psi, 1, &rescale);
		if (cosmo.m_ncdm[p] > 0.)
		{
			rescale = pow(cosmo.Omega_g * cosmo.h * cosmo.h / 4.48147e-7, 0.25) * cosmo.T_ncdm[p] * 8.61733e-5 / cosmo.m_ncdm[p];
			mean_q = applyMomentumDistribution(pcls_ncdm+p, (unsigned int) (ic.seed + p), rescale);
			parallel.sum(mean_q);
			COUT << " species " << p+1 << " Fermi-Dirac distribution had mean q/m = " << mean_q / sim.numpcl[p+1] << endl;
		}
		max_vel_ncdm[p] = pcls_ncdm[p].updateVel(update_q, 0., &psi, 1, &a);
	}
	a = 1. / (1. + sim.z_in);
	
	projection_init(Bi);
	vectorProjectionCICNGP_project(pcls_cdm, Bi);
	for (p = 0; p < cosmo.num_ncdm; p++)
		vectorProjectionCICNGP_project(pcls_ncdm+p, Bi);
	vectorProjectionCICNGP_comm(Bi);
	plan_Bi->execute(FFT_FORWARD);
	projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);	
	plan_Bi->execute(FFT_BACKWARD);	
	Bi->updateHalo();
	
	projection_init(Sij);
	projection_Tij_project(pcls_cdm, Sij, a, phi);
	for (p = 0; p < cosmo.num_ncdm; p++)
		projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
	projection_Tij_comm(Sij);
	
	prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);	
	plan_Sij->execute(FFT_FORWARD);	
	projectFTscalar(*SijFT, *scalarFT);
	plan_chi->execute(FFT_BACKWARD);		
	chi->updateHalo();
	
	return max_vel;
}

#endif

#endif
