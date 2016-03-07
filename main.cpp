//////////////////////////
// Copyright (c) 2015-2016 Julian Adamek (Université de Genève)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//  
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//  
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////

//////////////////////////
// main.cpp
//////////////////////////
// 
// main control sequence of Geneva N-body code with evolution of metric perturbations (gevolution)
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: March 2016
//
//////////////////////////

#include <stdlib.h>
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "parser.hpp"
#include "tools.hpp"

#ifdef EXTERNAL_IO
#ifndef NUMBER_OF_IO_FILES
#define NUMBER_OF_IO_FILES 4
#endif
#endif

using namespace std;

using namespace LATfield2;

int main(int argc, char **argv)
{
	
#ifdef BENCHMARK
	//benchmarking variables
	
	double start_time, ref_time, ref2_time, cycle_start_time;
	
	double initialization_time;
	double run_time;
	double cycle_time=0;
	double projection_time = 0;
	double snapshot_output_time = 0;
	double spectra_output_time = 0;
	double gravity_solver_time = 0;
	double fft_time = 0;
	int fft_count = 0;   
	double update_q_time = 0;
	int update_q_count = 0;
	double moveParts_time = 0;
	int  moveParts_count =0;
	
#endif  //BENCHMARK
	
	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;
	
	int i, j, cycle = 0, snapcount = 0, pkcount = 0, usedparams, numparam = 0, doneTT, redo_psi = 0, numsteps;
	int numsteps_ncdm[MAX_PCL_SPECIES];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dx, tau, a, fourpiG, tau_Lambda, maxvel;
	double maxvel_ncdm[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	string h5filename;
	char * settingsfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	gadget2_header hdr;
	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;
	int numbins;
	Real divB, curlB, divh, traceh, normh, T00hom;
	Cplx tempk;
	
	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
				break;
			case 'i':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_size =  atoi(argv[++i]);
				break;
			case 'g':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_group_size = atoi(argv[++i]);
		}
	}

#ifndef EXTERNAL_IO
	parallel.initialize(n,m);
#else
	parallel.initialize(n,m,io_size,io_group_size);
	if(parallel.isIO()) ioserver.start();
	else
	{
#endif
		
	COUT << "  _   _      _         __ ,  _" << endl;
	COUT << " (_| (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 1.0    running on " << n*m << " cores." << endl;
	COUT << "  -'" << endl << endl;
	
	if (settingsfile == NULL)
	{
		COUT << " error: no settings file specified!" << endl;
		parallel.abortForce();
	}
	
	COUT << " initializing..." << endl;
	
#ifdef BENCHMARK
	start_time = MPI_Wtime();
#endif
	
	
	numparam = loadParameterFile(settingsfile, params);
	
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);
	
	COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;
	
	sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
	saveParameterFile(filename, params, numparam);
	
	free(params);
	
	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);
	h5filename += sim.basename_snapshot;
	
	box[0] = sim.numpts;
	box[1] = sim.numpts;
	box[2] = sim.numpts;
	
	Lattice lat(3,box,1);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);
	
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES];
	Field<Real> * update_cdm_fields[3];
	Field<Real> * update_ncdm_fields[3];
	double f_params[5];

	Field<Real> phi;
	Field<Real> psi;
	Field<Real> source;
	Field<Real> chi;
	Field<Real> Sij;
	Field<Real> Bi;
	Field<Cplx> scalarFT;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;
	source.initialize(lat,1);
	psi.initialize(lat,1);
	phi.initialize(lat,1);
	chi.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_source(&source, &scalarFT);
	PlanFFT<Cplx> plan_psi(&psi, &scalarFT);
	PlanFFT<Cplx> plan_phi(&phi, &scalarFT);
	PlanFFT<Cplx> plan_chi(&chi, &scalarFT);
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	PlanFFT<Cplx> plan_Sij(&Sij, &SijFT);
	Bi.initialize(lat,3);
	BiFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi(&Bi, &BiFT);
#ifdef CHECK_B
	Field<Real> Bi_check;
	Field<Cplx> BiFT_check;
	Bi_check.initialize(lat,3);
	BiFT_check.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi_check(&Bi_check, &BiFT_check);
#endif

	update_cdm_fields[0] = &psi;
	update_cdm_fields[1] = &phi;
	update_cdm_fields[2] = &Bi;
	
	update_ncdm_fields[0] = &psi;
	update_ncdm_fields[1] = &phi;
	update_ncdm_fields[2] = &Bi;
	
	Site x(lat);
	rKSite kFT(latFT);
	
	dx = 1.0 / (double) sim.numpts;
	numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
	
	for (i = 0; i < 3; i++) // particles may never move farther than to the adjacent domain
	{
		if (lat.sizeLocal(i)-1 < sim.movelimit)
			sim.movelimit = lat.sizeLocal(i)-1;
	}
	parallel.min(sim.movelimit);

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / 2997.92458 / 2997.92458;
	a = 1. / (1. + sim.z_in);
	tau = 2. / Hconf(a, fourpiG, cosmo);
	tau_Lambda = -1.0;
	
	if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
		dtau = sim.Cf * dx;
	else
		dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
		
	dtau_old = dtau;
	
	if (ic.generator == ICGEN_BASIC)
		maxvel = generateIC_basic(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel_ncdm, &phi, &psi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_psi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij); // generates ICs on the fly
#ifdef ICGEN_PREVOLUTION
	else if (ic.generator == ICGEN_PREVOLUTION)
	{
		maxvel = generateIC_prevolution(sim, ic, cosmo, fourpiG, a, &pcls_cdm, pcls_ncdm, maxvel_ncdm, &phi, &psi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_psi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
	}
#endif
#ifdef ICGEN_FALCONIC
	else if (ic.generator == ICGEN_FALCONIC)
	{
		maxvel = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel_ncdm, &phi, &psi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_psi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
	}
#endif
	else
	{
		COUT << " error: IC generator not implemented!" << endl;
		parallel.abortForce();
	}
	
	parallel.max<double>(maxvel);	
	if (cosmo.num_ncdm > 0) parallel.max<double>(maxvel_ncdm, cosmo.num_ncdm);
	
	if (sim.gr_flag > 0)
	{
		maxvel /= sqrt(maxvel * maxvel + 1.0);
		for (i = 0; i < cosmo.num_ncdm; i++)
			maxvel_ncdm[i] /= sqrt(maxvel_ncdm[i] * maxvel_ncdm[i] + 1.0);
	}
	
	kbin = (Real *) malloc(sim.numbins * sizeof(Real));
	power = (Real *) malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	occupation = (int *) malloc(sim.numbins * sizeof(int));
	
	for (i = 0; i < 6; i++)
	{
		hdr.npart[i] = 0;
		hdr.npartTotal[i] = 0;
		hdr.mass[i] = 0.;
	}
	hdr.num_files = 1;
	hdr.Omega0 = cosmo.Omega_m;
	hdr.OmegaLambda = 1. - cosmo.Omega_m;
	hdr.HubbleParam = cosmo.h;
	hdr.BoxSize = sim.boxsize * 1000.;
	hdr.flag_sfr = 0;
	hdr.flag_cooling = 0;
	hdr.flag_feedback = 0;
	for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; i++)
		hdr.fill[i] = 0;
	
#ifdef BENCHMARK
	initialization_time = MPI_Wtime() - start_time;
	parallel.sum(initialization_time);
	COUT << " initialization complete. BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
	COUT << " initialization complete." << endl << endl;
#endif
	
	while (true)    // main loop
	{
#ifdef BENCHMARK		
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor
		projection_init(&source);
		if (sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_T00_project(pcls_ncdm+i, &source, a, &phi);
		}
		else
		{
			scalarProjectionCIC_project(&pcls_cdm, &source);
			for (i = 0; i < cosmo.num_ncdm; i++)
				scalarProjectionCIC_project(pcls_ncdm+i, &source);
		}
		projection_T00_comm(&source);
		
		projection_init(&Sij);
		projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
			projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
		projection_Tij_comm(&Sij);
		
		doneTT = 0;
		
#ifdef BENCHMARK 
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif
		// snapshot output
		if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			sprintf(filename, "%03d", snapcount);
			snapcount++;
			
			COUT << " writing snapshot at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			
#ifdef EXTERNAL_IO
			while (ioserver.openOstream()== OSTREAM_FAIL);
			
			if (sim.out_snapshot & MASK_PCLS)
				pcls_cdm.saveHDF5_server_open(h5filename + filename + "_part");
				
			if (sim.out_snapshot & MASK_DBARE && sim.gr_flag > 0)
				psi.saveHDF5_server_open(h5filename + filename + "_deltaN");
			
			if (sim.out_snapshot & MASK_T00)
				source.saveHDF5_server_open(h5filename + filename + "_T00");
				
			if (sim.out_snapshot & MASK_B)
				Bi.saveHDF5_server_open(h5filename + filename + "_B");
			
			if (sim.out_snapshot & MASK_PHI)
				phi.saveHDF5_server_open(h5filename + filename + "_phi");
				
			if (sim.out_snapshot & MASK_CHI)
				chi.saveHDF5_server_open(h5filename + filename + "_chi");
			
			if (sim.out_snapshot & MASK_HIJ)
				Sij.saveHDF5_server_open(h5filename + filename + "_hij");
				
#ifdef CHECK_B
			if (sim.out_snapshot & MASK_B)
				Bi_check.saveHDF5_server_open(h5filename + filename + "_B_check");
#endif
#endif		
			
			if (sim.out_snapshot & MASK_DBARE)
			{
				if (sim.gr_flag > 0)
				{   
					projection_init(&psi);
					scalarProjectionCIC_project(&pcls_cdm, &psi);
					for (i = 0; i < cosmo.num_ncdm; i++)
						scalarProjectionCIC_project(pcls_ncdm+i, &psi);
					scalarProjectionCIC_comm(&psi);
#ifdef EXTERNAL_IO
					psi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
					psi.saveHDF5(h5filename + filename + "_deltaN.h5");
#endif
					redo_psi = 1;
				}
				else
					source.saveHDF5(h5filename + filename + "_deltaN.h5");
			}
			
			if (sim.out_snapshot & MASK_POT)
			{
				if (!redo_psi)
				{
					projection_init(&psi);
					scalarProjectionCIC_project(&pcls_cdm, &psi);
					for (i = 0; i < cosmo.num_ncdm; i++)
						scalarProjectionCIC_project(pcls_ncdm+i, &psi);
					scalarProjectionCIC_comm(&psi);
					redo_psi = 1;
				}
				plan_psi.execute(FFT_FORWARD);				
				solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);
				plan_psi.execute(FFT_BACKWARD);
				psi.saveHDF5(h5filename + filename + "_psiN.h5");
			}
				
			if (sim.out_snapshot & MASK_T00)
#ifdef EXTERNAL_IO
				source.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
				source.saveHDF5(h5filename + filename + "_T00.h5");
#endif
				
			if (sim.out_snapshot & MASK_B)
			{
				if (sim.gr_flag == 0)
				{
					plan_Bi.execute(FFT_BACKWARD);
				}
				for (x.first(); x.test(); x.next())
				{
					Bi(x,0) /= a * a * sim.numpts;
					Bi(x,1) /= a * a * sim.numpts;
					Bi(x,2) /= a * a * sim.numpts;
				}
				Bi.updateHalo();
				
				computeVectorDiagnostics(Bi, divB, curlB);			
				COUT << " B diagnostics: max |divB| = " << divB << ", max |curlB| = " << curlB << endl;

#ifdef EXTERNAL_IO
				Bi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else				
				Bi.saveHDF5(h5filename + filename + "_B.h5");
#endif
				
				if (sim.gr_flag > 0)
				{
					plan_Bi.execute(FFT_BACKWARD);
					Bi.updateHalo();
				}
			}
			
			if (sim.out_snapshot & MASK_PHI)
#ifdef EXTERNAL_IO
				phi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
				phi.saveHDF5(h5filename + filename + "_phi.h5");
#endif
				
			if (sim.out_snapshot & MASK_CHI)
#ifdef EXTERNAL_IO
				chi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else	
				chi.saveHDF5(h5filename + filename + "_chi.h5");
#endif
				
			if (sim.out_snapshot & MASK_TIJ)	 
				Sij.saveHDF5(h5filename + filename + "_Tij.h5");
				
			if (sim.out_snapshot & MASK_HIJ)
			{
				projectFTtensor(SijFT, SijFT);
				doneTT = 1;
				plan_Sij.execute(FFT_BACKWARD);
				Sij.updateHalo();
				
				computeTensorDiagnostics(Sij, divh, traceh, normh);
				COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
				Sij.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else				
				Sij.saveHDF5(h5filename + filename + "_hij.h5");
#endif
				
				projection_init(&Sij);
				projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
				for (i = 0; i < cosmo.num_ncdm; i++)
					projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
				projection_Tij_comm(&Sij);
			}
			
			if (sim.out_snapshot & MASK_P)
			{
#ifdef CHECK_B
				projection_init(&Bi_check);
				vectorProjectionCICNGP_project(&pcls_cdm, &Bi_check);
				for (i = 0; i < cosmo.num_ncdm; i++)
					vectorProjectionCICNGP_project(pcls_ncdm+i, &Bi_check);
				vectorProjectionCICNGP_comm(&Bi_check);
				Bi_check.saveHDF5(h5filename + filename + "_p.h5");
#else
				projection_init(&Bi);
				vectorProjectionCICNGP_project(&pcls_cdm, &Bi);
				for (i = 0; i < cosmo.num_ncdm; i++)
					vectorProjectionCICNGP_project(pcls_ncdm+i, &Bi);
				vectorProjectionCICNGP_comm(&Bi);
				Bi.saveHDF5(h5filename + filename + "_p.h5");
				if (sim.gr_flag > 0)
				{
					plan_Bi.execute(FFT_BACKWARD);
					Bi.updateHalo();
				}
#endif
			}
				
#ifdef CHECK_B
			if (sim.out_snapshot & MASK_B)
			{
				if (!(sim.out_snapshot & MASK_P))
				{
					projection_init(&Bi_check);
					vectorProjectionCICNGP_project(&pcls_cdm, &Bi_check);
					for (i = 0; i < cosmo.num_ncdm; i++)
						vectorProjectionCICNGP_project(pcls_ncdm+i, &Bi_check);
					vectorProjectionCICNGP_comm(&Bi_check);
				}
				plan_Bi_check.execute(FFT_FORWARD);
				projectFTvector(BiFT_check, BiFT_check, fourpiG * dx * dx);
				plan_Bi_check.execute(FFT_BACKWARD);
			
				for (x.first(); x.test(); x.next())
				{
					Bi_check(x,0) /= a * a * sim.numpts;
					Bi_check(x,1) /= a * a * sim.numpts;
					Bi_check(x,2) /= a * a * sim.numpts;
				}
#ifdef EXTERNAL_IO
				Bi_check.saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
				Bi_check.saveHDF5(h5filename + filename + "_B_check.h5");
#endif
			}
#endif

			if (sim.out_snapshot & MASK_GADGET)
			{
				hdr.time = a;
				hdr.redshift = (1./a) - 1.;
				
				hdr.npart[1] = (unsigned int) (sim.numpcl[0] / sim.tracer_factor[0]);
				hdr.npartTotal[1] = hdr.npart[1];
				hdr.mass[1] = (double) sim.tracer_factor[0] * 27.7459457 * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0];
				pcls_cdm.saveGadget2(h5filename + filename + "_part.dat", hdr, sim.tracer_factor[0]);
				for (i = 0; i < cosmo.num_ncdm; i++)
				{
					sprintf(buffer, "_ncdm%d.dat", i);
					hdr.npart[1] = (unsigned int) (sim.numpcl[i+1] / sim.tracer_factor[i+1]);
					hdr.npartTotal[1] = hdr.npart[1];
					hdr.mass[1] = (double) sim.tracer_factor[i+1] * 27.7459457 * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1];
					pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i+1]);
				}
			}
			
			if (sim.out_snapshot & MASK_PCLS)
			{
#ifdef EXTERNAL_IO
				pcls_cdm.saveHDF5_server_write();
#else
				pcls_cdm.saveHDF5(h5filename + filename + "_part", 1);
#endif
			}
			
#ifdef EXTERNAL_IO
			ioserver.closeOstream();
#endif	
		}   // snapshot output done
		
#ifdef BENCHMARK
		snapshot_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif
		
		// power spectra
		if (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		{
			pkcount++;
			
			COUT << " writing power spectra at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			
			if (sim.out_pk & MASK_DBARE || sim.out_pk & MASK_POT || (sim.out_pk & MASK_T00 && sim.gr_flag == 0))
			{
				if (!redo_psi && sim.gr_flag > 0)
				{
					projection_init(&psi);
					scalarProjectionCIC_project(&pcls_cdm, &psi);
					for (i = 0; i < cosmo.num_ncdm; i++)
						scalarProjectionCIC_project(pcls_ncdm+i, &psi);
					scalarProjectionCIC_comm(&psi);
					redo_psi = 1;
					plan_psi.execute(FFT_FORWARD);
				}
				else if (sim.gr_flag > 0)
				{
					plan_psi.execute(FFT_FORWARD);
				}
				else
				{
					plan_source.execute(FFT_FORWARD);
				}
				
				if (sim.out_pk & MASK_DBARE || (sim.out_pk & MASK_T00 && sim.gr_flag == 0))
					extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				
				if (sim.out_pk & MASK_DBARE)
				{
					sprintf(filename, "%s%s%03d_deltaN.dat", sim.output_path, sim.basename_pk, pkcount-1);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of delta", a);
				}
				
				if (sim.out_pk & MASK_T00 && sim.gr_flag == 0)
				{
					sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount-1);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of T00", a);
				}
				
				if (sim.out_pk & MASK_POT)
				{
					solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);
					extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_psiN.dat", sim.output_path, sim.basename_pk, pkcount-1);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of psi_N", a);
				}
				
				if (cosmo.num_ncdm > 0 && (sim.out_pk & MASK_DBARE || (sim.out_pk & MASK_T00 && sim.gr_flag == 0)))
				{
					redo_psi = 1;
					projection_init(&psi);
					scalarProjectionCIC_project(&pcls_cdm, &psi);
					scalarProjectionCIC_comm(&psi);
					plan_psi.execute(FFT_FORWARD);
					extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_cdm.dat", sim.output_path, sim.basename_pk, pkcount-1);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of delta for cdm", a);
					for (i = 0; i < cosmo.num_ncdm; i++)
					{
						projection_init(&psi);
						scalarProjectionCIC_project(pcls_ncdm+i, &psi);
						scalarProjectionCIC_comm(&psi);
						plan_psi.execute(FFT_FORWARD);
						extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
						sprintf(filename, "%s%s%03d_ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount-1, i);
						sprintf(buffer, "power spectrum of delta for ncdm %d", i);
						writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, buffer, a);
#ifdef CHECK_B
						// store k-space information for cross-spectra using BiFT_check as temporary array
						if (cosmo.num_ncdm > 1 && i < cosmo.num_ncdm-1 && i < 3)
						{
							for (kFT.first(); kFT.test(); kFT.next())
								BiFT_check(kFT, i) = scalarFT(kFT);
						}
#endif						
					}
#ifdef CHECK_B
					for (i = 0; i < cosmo.num_ncdm-1 && i < 3; i++)
					{
						if (i > 0)
						{
							for (kFT.first(); kFT.test(); kFT.next())
								scalarFT(kFT) = BiFT_check(kFT, i);
						}
						for (j = i+1; j < cosmo.num_ncdm && j <= 3; j++)
						{   
							if (j > i+1)
							{
								for (kFT.first(); kFT.test(); kFT.next())
								{
									tempk = BiFT_check(kFT, 0);
									BiFT_check(kFT, 0) = BiFT_check(kFT, j-1);
									BiFT_check(kFT, j-1) = tempk;
								}
							}
							
							extractCrossSpectrum(scalarFT, BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
							sprintf(filename, "%s%s%03d_ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount-1, i, j);
							if (cosmo.num_ncdm < 4)
								sprintf(buffer, "cross power spectrum of delta for ncdm %d x %d", (i+cosmo.num_ncdm-1)%cosmo.num_ncdm, (j+cosmo.num_ncdm-1)%cosmo.num_ncdm);
							else
								sprintf(buffer, "cross power spectrum of delta for ncdm %d x %d", (6-11*i+5*i*i)/2, i ? 4-j : (j+3)%4);
							writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, buffer, a);
						}
					}
#endif
				}
			}
			
			if (sim.out_pk & MASK_PHI)
			{
				plan_phi.execute(FFT_FORWARD);
				extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_phi.dat", sim.output_path, sim.basename_pk, pkcount-1);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a);
			}
			
			if (sim.out_pk & MASK_CHI)
			{
				plan_chi.execute(FFT_FORWARD);
				extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_chi.dat", sim.output_path, sim.basename_pk, pkcount-1);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of chi", a);
			}
			
			if (sim.out_pk & MASK_HIJ)
			{
				if (!doneTT)
					projectFTtensor(SijFT, SijFT);
				extractPowerSpectrum(SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_hij.dat", sim.output_path, sim.basename_pk, pkcount-1);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij", a);
			}
			
			if (sim.out_pk & MASK_T00 && sim.gr_flag > 0)
			{
				plan_source.execute(FFT_FORWARD);
				extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount-1);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of T00", a);
				
				if (cosmo.num_ncdm > 0)
				{
					redo_psi = 1;
					projection_init(&psi);
					projection_T00_project(&pcls_cdm, &psi, a, &phi);
					projection_T00_comm(&psi);
					plan_psi.execute(FFT_FORWARD);
					extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_T00cdm.dat", sim.output_path, sim.basename_pk, pkcount-1);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of T00 for cdm", a);
					for (i = 0; i < cosmo.num_ncdm; i++)
					{
						projection_init(&psi);
						projection_T00_project(pcls_ncdm+i, &psi, a, &phi);
						projection_T00_comm(&psi);
						plan_psi.execute(FFT_FORWARD);
						extractPowerSpectrum(scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
						sprintf(filename, "%s%s%03d_T00ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount-1, i);
						sprintf(buffer, "power spectrum of T00 for ncdm %d", i);
						writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, buffer, a);						
#ifdef CHECK_B
						// store k-space information for cross-spectra using BiFT_check as temporary array
						if (cosmo.num_ncdm > 1 && i < 3)
						{
							for (kFT.first(); kFT.test(); kFT.next())
								BiFT_check(kFT, i) = scalarFT(kFT);
						}
#endif
					}
#ifdef CHECK_B
					for (i = 0; i < cosmo.num_ncdm-1 && i < 3; i++)
					{
						if (i > 0)
						{
							for (kFT.first(); kFT.test(); kFT.next())
								scalarFT(kFT) = BiFT_check(kFT, i);
						}
						for (j = i+1; j < cosmo.num_ncdm && j <=3; j++)
						{   
							if (j > i+1)
							{
								for (kFT.first(); kFT.test(); kFT.next())
								{
									tempk = BiFT_check(kFT, 0);
									BiFT_check(kFT, 0) = BiFT_check(kFT, j-1);
									BiFT_check(kFT, j-1) = tempk;
								}
							}
							
							extractCrossSpectrum(scalarFT, BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
							sprintf(filename, "%s%s%03d_T00ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount-1, i, j);
							if (cosmo.num_ncdm < 4)
								sprintf(buffer, "cross power spectrum of T00 for ncdm %d x %d", (i+cosmo.num_ncdm-1)%cosmo.num_ncdm, (j+cosmo.num_ncdm-1)%cosmo.num_ncdm);
							else
								sprintf(buffer, "cross power spectrum of T00 for ncdm %d x %d", (6-11*i+5*i*i)/2, i ? 4-j : (j+3)%4);
							writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, buffer, a);
						}
					}
#endif
				}
			}
			
			if (sim.out_pk & MASK_B)
			{
				extractPowerSpectrum(BiFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_B.dat", sim.output_path, sim.basename_pk, pkcount-1);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a);
			
#ifdef CHECK_B
				projection_init(&Bi_check);
				vectorProjectionCICNGP_project(&pcls_cdm,&Bi_check);
				for (i = 0; i < cosmo.num_ncdm; i++)
					vectorProjectionCICNGP_project(pcls_ncdm+i, &Bi_check);
				vectorProjectionCICNGP_comm(&Bi_check);
				plan_Bi_check.execute(FFT_FORWARD);
				projectFTvector(BiFT_check, BiFT_check, fourpiG * dx * dx);
				extractPowerSpectrum(BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_B_check.dat", sim.output_path, sim.basename_pk, pkcount-1);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a);
#endif
			}		
		}   // power spectra done
		
#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif	  
		
		if (sim.gr_flag > 0)
		{
			if (redo_psi)
			{
				for (x.first(); x.test(); x.next())
					psi(x) = phi(x) - chi(x);					
				redo_psi = 0;
			}
			
			T00hom = 0.;
			for (x.first(); x.test(); x.next())
				T00hom += source(x);
			parallel.sum<Real>(T00hom);
			T00hom /= (Real) numpts3d;
			
			if (cycle % 10 == 0)
			{
				COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
			}
			
			prepareFTsource<Real>(phi, psi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update

#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_source.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
		
			solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);  // phi update (k-space)
		}
		else
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
		
			solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)
		}

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif		
		plan_phi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif	
		phi.updateHalo();  // communicate halo values

		// record some background data
		if (kFT.setCoord(0, 0, 0))
		{
			sprintf(filename, "%s%s_background.dat", sim.output_path, sim.basename_generic);
			outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for background output!" << endl;
			}
			else
			{
				if (cycle == 0)
					fprintf(outfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
				fprintf(outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau, a, Hconf(a, fourpiG, cosmo) / Hconf(1., fourpiG, cosmo), scalarFT(kFT).real(), T00hom);
				fclose(outfile);
			}
		}
		// done recording background data

		if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot) break; // simulation complete   
		
		prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif		
		plan_Sij.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count += 6;
#endif
		
		projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)
		
		evolveFTvector(SijFT, BiFT, a * a * dtau_old);  // evlolve B using vector projection (k-space)

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif		
		plan_chi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif	
		chi.updateHalo();  // communicale halo values
		
		if (sim.gr_flag > 0)
		{
			for (x.first(); x.test(); x.next())  // update psi
				psi(x) = phi(x) - chi(x);

#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif				
			plan_Bi.execute(FFT_BACKWARD);  // go back to position space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count += 3;
#endif
			Bi.updateHalo();  // communicate halo values
		}
		else
		{
			for (x.first(); x.test(); x.next())  // Newton: update psi
				psi(x) = phi(x);
		}
		
		psi.updateHalo();  // communicate halo values
		
		// compute number of step subdivisions for particle updates
		numsteps = 1;
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (dtau * maxvel_ncdm[i] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau * maxvel_ncdm[i] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;
			
			if (numsteps < numsteps_ncdm[i]) numsteps = numsteps_ncdm[i];
		}
		if (numsteps > 1 && numsteps % 2 > 0) numsteps++;   // if >1, make it an even number
		
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (numsteps / numsteps_ncdm[i] <= 1) numsteps_ncdm[i] = numsteps;
			else if (numsteps_ncdm[i] > 1) numsteps_ncdm[i] = numsteps / 2;
		}
		
		if (cycle % 10 == 0)
		{
			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel << " (cdm Courant factor = " << maxvel * sim.Cf << "), time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau;
			
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (i == 0)
				{
					COUT << endl << " time step subdivision for ncdm species: ";
				}
				COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel_ncdm[i] << ")";
				if (i < cosmo.num_ncdm-1)
				{
					COUT << ", ";
				}
			}
			
			COUT << endl;
		}
		
		for (j = 0; j < numsteps; j++) // particle update
		{
#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			if (j == 0)
			{
				if (sim.gr_flag > 0)
					maxvel = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 1), f_params);
				else
					maxvel = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, 1, f_params);

#ifdef BENCHMARK
				update_q_count++;
#endif
			}
				
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (j % (numsteps / numsteps_ncdm[i]) == 0)
				{
					if (sim.gr_flag > 0)
						maxvel_ncdm[i] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						maxvel_ncdm[i] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, 1, f_params);

#ifdef BENCHMARK
					update_q_count++;
#endif
				}
			}
#ifdef BENCHMARK
			update_q_time += MPI_Wtime() - ref2_time;
			ref2_time = MPI_Wtime();
#endif
			if (numsteps > 1 && j == numsteps / 2)
			{
				if (sim.gr_flag > 0)
					pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
				else
					pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);

#ifdef BENCHMARK
				moveParts_count++;
#endif
			}
		
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (numsteps > 1 && ((numsteps_ncdm[i] == 1 && j == numsteps / 2) || (numsteps_ncdm[i] == numsteps / 2 && j % 2 > 0)))
				{
					if (sim.gr_flag > 0)
						pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);

#ifdef BENCHMARK
					moveParts_count++;
#endif
				}	  
			}
#ifdef BENCHMARK
			moveParts_time += MPI_Wtime() - ref2_time;
#endif  
			rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step
#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			if (numsteps == 1)
			{
				if (sim.gr_flag > 0)
					pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
				else
					pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);

#ifdef BENCHMARK
				moveParts_count++;
#endif
			}
				
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (numsteps_ncdm[i] == numsteps)
				{
					if (sim.gr_flag > 0)
						pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);

#ifdef BENCHMARK
					moveParts_count++;
#endif
				}
			}
#ifdef BENCHMARK
			moveParts_time += MPI_Wtime() - ref2_time;
#endif
			rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step
		}   // particle update done
		
		parallel.max<double>(maxvel);
		if (cosmo.num_ncdm > 0) parallel.max<double>(maxvel_ncdm, cosmo.num_ncdm);
		
		if (sim.gr_flag > 0)
		{
			maxvel /= sqrt(maxvel * maxvel + 1.0);
			for (i = 0; i < cosmo.num_ncdm; i++)
				maxvel_ncdm[i] /= sqrt(maxvel_ncdm[i] * maxvel_ncdm[i] + 1.0);
		}
		
		tau += dtau;
		
		if (tau_Lambda < 0. && (cosmo.Omega_m / a / a / a) < cosmo.Omega_Lambda)
		{
			tau_Lambda = tau;
			COUT << "matter-dark energy equality at z=" << ((1./a) - 1.) << endl;
		}
		
		dtau_old = dtau;
		
		if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
			dtau = sim.Cf * dx;
		else
			dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);
		   
		cycle++;
		
#ifdef BENCHMARK
		gravity_solver_time += MPI_Wtime() - ref_time;
		cycle_time += MPI_Wtime()-cycle_start_time;
#endif
	}
	
	COUT << " simulation complete." << endl;
	
#ifdef BENCHMARK
	run_time = MPI_Wtime() - start_time;

	parallel.sum(run_time);
	parallel.sum(cycle_time);
	parallel.sum(projection_time);
	parallel.sum(snapshot_output_time);
	parallel.sum(spectra_output_time);
	parallel.sum(gravity_solver_time);
	parallel.sum(fft_time);
	parallel.sum(update_q_time);
	parallel.sum(moveParts_time);
	
	COUT << endl << "BENCHMARK" << endl;   
	COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
	COUT << "total number of cycles: "<< cycle << endl;
	COUT << "time consumption breakdown:" << endl;
	COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
	COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;
	
	COUT << "----------- main loop: components -----------"<<endl;
	COUT << "projections          : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
	COUT << "snapshot outputs     : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
	COUT << "power spectra outputs: "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
	COUT << "gravity solver       : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
	
	COUT << "----------- gravity solver: components ------------"<<endl;
	COUT << "Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
	COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/gravity_solver_time <<"%."<<endl;
	COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/gravity_solver_time <<"%."<<endl;
#endif
	
	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);

#ifdef EXTERNAL_IO	
		ioserver.stop();
	}
#endif
}

