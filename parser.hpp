//////////////////////////
// parser.hpp
//////////////////////////
// 
// Parser for settings file
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: March 2016
//
//////////////////////////

#ifndef PARSER_HEADER
#define PARSER_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "metadata.hpp"

#ifndef PARAM_MAX_LENGTH
#define PARAM_MAX_LENGTH 64
#endif
#ifndef PARAM_MAX_LINESIZE
#define PARAM_MAX_LINESIZE 256
#endif

using namespace std;

struct parameter
{
	char name[PARAM_MAX_LENGTH];
	char value[PARAM_MAX_LENGTH];
	bool used;
};


int compare_redshifts (const void * z1, const void * z2)
{
	if (* (double *) z1 > * (double *) z2)
		return -1;
    else if (* (double *) z1 < * (double *) z2)
		return 1;
	else
		return 0;
}


//////////////////////////
// readline
//////////////////////////
// Description:
//   reads a line of characters and checks if it declares a parameter; if yes,
//   i.e. the line has the format "<parameter name> = <parameter value>" (with
//   an optional comment added, preceded by a hash-symbol, '#'), the parameter
//   name and value are copied to the corresponding arrays, and 'true' is returned.
//   If the format is not recognized (or the line is commented using the hash-symbol)
//   'false' is returned instead.
// 
// Arguments:
//   line       string containing the line to be read
//   pname      will contain the name of the declared parameter (if found)
//   pvalue     will contain the value of the declared parameter (if found)
//
// Returns:
//   'true' if a parameter is declared in the line, 'false' otherwise.
// 
//////////////////////////

bool readline(char * line, char * pname, char * pvalue)
{
	char * pequal;
	char * phash;
	char * l;
	char * r;
	
	pequal = strchr(line, '=');
	
	if (pequal == NULL || pequal == line) return false;
	
	phash = strchr(line, '#');
	
	if (phash != NULL && phash < pequal) return false;
	
	l = line;
	while (*l == ' ' || *l == '\t') l++;
	
	r = pequal-1;
	while ((*r == ' ' || *r == '\t') && r > line) r--;
	
	if (r < l) return false;
	
	if (r-l+1 >= PARAM_MAX_LENGTH) return false;
	
	strncpy(pname, l, r-l+1);
	pname[r-l+1] = '\0';
	
	l = pequal+1;
	while (*l == ' ' || *l == '\t') l++;
	
	if (phash == NULL)
		r = line+strlen(line)-1;
	else
		r = phash-1;
	  
	while (*r == ' ' || *r == '\t' || *r == '\n' || *r == '\r') r--;
	
	if (r < l) return false;
	
	if (r-l+1 >= PARAM_MAX_LENGTH) return false;
	
	strncpy(pvalue, l, r-l+1);
	pvalue[r-l+1] = '\0';
	
	return true;
}


//////////////////////////
// loadParameterFile
//////////////////////////
// Description:
//   loads a parameter file and creates an array of parameters declared therein
// 
// Arguments:
//   filename   string containing the path to the parameter file
//   params     will contain the array of parameters (memory will be allocated)
//
// Returns:
//   number of parameters defined in the parameter file (= length of parameter array)
// 
//////////////////////////

int loadParameterFile(const char * filename, parameter * & params)
{
	int numparam = 0;
	int i = 0;
	
	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * paramfile;
		char line[PARAM_MAX_LINESIZE];
		char pname[PARAM_MAX_LENGTH];
		char pvalue[PARAM_MAX_LENGTH];
		
		paramfile = fopen(filename, "r");
		
		if (paramfile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Unable to open parameter file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		while (!feof(paramfile) && !ferror(paramfile))
		{
			fgets(line, PARAM_MAX_LINESIZE, paramfile);
			
			if (readline(line, pname, pvalue) == true) numparam++;
		}
		
		if (numparam == 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! No valid data found in file " << filename << "." << endl;
			fclose(paramfile);
			parallel.abortForce();
		}
		
		params = (parameter *) malloc(sizeof(parameter) * numparam);
		
		if (params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			fclose(paramfile);
			parallel.abortForce();
		}
		
		rewind(paramfile);
		
		while (!feof(paramfile) && !ferror(paramfile) && i < numparam)
		{
			fgets(line, PARAM_MAX_LINESIZE, paramfile);
			
			if (readline(line, params[i].name, params[i].value) == true)
			{
				params[i].used = false;
				i++;
			}
		}
		
		fclose(paramfile);
		
		if (i < numparam)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! File may have changed or file pointer corrupted." << endl;
			free(params);
			parallel.abortForce();
		}
		
		parallel.broadcast_dim0<int>(numparam, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numparam, 0);
		
		params = (parameter *) malloc(sizeof(parameter) * numparam);
		
		if (params == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadParameterFile! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<parameter>(params, numparam, 0);
	
	return numparam;
}


//////////////////////////
// saveParameterFile
//////////////////////////
// Description:
//   saves a parameter file
// 
// Arguments:
//   filename   string containing the path to the parameter file
//   params     array of parameters
//   numparam   length of parameter array
//   used_only  if 'true', only the used parameters will be written (default)
//
// Returns:
// 
//////////////////////////

void saveParameterFile(const char * filename, parameter * params, const int numparam, bool used_only = true)
{
	if (parallel.isRoot())
	{
		FILE * paramfile;
		
		paramfile = fopen(filename, "w");
		
		if (paramfile == NULL)
		{
			cout << " error in saveParameterFile! Unable to open file " << filename << "." << endl;
		}
		else
		{
			for (int i = 0; i < numparam; i++)
			{
				if (!used_only || params[i].used)
					fprintf(paramfile, "%s = %s\n", params[i].name, params[i].value);
			}
			
			fclose(paramfile);
		}
	}
}


//////////////////////////
// parseParameter (int)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, int & pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%d", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (long)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as integer
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, long & pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%ld", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (double)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value as double
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to double which will contain the parsed parameter value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, double & pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			if (sscanf(params[i].value, "%lf", &pvalue) == 1)
			{
				params[i].used = true;
				return true;
			}
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (char *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and retrieves its value as string
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     character string which will contain a copy of the parameter value (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, char * pvalue)
{   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			strcpy(pvalue, params[i].value);
			params[i].used = true;
			return true;
		}
	}
	
	return false;
}


//////////////////////////
// parseParameter (double *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated double values
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of double values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, double * pvalue, int & nmax)
{
	char * start;
	char * comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;
	   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma-start);
					item[comma-start] = '\0';
					if (sscanf(item, " %lf ", pvalue+n) != 1)
					{
						nmax = n;
						return false;
					}
					if (++n > nmax-2)
						break;
					start = comma+1;
				}
			}   
			if (sscanf(start, " %lf ", pvalue+n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}
	
	nmax = 0;
	return false;
}


//////////////////////////
// parseParameter (int *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated integer values
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of integer values which will contain the list of parsed parameters (if found)
//   nmax       maximum size of array; will be set to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseParameter(parameter * & params, const int numparam, const char * pname, int * pvalue, int & nmax)
{
	char * start;
	char * comma;
	char item[PARAM_MAX_LENGTH];
	int n = 0;
	   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			start = params[i].value;
			if (nmax > 1)
			{
				while ((comma = strchr(start, ',')) != NULL)
				{
					strncpy(item, start, comma-start);
					item[comma-start] = '\0';
					if (sscanf(item, " %d ", pvalue+n) != 1)
					{
						nmax = n;
						return false;
					}
					if (++n > nmax-2)
						break;
					start = comma+1;
				}
			}   
			if (sscanf(start, " %d ", pvalue+n) != 1)
			{
				nmax = n;
				return false;
			}
			nmax = ++n;
			params[i].used = true;
			return true;
		}
	}
	
	nmax = 0;
	return false;
}


//////////////////////////
// parseFieldSpecifiers
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a list of comma-separated field specifiers
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     integer which will contain the binary-encoded list of parsed specifiers (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
// 
//////////////////////////

bool parseFieldSpecifiers(parameter * & params, const int numparam, const char * pname, int & pvalue)
{
	char * start;
	char * comma;
	int pos;
	char item[PARAM_MAX_LENGTH];
		   
	for (int i = 0; i < numparam; i++)
	{
		if (strcmp(params[i].name, pname) == 0)
		{
			pvalue = 0;
			start = params[i].value;
			while ((comma = strchr(start, ',')) != NULL)
			{
				strncpy(item, start, comma-start);
				for (pos = comma-start; pos > 0; pos--)
				{
					if (item[pos-1] != ' ' && item[pos-1] != '\t') break;
				}
				item[pos] = '\0';
				
				if (strcmp(item, "Phi") == 0 || strcmp(item, "phi") == 0)
					pvalue |= MASK_PHI;
				else if (strcmp(item, "Chi") == 0 || strcmp(item, "chi") == 0)
					pvalue |= MASK_CHI;
				else if (strcmp(item, "Pot") == 0 || strcmp(item, "pot") == 0 || strcmp(item, "Psi_N") == 0 || strcmp(item, "psi_N") == 0 || strcmp(item, "PsiN") == 0 || strcmp(item, "psiN") == 0)
					pvalue |= MASK_POT;
				else if (strcmp(item, "B") == 0 || strcmp(item, "Bi") == 0)
					pvalue |= MASK_B;
				else if (strcmp(item, "P") == 0 || strcmp(item, "p") == 0 || strcmp(item, "v") == 0)
					pvalue |= MASK_P;
				else if (strcmp(item, "T00") == 0 || strcmp(item, "rho") == 0)
					pvalue |= MASK_T00;
				else if (strcmp(item, "Tij") == 0)
					pvalue |= MASK_TIJ;
				else if (strcmp(item, "delta_N") == 0 || strcmp(item, "deltaN") == 0)
					pvalue |= MASK_DBARE;
				else if (strcmp(item, "hij") == 0 || strcmp(item, "GW") == 0)
					pvalue |= MASK_HIJ;
				else if (strcmp(item, "Gadget") == 0 || strcmp(item, "Gadget2") == 0 || strcmp(item, "gadget") == 0 || strcmp(item, "gadget2") == 0)
					pvalue |= MASK_GADGET;
				else if (strcmp(item, "Particles") == 0 || strcmp(item, "particles") == 0 || strcmp(item, "pcls") == 0 || strcmp(item, "part") == 0)
					pvalue |= MASK_PCLS;
					
				start = comma+1;
				while (*start == ' ' || *start == '\t') start++;
			}  
			
			if (strcmp(start, "Phi") == 0 || strcmp(start, "phi") == 0)
				pvalue |= MASK_PHI;
			else if (strcmp(start, "Chi") == 0 || strcmp(start, "chi") == 0)
				pvalue |= MASK_CHI;
			else if (strcmp(start, "Pot") == 0 || strcmp(start, "pot") == 0 || strcmp(start, "Psi_N") == 0 || strcmp(start, "psi_N") == 0 || strcmp(start, "PsiN") == 0 || strcmp(start, "psiN") == 0)
				pvalue |= MASK_POT;
			else if (strcmp(start, "B") == 0 || strcmp(start, "Bi") == 0)
				pvalue |= MASK_B;
			else if (strcmp(start, "P") == 0 || strcmp(start, "p") == 0 || strcmp(start, "v") == 0)
				pvalue |= MASK_P;
			else if (strcmp(start, "T00") == 0 || strcmp(start, "rho") == 0)
				pvalue |= MASK_T00;
			else if (strcmp(start, "Tij") == 0)
				pvalue |= MASK_TIJ;
			else if (strcmp(start, "delta_N") == 0 || strcmp(start, "deltaN") == 0)
				pvalue |= MASK_DBARE;
			else if (strcmp(start, "hij") == 0 || strcmp(start, "GW") == 0)
				pvalue |= MASK_HIJ;
			else if (strcmp(start, "Gadget") == 0 || strcmp(start, "Gadget2") == 0 || strcmp(start, "gadget") == 0 || strcmp(start, "gadget2") == 0)
				pvalue |= MASK_GADGET;
			else if (strcmp(start, "Particles") == 0 || strcmp(start, "particles") == 0 || strcmp(start, "pcls") == 0 || strcmp(start, "part") == 0)
				pvalue |= MASK_PCLS; 
			
			params[i].used = true;
			return true;
		}
	}
	
	return false;
}


//////////////////////////
// parseMetadata
//////////////////////////
// Description:
//   parses all metadata from the parameter array
// 
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   sim        reference to metadata stucture (holds simulation parameters)
//   cosmo      reference to cosmology structure (holds cosmological parameters)
//   ic         reference to icsettings structure (holds settings for IC generation)
//
// Returns:
//   number of parameters parsed
// 
//////////////////////////

int parseMetadata(parameter * & params, const int numparam, metadata & sim, cosmology & cosmo, icsettings & ic)
{
	char par_string[PARAM_MAX_LENGTH];
	char * ptr;
	int usedparams = 0;
	int i;

	// parse settings for IC generator
	
	ic.seed = 0;
	ic.flags = 0;
	ic.z_relax = -2.;
	ic.steplimit = 0.2;
	ic.A_s = 2.215e-9;
	ic.n_s = 0.9619;
	
	parseParameter(params, numparam, "seed", ic.seed);
	
	if (parseParameter(params, numparam, "IC generator", par_string))
	{
		if (par_string[0] == 'B' || par_string[0] == 'b')
			ic.generator = ICGEN_BASIC;
		else if (par_string[0] == 'R' || par_string[0] == 'r')
			ic.generator = ICGEN_READ_FROM_DISK;
#ifdef ICGEN_PREVOLUTION
		else if (par_string[0] == 'P' || par_string[0] == 'p')
			ic.generator = ICGEN_PREVOLUTION;
#endif
#ifdef ICGEN_SONG
		else if (par_string[0] == 'S' || par_string[0] == 's')
			ic.generator = ICGEN_SONG;
#endif
#ifdef ICGEN_FALCONIC
		else if (par_string[0] == 'F' || par_string[0] == 'f')
			ic.generator = ICGEN_FALCONIC;
#endif
		else
		{
			COUT << " error: IC generator not recognized!" << endl;
			parallel.abortForce();
		}
	}
	else
	{
		COUT << " warning: IC generator not specified, selecting default (basic)" << endl;
		ic.generator = ICGEN_BASIC;
	}
	
	if (!parseParameter(params, numparam, "template file", ic.pclfile))
	{
		COUT << " error: no template file specified!" << endl;
		parallel.abortForce();
	}
	
	if (!parseParameter(params, numparam, "mPk file", ic.pkfile))
	{
		COUT << " error: no power spectrum file specified!" << endl;
		parallel.abortForce();
	}
	
	if (parseParameter(params, numparam, "correct displacement", par_string))
	{
		if (par_string[0] == 'Y' || par_string[0] == 'y')
			ic.flags |= ICFLAG_CORRECT_DISPLACEMENT;
		else if (par_string[0] != 'N' && par_string[0] != 'n')
			COUT << " warning: setting chosen for deconvolve displacement option not recognized, using default (no)" << endl;
	}
	
	if (parseParameter(params, numparam, "k-domain", par_string))
	{
		if (par_string[0] == 'S' || par_string[0] == 's')
			ic.flags |= ICFLAG_KSPHERE;
		else if (par_string[0] != 'C' && par_string[0] != 'c')
			COUT << " warning: setting chosen for k-domain option not recognized, using default (cube)" << endl;
	}
	
	parseParameter(params, numparam, "tiling factor", ic.numtile);
	if (ic.numtile <= 0)
	{
		COUT << " warning: tiling number for particle template not set properly; using default value (1)" << endl;
		ic.numtile = 1;
	}
	
	parseParameter(params, numparam, "relaxation redshift", ic.z_relax);
	
	cosmo.num_ncdm = MAX_PCL_SPECIES;
	parseParameter(params, numparam, "PSD samples", cosmo.m_ncdm, cosmo.num_ncdm);
	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		ic.psd_samples[i] = (int) floor(cosmo.m_ncdm[i]);
		if (ic.psd_samples[i] < 1) ic.psd_samples[i] = 1;
	}
	for(; i < MAX_PCL_SPECIES; i++) ic.psd_samples[i] = 1;

#ifdef ICGEN_PREVOLUTION	
	if (ic.generator == ICGEN_PREVOLUTION)
	{
		if (parseParameter(params, numparam, "k-mode files", ic.kmodefiles))
		{
			if ((ptr = strchr(ic.kmodefiles, '*')) == NULL)
			{
				COUT << " error: missing wildcard character (*) in k-mode files = " << ic.kmodefiles << endl;
				parallel.abortForce();
			}
			
			if (strchr(++ptr, '*') != NULL)
			{
				COUT << " error: two or more wildcard characters (*) detected in k-mode files = " << ic.kmodefiles << endl;
				parallel.abortForce();
			}
		}
		else
		{
			COUT << " error: no k-mode files specified for IC generator = prevolution" << endl;
			parallel.abortForce();
		}
		
		if (!parseParameter(params, numparam, "background file", ic.bgfile))
		{
			COUT << " error: no background file specified for IC generator = prevolution" << endl;
			parallel.abortForce();
		}
		
		if (!parseParameter(params, numparam, "IC redshift", ic.z_ic))
		{
			COUT << " error: no IC redshift specified for IC generator = prevolution" << endl;
			parallel.abortForce();
		}
		
		if (!parseParameter(params, numparam, "IC Courant factor", ic.Cf))
		{
			COUT << " warning: no IC Courant factor specified for IC generator = prevolution; using default value (1)" << endl;
			ic.Cf = 1.;
		}
	}
#endif

#ifdef ICGEN_FALCONIC
	if (ic.generator == ICGEN_FALCONIC)
	{
		parseParameter(params, numparam, "A_s", ic.A_s);
		
		parseParameter(params, numparam, "n_s", ic.n_s);
	}
#endif
	
	// parse metadata
	
	sim.numpts = 0;
	for (i = 0; i <= MAX_PCL_SPECIES; i++) sim.numpcl[i] = 0;
	sim.gr_flag = 0;
	sim.out_pk = 0;
	sim.out_snapshot = 0;
	sim.num_pk = MAX_OUTPUTS;
	sim.numbins = 0;
	sim.num_snapshot = MAX_OUTPUTS;
	for (i = 0; i <= MAX_PCL_SPECIES; i++) sim.tracer_factor[i] = 1;
	sim.Cf = 1.;
	sim.steplimit = 1.;
	sim.boxsize = -1.;
	sim.z_in = 0.;
	
	if (!parseParameter(params, numparam, "generic file base", sim.basename_generic))
		sim.basename_generic[0] = '\0';
	
	if (!parseParameter(params, numparam, "snapshot file base", sim.basename_snapshot))
		strcpy(sim.basename_snapshot, "snapshot");
		
	if (!parseParameter(params, numparam, "Pk file base", sim.basename_pk))
		strcpy(sim.basename_pk, "pk");
		
	if (!parseParameter(params, numparam, "output path", sim.output_path))
		sim.output_path[0] = '\0';
		
	parseParameter(params, numparam, "boxsize", sim.boxsize);
	if (sim.boxsize <= 0. || !isfinite(sim.boxsize))
	{
		COUT << " error: simulation box size not set properly!" << endl;
		parallel.abortForce();
	}
	
	parseParameter(params, numparam, "Ngrid", sim.numpts);
	if (sim.numpts < 2 || !isfinite(sim.numpts))
	{
		COUT << " error: number of grid points not set properly!" << endl;
		parallel.abortForce();
	}
	
	if (!parseParameter(params, numparam, "move limit", sim.movelimit))
		sim.movelimit = sim.numpts;
	
	if (!parseParameter(params, numparam, "initial redshift", sim.z_in))
	{
		COUT << " error: initial redshift not specified!" << endl;
		parallel.abortForce();
	}
	
	if (ic.z_relax < -1.) ic.z_relax = sim.z_in;
	
	parseParameter(params, numparam, "snapshot redshifts", sim.z_snapshot, sim.num_snapshot);
	if (sim.num_snapshot > 0)
		qsort((void *) sim.z_snapshot, (size_t) sim.num_snapshot, sizeof(double), compare_redshifts);
	
	parseParameter(params, numparam, "Pk redshifts", sim.z_pk, sim.num_pk);
	if (sim.num_pk > 0)
		qsort((void *) sim.z_pk, (size_t) sim.num_pk, sizeof(double), compare_redshifts);
		
	parseFieldSpecifiers(params, numparam, "snapshot outputs", sim.out_snapshot);
	parseFieldSpecifiers(params, numparam, "Pk outputs", sim.out_pk);
	
	i = MAX_PCL_SPECIES+1;
	parseParameter(params, numparam, "tracer factor", sim.tracer_factor, i);
	for (; i > 0; i--)
	{
		if (sim.tracer_factor[i-1] < 1)
		{
			COUT << " warning: tracer factor not set properly; using default value (1)" << endl;
			sim.tracer_factor[i-1] = 1;
		}
	}
	
	if ((sim.num_snapshot <= 0 || sim.out_snapshot == 0) && (sim.num_pk <= 0 || sim.out_pk == 0))
	{
		COUT << " warning: no output specified!" << endl;
	}
	
	if (!parseParameter(params, numparam, "Pk bins", sim.numbins))
	{
		COUT << " warning: number of Pk bins not set properly; using default value (64)" << endl;
		sim.numbins = 64;
	}
	
	parseParameter(params, numparam, "Courant factor", sim.Cf);
	
	parseParameter(params, numparam, "time step limit", sim.steplimit);
	
	if (parseParameter(params, numparam, "gravity theory", par_string))
	{
		if (par_string[0] == 'N' || par_string[0] == 'n')
		{
			COUT << " gravity theory set to: Newtonian" << endl;
			sim.gr_flag = 0;
		}
		else if (par_string[0] == 'G' || par_string[0] == 'g')
		{
			COUT << " gravity theory set to: General Relativity" << endl;
			sim.gr_flag = 1;
		}
		else
		{
			COUT << " warning: gravity theory unknown, using default (General Relativity)" << endl;
			sim.gr_flag = 1;
		}
	}
	else
	{
		COUT << " warning: gravity theory not selected, using default (General Relativity)" << endl;
		sim.gr_flag = 1;
	}
	
	
	// parse cosmological parameters
	
	if (!parseParameter(params, numparam, "h", cosmo.h))
	{
		cosmo.h = 0.67556;
	}
	
	cosmo.num_ncdm = MAX_PCL_SPECIES;
	if (!parseParameter(params, numparam, "m_ncdm", cosmo.m_ncdm, cosmo.num_ncdm))
	{
		for (i = 0; i < MAX_PCL_SPECIES; i++) cosmo.m_ncdm[i] = 0.;
		cosmo.num_ncdm = 0;
	}
	
	if (parseParameter(params, numparam, "N_ncdm", i))
	{
		if (i < 0 || !isfinite(i))
		{
			COUT << " error: number of ncdm species not set properly!" << endl;
			parallel.abortForce();
		}
		if (i > cosmo.num_ncdm)
		{
			COUT << " error: N_ncdm = " << i << " is larger than the number of mass parameters specified (" << cosmo.num_ncdm << ")!" << endl;
			parallel.abortForce();
		}
		cosmo.num_ncdm = i;
	}
	else if (cosmo.num_ncdm > 0)
	{
		COUT << " warning: N_ncdm not specified, inferring from number of mass parameters in m_ncdm (" << cosmo.num_ncdm << ")!" << endl;
	}
	
	for (i = 0; i < MAX_PCL_SPECIES; i++)
	{
		cosmo.T_ncdm[i] = 0.71611;
		cosmo.deg_ncdm[i] = 1.0;	
	}
	parseParameter(params, numparam, "T_ncdm", cosmo.T_ncdm, i);
	i = MAX_PCL_SPECIES;
	parseParameter(params, numparam, "deg_ncdm", cosmo.deg_ncdm, i);
	
	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		cosmo.Omega_ncdm[i] = cosmo.m_ncdm[i] * cosmo.deg_ncdm[i] / 93.14 / cosmo.h / cosmo.h;
	}
	
	if (parseParameter(params, numparam, "T_cmb", cosmo.Omega_g))
	{
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g / cosmo.h;
		cosmo.Omega_g = cosmo.Omega_g * cosmo.Omega_g * 4.48147e-7; // Planck's law
	}
	else if (parseParameter(params, numparam, "omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_g", cosmo.Omega_g))
	{
		cosmo.Omega_g = 0.;
	}
	
	if (parseParameter(params, numparam, "N_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	else if (parseParameter(params, numparam, "N_eff", cosmo.Omega_ur))
	{
		cosmo.Omega_ur *= (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	else if (parseParameter(params, numparam, "omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_ur", cosmo.Omega_ur))
	{
		cosmo.Omega_ur = 3.046 * (7./8.) * pow(4./11., 4./3.) * cosmo.Omega_g;
	}
	
	cosmo.Omega_rad = cosmo.Omega_g + cosmo.Omega_ur;
	
	if (parseParameter(params, numparam, "omega_b", cosmo.Omega_b))
	{
		cosmo.Omega_b /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_b", cosmo.Omega_b))
	{
		COUT << " warning: Omega_b not found in settings file, setting to default (0)." << endl;
		cosmo.Omega_b = 0.;
	}
	
	if (parseParameter(params, numparam, "omega_cdm", cosmo.Omega_cdm))
	{
		cosmo.Omega_cdm /= cosmo.h * cosmo.h;
	}
	else if (!parseParameter(params, numparam, "Omega_cdm", cosmo.Omega_cdm))
	{
		COUT << " warning: Omega_cdm not found in settings file, setting to default (1)." << endl;
		cosmo.Omega_cdm = 1.;
	}
	
	cosmo.Omega_m = cosmo.Omega_cdm + cosmo.Omega_b;
	for (i = 0; i < cosmo.num_ncdm; i++) cosmo.Omega_m += cosmo.Omega_ncdm[i];
	
	if (cosmo.Omega_m <= 0. || cosmo.Omega_m > 1.)
	{
		COUT << " error: total matter density out of range!" << endl;
		parallel.abortForce();
	}
	else if (cosmo.Omega_rad < 0. || cosmo.Omega_rad > 1. - cosmo.Omega_m)
	{
		COUT << " error: total radiation energy density out of range!" << endl;
		parallel.abortForce();
	}
	else
	{
		COUT << " cosmological parameters are: Omega_m0 = " << cosmo.Omega_m << ", Omega_rad0 = " << cosmo.Omega_rad << ", h = " << cosmo.h << endl;
		cosmo.Omega_Lambda = 1. - cosmo.Omega_m - cosmo.Omega_rad;
	}
	
	for (i = 0; i < numparam; i++)
	{
		if (params[i].used) usedparams++;
	}
	
	return usedparams;
}

#endif
