//////////////////////////
// metadata.hpp
//////////////////////////
// 
// Metadata structure
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: March 2016
//
//////////////////////////

#ifndef METADATA_HEADER
#define METADATA_HEADER

#ifndef MAX_OUTPUTS
#define MAX_OUTPUTS 32
#endif

#ifndef PARAM_MAX_LENGTH
#define PARAM_MAX_LENGTH 64
#endif

#ifndef MAX_PCL_SPECIES
#define MAX_PCL_SPECIES 4
#endif

#define MASK_PHI    1
#define MASK_CHI    2
#define MASK_POT    4
#define MASK_B      8
#define MASK_T00    16
#define MASK_TIJ    32
#define MASK_DBARE  64
#define MASK_HIJ    128
#define MASK_P      256
#define MASK_GADGET 512
#define MASK_PCLS   1024

#define ICFLAG_CORRECT_DISPLACEMENT 1
#define ICFLAG_KSPHERE              2

#define ICGEN_BASIC                 0
#define ICGEN_READ_FROM_DISK        1
#ifdef ICGEN_PREVOLUTION
#undef ICGEN_PREVOLUTION
#define ICGEN_PREVOLUTION           2
#endif
#ifdef ICGEN_SONG
#undef ICGEN_SONG
#define ICGEN_SONG                  3
#endif
#ifdef ICGEN_FALCONIC
#undef ICGEN_FALCONIC
#define ICGEN_FALCONIC              4
#endif

struct metadata
{
	int numpts;
	long numpcl[MAX_PCL_SPECIES+1];
	int tracer_factor[MAX_PCL_SPECIES+1];
	int gr_flag;
	int out_pk;
	int out_snapshot;
	int num_pk;
	int numbins;
	int num_snapshot;
	int movelimit;
	double Cf;
	double steplimit;
	double boxsize;
	double z_in;
	double z_snapshot[MAX_OUTPUTS];
	double z_pk[MAX_OUTPUTS];
	char basename_snapshot[PARAM_MAX_LENGTH];
	char basename_pk[PARAM_MAX_LENGTH];
	char basename_generic[PARAM_MAX_LENGTH];
	char output_path[PARAM_MAX_LENGTH];
};

struct icsettings
{
	int numtile;
	int seed;
	int flags;
	int generator;
	int psd_samples[MAX_PCL_SPECIES];
	char pclfile[PARAM_MAX_LENGTH];
	char pkfile[PARAM_MAX_LENGTH];
	char bgfile[PARAM_MAX_LENGTH];
	char kmodefiles[PARAM_MAX_LENGTH];
	double z_ic;
	double z_relax;
	double Cf;
	double steplimit;
	double A_s;
	double n_s;
};

struct cosmology
{
	double Omega_cdm;
	double Omega_b;
	double Omega_m;
	double Omega_Lambda;
	double Omega_g;
	double Omega_ur;
	double Omega_rad;
	double Omega_ncdm[MAX_PCL_SPECIES];
	double h;
	double m_ncdm[MAX_PCL_SPECIES];
	double T_ncdm[MAX_PCL_SPECIES];
	double deg_ncdm[MAX_PCL_SPECIES];
	int num_ncdm;
};

#endif
