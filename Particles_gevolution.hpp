//////////////////////////
// Particles_gevolution.hpp
//////////////////////////
//
// Author: Julian Adamek (Université de Genève)
//
// Last modified: April 2016
//
//////////////////////////

#ifndef PARTICLES_GEVOLUTION_HEADER
#define PARTICLES_GEVOLUTION_HEADER

using namespace LATfield2;

#ifndef GADGET2_HEADER
#define GADGET2_HEADER
struct gadget2_header
{
	unsigned int npart[6];
	double mass[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	unsigned int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	   /* fills to 256 Bytes */
};
#endif

template <typename part, typename part_info, typename part_dataType>
class Particles_gevolution: public Particles<part, part_info, part_dataType>
{
	public:
		void saveGadget2(string filename, gadget2_header & hdr, const int tracer_factor = 1);
};

template <typename part, typename part_info, typename part_dataType>
void Particles_gevolution<part,part_info,part_dataType>::saveGadget2(string filename, gadget2_header & hdr, const int tracer_factor)
{
#ifndef PCLBUFFER
#define PCLBUFFER 1048576
#endif

	float * posdata;
	float * veldata;
	long * IDs;
	MPI_File outfile;
	long count, npart;
	MPI_Offset offset_pos, offset_vel, offset_ID;
	MPI_Status status;
	unsigned int blocksize;
	unsigned int i;
	char fname[filename.length()+1];
	double rescale_vel = 299792.458 / sqrt(hdr.time);
	
	filename.copy(fname, filename.length());
	fname[filename.length()] = '\0';
	
	LATfield2::Site xPart(this->lat_part_);
	typename std::list<part>::iterator it;
	
	if (hdr.num_files != 1)
	{
		COUT << " warning: writing multiple Gadget2 files not currently supported!" << endl;
		return;
	}
	
	posdata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	IDs = (long *) malloc(3 * sizeof(long) * PCLBUFFER);
	
	npart = 0;
	for(xPart.first(); xPart.test(); xPart.next())
	{
		if(this->field_part_(xPart).size!=0)
		{
			for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
			{
				if ((*it).ID % tracer_factor == 0)
					npart++;
			}
		}
	}
	
	if (parallel.rank() == 0)
	{
		parallel.send<long>(npart, 1);
		parallel.receive<long>(count, parallel.size()-1);
		if (count != hdr.npart[1]) cout << " error: number of particles in saveGadget2 does not match request!" << endl;
		count = 0;
	}
	else
	{
		parallel.receive<long>(count, parallel.rank()-1);
		npart += count;
		parallel.send<long>(npart, (parallel.rank()+1)%parallel.size());
	}
	
	MPI_File_open(parallel.lat_world_comm(), fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &outfile);
	
	offset_pos = (MPI_Offset) hdr.npart[1];
	offset_pos *= (MPI_Offset) (6 * sizeof(float) + sizeof(long));
	offset_pos += (MPI_Offset) (8 * sizeof(unsigned int) + sizeof(hdr));
	MPI_File_set_size(outfile, offset_pos);
	
	offset_pos = (MPI_Offset) (3 * sizeof(unsigned int) + sizeof(hdr)) + ((MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float)));
	offset_vel = offset_pos + (MPI_Offset) (2 * sizeof(unsigned int)) + ((MPI_Offset) hdr.npart[1]) * ((MPI_Offset) (3 * sizeof(float)));
	offset_ID = offset_vel + (MPI_Offset) (2 * sizeof(unsigned int)) + ((MPI_Offset) hdr.npart[1] - (MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float))) + ((MPI_Offset) count) * ((MPI_Offset) sizeof(long));
	
	if (parallel.rank() == 0)
	{
		blocksize = sizeof(hdr);		
		MPI_File_write_at(outfile, 0, &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, sizeof(unsigned int), &hdr, sizeof(hdr), MPI_BYTE, &status);
		MPI_File_write_at(outfile, sizeof(hdr) + sizeof(unsigned int), &blocksize, 1, MPI_UNSIGNED, &status);
		blocksize = 3 * sizeof(float) * hdr.npart[1];
		MPI_File_write_at(outfile, sizeof(hdr) + 2*sizeof(unsigned int), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_vel - 2*sizeof(unsigned int), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_vel - sizeof(unsigned int), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_ID - 2*sizeof(unsigned int), &blocksize, 1, MPI_UNSIGNED, &status);
		blocksize = sizeof(long) * hdr.npart[1];
		MPI_File_write_at(outfile, offset_ID - sizeof(unsigned int), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_ID + blocksize, &blocksize, 1, MPI_UNSIGNED, &status);
	}
	
	count = 0;
	for(xPart.first(); xPart.test(); xPart.next())
	{
		if(this->field_part_(xPart).size!=0)
		{
			for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
			{
				if ((*it).ID % tracer_factor == 0)
				{
					for (i = 0; i < 3; i++)
						posdata[3*count+i] = (*it).pos[i] * hdr.BoxSize;
					
					for (i = 0; i < 3; i++)
						veldata[3*count+i] = (*it).vel[i] * rescale_vel / sqrt(hdr.time * hdr.time + (*it).vel[0] * (*it).vel[0] + (*it).vel[1] * (*it).vel[1] + (*it).vel[2] * (*it).vel[2]);
						
					IDs[count] = (*it).ID;
					
					count++;
						
					if (count == PCLBUFFER)
					{
						MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
						offset_pos += 3 * PCLBUFFER * sizeof(float);
						MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
						offset_vel += 3 * PCLBUFFER * sizeof(float);
						MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_LONG, &status);
						offset_ID += PCLBUFFER * sizeof(long);
						count = 0;
					}
				}
			}
		}
	}
		
	if (count > 0)
	{
			MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
			MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
			MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_LONG, &status);
	}
	
	MPI_File_close(&outfile);
	
	free(posdata);
	free(veldata);
	free(IDs);
}

#endif
