//cuda includes
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/tuple.h>
#include <thrust/iterator/counting_iterator.h>
#include <curand_kernel.h>


//standard includes
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>


using namespace std;

//cuMesh includes

#include "cuMesh_define.h"
#include "cuMesh_structs.h"
#include "cuMesh_IO.h"
#include "cuMesh_init.h"
#include "cuMesh_functor_simple_ising_mmc.h"
#include "cuMesh_sweep.h"
#include "cuMesh_measurement.h"

//time measure functions
cudaEvent_t start;
cudaEvent_t end;

void timer_start()
{
	cudaEventCreate(&start); 
	cudaEventCreate(&end);
	cudaEventRecord(start,0);
}

float timer_stop_and_display()
{
	float elapsed_time;
	cudaEventRecord(end, 0);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&elapsed_time, start, end);

	std::cerr << "  ( "<< elapsed_time << "ms )" << std::endl;

	return elapsed_time;
}

int main(int argc, char **argv)
{

	if(argc == 1) {
		cerr << "Usage: spinGlassThrust paramFile or spinGlassThrust paramFile (Sfile) (Jfile) (Hfile)" << endl << endl << "All-to-all correlation in each epoch" << endl;
		exit(-1);
	}

	FILE * pfile;
	pfile = fopen (argv[1], "r");
	if (pfile == NULL){
	cout << "Unable to open file" << endl; 
	exit (EXIT_FAILURE);
	}
	
	int start = time(NULL);
	
	int t_pre  = atoi(load_Parameter("t_pre",pfile).c_str());		//number of preliminary MC steps
	int t_meas = atoi(load_Parameter("t_meas",pfile).c_str());	//number of measurement MC steps		
	int f_meas = atoi(load_Parameter("f_meas",pfile).c_str());	//frequency of measurements (in MC steps)
#ifdef SAVEOUT
	int f_save = atoi(load_Parameter("f_save",pfile).c_str());	//frequency of saveouts (in MC steps)
#endif
	
	int m = atoi(load_Parameter("m_size",pfile).c_str());		//number of rows	
	int n = atoi(load_Parameter("n_size",pfile).c_str());		//number of columns
#ifdef MESH_3D 
	int o = atoi(load_Parameter("o_size",pfile).c_str());	;	//number of layers
#endif	
	float kT = atoi(load_Parameter("kT",pfile).c_str());	;	//temperature
	float PJ = atof(load_Parameter("PJ",pfile).c_str());		//in case of J=+-1 the probability of +1 edges, in case of Gaussian J the sigma of the distribution (only used if no Jfile added)
	
#ifdef FIX_EDGE
	#ifdef PER_BOUND
	int fixedEdge = 1;		//one fixed layer needed for periodic boundary conditions
	cerr<< "Periodic boundary contitions set." << endl;
	#else
	int fixedEdge = atoi(load_Parameter("fixed_edge",pfile).c_str());		//layers to fix
	cerr<< fixedEdge <<" edges will be fixed." << endl;
	#endif
#endif
#ifdef FIX_SPIN
	int fixSpin = atoi(load_Parameter("fix_spin",pfile).c_str());		//spin to fix	
	cerr<< fixSpin <<" spin will be fixed." << endl;	
#endif
#ifdef GLOBAL_FIELD
	float h = atof(load_Parameter("h_field",pfile).c_str());		//magnetic field (only used if no Hfile added)
	cerr<< h <<" will be the uniform magnetic field." << endl;
#endif
#ifdef J_FIX
	float JF = atof(load_Parameter("h_field",pfile).c_str());		//fix value of J 
	cerr<< JF <<" will be the uniform Jij coupling." << endl;
#endif

//generating infoline

	stringstream infoline;
	
	infoline << "#cuMesh spinglass parameter file: ";
	infoline << " m " << m << " n " << n;
#ifdef MESH_3D 
	infoline << " o " << o;
	#ifdef MESH_DIA
	infoline << " Diamond 3D mesh";
	#else
	infoline << " Simple 3D mesh";
	#endif
#else
	#ifdef MESH_TRI
	infoline << " Triangular 2D mesh";
	#else
		#ifdef MESH_HEX
	infoline << " Hexagonal 2D mesh";	
		#else 
	infoline << " Simple 2D mesh";	
		#endif
	#endif
#endif
#ifdef PER_BOUND
	infoline << " with pediodic boundary conditions";	
#else
	infoline << " with non-periodic boundary conditions";
#endif


//initializing J, spin, and h vectors

	int size = m*n
#ifdef MESH_3D 
	*o
#endif		
	;
	
	//initializing cuRand gererator based on CURAND_GEN_TYPE and INIT_SEED define parameters, and heating it up
	curandGenerator_t gen;
	curandCreateGenerator(&gen,CURAND_GEN_TYPE);
	curandSetPseudoRandomGeneratorSeed(gen,INIT_SEED);
	gen = heatup_cuRand(gen,time(NULL));
	
	
	//spin vectors and configuration
	thrust::host_vector<spint> spin_host(size); //for displaying output
	thrust::device_vector<spint> spin1(size);
	thrust::device_vector<spint> spin2(size);
#if ( defined MESH_3D || defined MESH_HEX)
	thrust::device_vector<spint> spin3(size);
#endif

	if (argc > 2){
		load_Svector_int(argv[2], spin_host, size);
	}else{
		init_ising_spin(spin1, gen, size);
	}
#ifdef PER_BOUND
	per_bound_copy(spin_host, m, n
	#ifdef MESH_3D
		,o
	#endif
		)
#endif
	
	//J vectors and configuration
	
	thrust::host_vector<J_str> J_host(size);
	thrust::device_vector<J_str> J(size);
	
#ifndef J_FIX
	
	init_J_idx(J_host, size);
	
	if (argc > 3){
		load_Jvector_int(argv[3], J_host, size);
	}else{
		itnit_ising_J_ij(J_host, gen, size, m, n
	#ifdef MESH_3D
					,o
	#endif
		)
	}
#endif
	init_seed_thread(J_host,size
#ifdef RAND_OWN
	, gen
#endif	
	);

#ifdef LOCAL_FIELD
	//H value configuration
	if (argc > 4){
		load_Hvector(argv[4], J_host, size);
	}
#endif
#if !( defined J_GAUSS || defined GLOBAL_FIELD || defined LOCAL_FIELD)
	//setup subsidiary exponential vectors
	init_expE(J_host, kT, size);
#endif

#ifdef MEAS_MAGN	
	//for magnetism calculation
	float mm = 0.0f;
	float m_avg = 0.0f;
	float m2_avg = 0.0f;
	float m4_avg = 0.0f;	
	//for spin average calculation
	thrust::host_vector<float> sumSpin_host(size);
	thrust::device_vector<float> sumSpin(size);
	thrust::fill(sumSpin.begin(), sumSpin.end(), 0.0f);
	#ifdef MEAS_CORR	
	//for spin correlation calculation
	thrust::host_vector<float> corrSpin_host;
	thrust::device_vector<float> corrSpin(size);
	thrust::fill(corrSpin.begin(), corrSpin.end(), 0.0f);
	#endif
#endif

//Starting thermal sweeps

	J = J_host;
	spin1 = spin_host;

	timer_start();
	cerr<<"[initial sweeps]" << endl;
	
	for(int t = 0; t < t_pre; t++){
	
	//checkerboard update
#if ( defined MESH_3D || defined MESH_HEX)
		sweep_3sub_ising_metropolis(spin1, spin2, spin3, J, kT, m, n
	#ifdef MESH_3D
		,o
	#endif
	#ifdef FIX_EDGE
		,fixedEdge
	#endif
	#ifdef FIX_SPIN
		,fixSpin
	#endif
	#ifdef GLOBAL_FIELD
		,h
	#endif
#ifdef J_FIX
		,JF
	#endif
		);
#else
		sweep_2sub_ising_metropolis(spin1, spin2, J, kT, m, n
	#ifdef FIX_EDGE
		,fixedEdge
	#endif
	#ifdef FIX_SPIN
		,fixSpin
	#endif
	#ifdef GLOBAL_FIELD
		,h
	#endif
	#ifdef J_FIX
		,JF
	#endif
		);
#endif
#ifdef PER_BOUND
		per_bound_copy(spin1, m, n
	#ifdef MESH_3D
		,o
	#endif
		)
#endif
	}
	timer_stop_and_display();

//Starting measure sweeps

	timer_start();
	cerr<<"[averaging sweeps]" << endl;

	for(int t = 0; t < t_meas; t++){	
	
	//checkerboard update
#if ( defined MESH_3D || defined MESH_HEX)
		sweep_3sub_ising_metropolis(spin1, spin2, spin3, J, kT, m, n
	#ifdef MESH_3D
		,o
	#endif
	#ifdef FIX_EDGE
		,fixedEdge
	#endif
	#ifdef FIX_SPIN
		,fixSpin
	#endif
	#ifdef GLOBAL_FIELD
		,h
	#endif
#ifdef J_FIX
		,JF
	#endif
		);
#else
		sweep_2sub_ising_metropolis(spin1, spin2, J, kT, m, n
	#ifdef FIX_EDGE
		,fixedEdge
	#endif
	#ifdef FIX_SPIN
		,fixSpin
	#endif
	#ifdef GLOBAL_FIELD
		,h
	#endif
	#ifdef J_FIX
		,JF
	#endif
		);
#endif
#ifdef PER_BOUND
		per_bound_copy(spin1, m, n
	#ifdef MESH_3D
		,o
	#endif
		)
#endif

#ifdef MEAS_MAGN	
	//Measurement steps
		mm = meas_sample_avg(spin1.begin(), spin1.end(), thrust::dev_ptr<int> sum_begin);
		m_avg += mm;
		m2_avg += mm^2;
		m4_avg += mm^4;
	#ifdef MEAS_CORR		
		if( t%f_meas==0 ){
			meas_time_avg(spin1.begin(), spin1.end(), sumSpin.begin());
			meas_cross_corr(spin1.begin(), spin1.end(), corrSpin.begin(), size, n, m 
		#ifdef MESH_3D
			, o
		#endif
			)
		}
	#endif
#endif

#ifdef SAVEOUT		
	//Frequently saveouts
		if( t%f_meas==0 ){
			spin_host = spin1;
			saveout_spin(infoline, spin_host, start, t, size, kT, n, m 
	#ifdef MESH_3D
			, o
	#endif
			);
	#ifdef MEAS_MAGN
			sumSpin_host = sumSpin;
			saveout_magn(infoline, sumSpin_host, start, t, size, f_meas, kT, n, m 
			#ifdef MESH_3D
			, o
			#endif
			);
			corrSpin_host = corrSpin;
			saveout_corr(infoline, sumSpin_host, corrSpin_host, start, t, f_meas, size, kT, n, m 
			#ifdef MESH_3D
			, o
			#endif
			);
	#endif
		}
#endif
	}
	
	
//Final calculation and saveouts
	timer_start();
	cerr<<"[final calculations and saveout]" << endl;

		spin_host = spin1;
			saveout_spin(infoline, spin_host, start, t_meas, size, kT, n, m 
	#ifdef MESH_3D
			, o
	#endif
		);
	#ifdef MEAS_MAGN
	sumSpin_host = sumSpin;
		saveout_magn(infoline, sumSpin_host, start, t_meas, f_meas, size, kT, n, m 
		#ifdef MESH_3D
			, o
		#endif
		);
		#ifdef MEAS_CORR
	corrSpin_host = corrSpin;
		saveout_corr(infoline, sumSpin_host, corrSpin_host, start, t_meas, f_meas, size, kT, n, m 
			#ifdef MESH_3D
			, o
			#endif
		);
		#endif
		
	cout<<"Binder cummulant: " << 1-((m4_avg/t_meas)/(3*((m2_avg/t_meas)^2))) <<endl;
	
	#endif
	
	cerr<<"[Exiting... Bye!]" << endl;
	
	return 0;
}