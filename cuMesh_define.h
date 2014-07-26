//topology
					//without any define -> simple 2D mesh with coordination number 4
//#define MESH_TRI	//2D tridiagonal mesh with coordination number 3
//#define MESH_HEX	//2D hexagonal mesh with coordination number 6
//#define MESH_3D		//3D regular mesh with coordination number 6
//#define MESH_DIA	//3D diamond mesh with coordination number 4 !!works only with MESH_3D also defined!!

//model attributes (couplings and fields)
					//without any define -> integer J values
#define J_FIX 		//fix number of J for every interacion -> for fully ferromagnetic or antriferromagnetic (+/-1)
//#define J_GAUSS 	//float J values enabled
					//without any define -> no fixed spins or edges/layers
//#define FIX_SPIN	//fix spin named in the parameter file
//#define FIX_EDGE	//fix edges/layers named in the parameter file
//#define PER_BOUND	//fix one edge/layer on each side/face and copy opposit moving faces in every round, thus making periodic boundary conditions
					//without any define -> no fields added
//#define GLOBAL_FIELD//equal field on every spin defined in the parameter file
//#define LOCAL_FIELD	//identical field on every spin defined in the parameter file

//simulation attributes
					//must choose one
#define METROPOLIS 	//mmc steps evaluated with metripolis dinamics
//#define GLAUBER		//mmc steps evaluated with glauber dinamics

//measurement
//#define MEAS_MAGN	//measure and save magnetisation for every spin, also avarege magnetism and binder parameter
//#define MEAS_CORR	//measure and save correaltion and connected correlation values !!works only with MEAS_MAGN also defined!!
//#define SAVEOUT		//save spin configuration, magnetism and correlation with frequency defiend in the parameter file
//#define GNUPLOT_DATA	//save spin, magnetisation and correlation values also in gnuplot pm3d compatible form

//random generators
//#define RAND_OWN	//define if want to use built-in random generator instead of cuRand
#define INIT_SEED 55731ULL	// Initial seed of spin and interaction initialization
#define CURAND_GEN_TYPE CURAND_RNG_PSEUDO_XORWOW 	// Initial type of cuRand generator