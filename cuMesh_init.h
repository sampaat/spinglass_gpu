/*
 * Initialization functions
 */

// Heating up random generator
curandGenerator_t heatup_cuRand(curandGenerator_t gen, int time){
	int size = 1000;
	float *randTray;
	cudaMalloc((void**)&randTray, size*sizeof(float));
	for(int i=0; i<(time%527); i++){
	curandGenerateUniform(gen, randTray, size);
	}
	return gen;
}

// Initialize id part of interaction matrix
void init_J_idx(thrust::host_vector<J_str>& v,int size){
	for(unsigned int i = 0 ; i < size; i++) {
		v[i].idx = i;
	}
}

struct init_curand_functor
{
	int s;
			
	init_curand_functor(int _s) : s(_s) {}
	__host__ __device__
		void operator()(J_str &v)  {
			curand_init(s, v.idx, 0, &v.seed);
	}
			
};

// Initialize seeds for mmc
void init_seed_thread(thrust::host_vector<J_str>& v, int size
#ifdef RAND_OWN
	, curandGenerator_t gen
#endif
){
#ifdef RAND_OWN
	float *randTray;
	cudaMalloc((void**)&randTray, size*sizeof(float));
	curandGenerateUniform(gen, randTray, size);
	for(int i = 0 ; i < size; i++) {
		v[i].seed = unsigned long long int(1307737860*randTray[i]);
	}
#else
	thrust::transform(v.begin(),v.end(),v.begin(),init_curand_functor(INIT_SEED));
#endif
}


// Gives a spin value to a uniform [0,1[ float with proportion based on "a"
struct spingen_functor
{
	const float a;
	
	spingen_functor(float _a) : a(_a) {}
	
	__host__ __device__
		char operator()( float& x)  {
			return( char( x + a ) * 2 - 1 );
	}
};

	// Initialize spin vector with random boolean (treated as +/-1)
void init_ising_spin(thrust::device_vector<spint>& v, curandGenerator_t gen, int size)
{
	float *randTray;
	thrust::device_ptr<float> randTray_ptr;
	cudaMalloc((void**)&randTray, size*sizeof(float));
	curandGenerateUniform(gen, randTray, size);
	randTray_ptr = thrust::device_pointer_cast(randTray);
	thrust::transform(randTray_ptr, randTray_ptr+size, v.begin(), spingen_functor(0.0f));
}

#ifndef J_FIX
	// Initialize J_ij vector with random
void init_ising_J_ij(thrust::host_vector<J_str>& v, float a, curandGenerator_t gen, int size, int m, int n
#ifdef MESH_3D
					,int o
#endif
){
	a -= floor(a); //just to be sure, that a in [0,1[
	float randTray;
	thrust::device_ptr<float> randTray_ptr;
	cudaMalloc((void **)&randTray, 6*size);

	//generate random couplings of desired type
#ifdef J_GAUSS
	curandGenerateNormal(gen, randTray, 6*size, float 0.0, float 1.0);
#else
	curandGenerateUniform(gen, &randTray, 6*size);
	thrust::transform(randTray_ptr, randTray_ptr+size, randTray_ptr, spingen_functor(a));
	//TODO check if works
#endif

	// we have to go in a chessboard order and set all the 4 J_IJ to be symmetric

	for(int i=0; i<size; i++){
#ifdef MESH_3D
		int z = i / m / n;
		int y = (i - z * n * m) / m;
		int x = i - z * m * n - y * m;
#else
		int y = i / m;
		int x = i - y * m;
#endif
		//only set for one of the sub-meshes ("chessboard" logic)
#ifdef MESH_HEX
		//in this topology there are 3 separate sub-meshes
		if((x-y)%3 != 0)
#elif defined MESH_3D
		//in this topology there are 2 separate sub-meshes
		if ( ( (z%2)==0 &&( (x%2)==(y%2) ) ) || ( (z%2)!=0 &&( (x%2)!=(y%2) ) ) )
#else
		//in this topology there are 2 separate sub-meshes
		if ((x%2)==(y%2))
#endif
		{
				if( x < m-1 ){
					v[i].J[0] = v[i+1].J[2] = randTray[i]; // J(i,i+1) == J(i+1,i)
				}
				if( y < n-1){
					v[i].J[1] = v[i+m].J[3] = randTray[i+size]; // J(i,i+m) == J(i+m,i)
				}
				if( x > 0 ){
					v[i].J[2] = v[i-1].J[0] = randTray[i+2*size]; // J(i,i-1) == J(i-1,i)
				}
				if( y > 0 ){
					v[i].J[3] = v[i-m].J[1] = randTray[i+3*size]; // J(i,i-m) == J(i-m,i)
				}
#ifdef MESH_HEX
				if(y > 0 || x < m-1 ){
					v[i].J[4] = v[i-m+1].J[5] = randTray[i+4*size]; // J(i,i-m+1) == J(i-m+1,i)
				}
				if(x > 0 || y < n-1 ){
					v[i].J[5] = v[i+m-1].J[4] = randTray[i+5*size]; // J(i,i+m-1) == J(i+m-1,i)
				}
#elif defined MESH_3D
				if(z < o-1 ){
					v[i].J[4] = v[i+(m*n)].J[5] = randTray[i+4*size]; // J(i,i+(m*n)) == J(i+(m*n),i)
				}
				if(z > 0 ){
					v[i].J[5] = v[i-(m*n)].J[4] = randTray[i+5*size]; // J(i,i-(m*n)) == J(i-(m*n),i)
				}
#endif
		}
	}
	
#ifdef PER_BOUND
	//set simetric couplings on the closed sides (mind the fixed edges!)
	
	#ifdef MESH_3D
	
	int mn = m * n;
	
	for(int j = 1; j < ( o - 1 ); j++){
		for(int i = 1; i < ( m - 1 ); i++){
			v[j*mn+i+((n-2)*m)].J[1] = v[j*mn+i].J[3];			
		}
		for(int i = 1; i < ( n - 1 ); i++){
			v[j*mn+m-2+(i*m)].J[0] = v[j*mn+1+(i*m)].J[2];			
		}
	}
	for(int i = 1; i < ( m - 1 ); i++){
		for(int j = 1; j < ( n - 1 ); j++){
			v[mn+(j*m)+i].J[5] = v[(o-2)*mn+(j*m)+i].J[4];	
		}
	}
	#else
		for(int i = 1; i < ( m - 1 ); i++){
			v[i+((n-2)*m)].J[1] = v[i].J[3];			
		}
		for(int i = 1; i < ( n - 1 ); i++){
			v[m-2+(i*m)].J[0] = v[1+(i*m)].J[2];			
		}
		#ifdef MESH_HEX
		v[(n-2)*m+1].J[5] = v[2*m-2].J[4];
		#endif
	#endif
	
#endif
}
#endif

#ifdef PER_BOUND
	// Set periodic boundaries on spin vector
void per_bound_copy(thrust::device_vector<spint>& in, int m, int n
#ifdef MESH_3D
	, int o
#endif
		){
#ifdef MESH_3D
	int mn = m * n;
	
	for(int j = 1; j < ( o - 2 ); j++ ){
		for(int i = 1; i < ( n - 2 ); i++ ){
			thrust::copy_n(in.begin()+j*mn+i*m+1, 1, in.begin()+j*mn+(i*m+m); //creating side edges left to right
			thrust::copy_n(in.begin()+j*mn+i*m+(m-1), 1, in.begin()+j*mn+i*m); //creating side edges right to left
		}
		thrust::copy_n(in.begin()+j*mn+m+1, m-2, in.begin()+j*mn+((n-1)*m)+1);	//creating top and bottom edges
		thrust::copy_n(in.begin()+j*mn+((n-2)*m)+1, m-2, in.begin()+j*mn+1);
	}
		thrust::copy_n(in.begin()+mn, mn, in.begin()+(o-1)*mn);	//creating front and rear faces
		thrust::copy_n(in.begin()+(o-1)*mn, mn, in.begin());
#else
	for(int i = 1; i < ( n - 2 ); i++){
		thrust::copy_n(in.begin()+i*m+1, 1, in.begin()+(i*m+m)); //creating side edges left to right
		thrust::copy_n(in.begin()+i*m+(m-1), 1, in.begin()+i*m); //creating side edges right to left
	}
	thrust::copy_n(in.begin()+m+1, m-2, in.begin()+((n-1)*m)+1);	//creating top and bottom edges
	thrust::copy_n(in.begin()+((n-2)*m)+1, m-2, in.begin()+1);

	#ifdef MESH_HEX	
	thrust::copy_n(in.begin()+m+1, 1, in.begin()+(n*m)-1);	//creating corners
	thrust::copy_n(in.begin()+m+(m-2), 1, in.begin()+((n-1)*m);
	thrust::copy_n(in.begin()+((n-2)*m)+1, 1, in.begin());
	thrust::copy_n(in.begin()+(n*m)-m-2, 1, in.begin()+m-1);
	#endif
#endif
}
#endif

#if !( defined J_GAUSS || defined GLOBAL_FIELD || defined LOCAL_FIELD)
	// Calculate exponentials to help energy calcucation
void init_expE(thrust::host_vector<J_str>& v, float kT, int size){


	float e2 = exp(2/kT);
	float e4 = exp(4/kT);
	float e6 = exp(6/kT);
	float e8 = exp(8/kT);
	float e10 = exp(10/kT);
	float e12 = exp(12/kT);
	for(int i =0; i < size; i++){

	v[i].E[0] = e2;
	v[i].E[1] = e4;
	v[i].E[2] = e6;
	v[i].E[3] = e8;
	#if ( defined MESH_3D || defined MESH_HEX)
	v[i].E[4] = e10;	
	v[i].E[5] = e12;
	#endif
	}

}
#endif