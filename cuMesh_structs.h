/*
 * Spin and interaction (J) structures
 */

typedef signed char spint;
 
typedef struct {
		unsigned int idx;
#ifdef J_FIX

#elif defined J_GAUSS
	#if ( defined MESH_3D || defined MESH_HEX)
		float J[6];
	#else
		float J[4];
	#endif
#else
	#if ( defined MESH_3D || defined MESH_HEX)
		signed char J[6];
	#else
		signed char J[4];
	#endif
#endif
#ifdef RAND_OWN
	unsigned long long int seed;
#else
	curandState seed;
#endif	
#ifdef LOCAL_FIELD
		float h;
#endif
#if !( defined J_GAUSS || defined GLOBAL_FIELD || defined LOCAL_FIELD)
	#if ( defined MESH_3D || defined MESH_HEX)
		float E[6];
	#else
		float E[4];
	#endif
#endif
} J_str;