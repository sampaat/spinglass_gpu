/*
 * Sweep functions
 *
 * TODO Sweep for local energy, and energy count
 */

	// Executes metropolis monte carlo sweep on defined sublattice
	// Uses MMC_IsingspinFunctor
void subsweep_ising_metropolis(thrust::device_vector<spint>& in,thrust::device_vector<spint>& out, thrust::device_vector<J_str>& J, float kT, int iter, int m, int n
#ifdef MESH_3D
	,int o				//z dimension
#endif
#ifdef FIX_EDGE
	,int fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,int fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,float h			 //external magnetic field
#endif
#ifdef J_FIX
	,float JF			 //fix value of J
#endif
)
{
	// note: make_tuple can handle only max 10 elements
	// put spin index and seed into a struct


	thrust::transform(
		thrust::make_zip_iterator(  // starting indicies
				thrust::make_tuple(
						in.begin(),     //S(x,y(,z))=S(i)
						J.begin(),      // J_struct
						in.begin() + 1, //S(x+1,y(,z)) = S(i+1)
						in.begin() + m, //S(x,y+1(,z)) = S(i+m)
						in.begin() - 1, //S(x-1,y(,z)) = S(i-1)
						in.begin() - m //S(x,y+1(,z)) = S(i-m)
#ifdef MESH_3D
						,in.begin() - m*n //S(x,y,z-1) = S(i-m*n)
						,in.begin() + m*n //S(x,y,z+1) = S(i+m*n)
#endif
#ifdef MESH_HEX
						,in.begin() - m+1 //S(x+1,y-1) = S(i-m+1)
						,in.begin() + m-1 //S(x-1,y+1) = S(i+m-1)
#endif

				)
		//thrust::counting_iterator<int>(0))
		),
		thrust::make_zip_iterator(  // ending indicies
				thrust::make_tuple(
						in.end(),     //S(x,y) = S(i) // note the x,y indices are rows and columns not columns and rows like in standard C matrix notation
						J.end(),      // J(i,i+1) == J({x,y},{x+1,y}), J(i,i+m), J(i,i-1), J(i,i-m) // see above the spin indecies
						in.end() + 1, //S(x+1,y) = S(i+1)
						in.end() + m, //S(x,y+1) = S(i+m)
						in.end() - 1, //S(x-1,y) = S(i-1)
						in.end() - m //S(x,y+1) = S(i+m)
#ifdef MESH_3D
						,in.end() - m*n //S(x,y,z-1) = S(i-m*n)
						,in.end() + m*n //S(x,y,z+1) = S(i+m*n)
#endif
#ifdef MESH_HEX
						,in.end() - m+1 //S(x+1,y-1) = S(i-m+1)
						,in.end() + m-1 //S(x-1,y+1) = S(i+m-1)
#endif
				)
		)
		,
		out.begin(), // indicies for ouput
		functor_ising_mmc(m,n,kT,iter
#ifdef MESH_3D
				,o				//z dimension
#endif
#ifdef FIX_EDGE
				,fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
				,fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
				,h			 //external magnetic field
#endif
#ifdef J_FIX
				,JF			 //fix value of J
#endif
		) // functor to use
	);
}

	//Executes mmc sweeps on al the submeshes
 
void sweep_2sub_ising_metropolis(thrust::device_vector<spint>& spin1,thrust::device_vector<spint>& spin2, thrust::device_vector<J_str>& J, float kT, int m, int n
#ifdef FIX_EDGE
	,int fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,int fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,float h			 //external magnetic field
#endif
#ifdef J_FIX
	,float JF			 //fix value of J
#endif
){
//sublattice 1
subsweep_ising_metropolis(spin1, spin2, J, kT, 0, m, n
#ifdef FIX_EDGE
	, fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	, fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	, h			 //external magnetic field
#endif
#ifdef J_FIX
	, JF			 //fix value of J
#endif
);
//sublattice 2
subsweep_ising_metropolis(spin2, spin1, J, kT, 1, m, n
#ifdef FIX_EDGE
	,fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,h			 //external magnetic field
#endif
#ifdef J_FIX
	,JF			 //fix value of J
#endif
);
}

void sweep_3sub_ising_metropolis(thrust::device_vector<spint>& spin1, thrust::device_vector<spint>& spin2, thrust::device_vector<spint>& spin3, thrust::device_vector<J_str>& J, float kT, int m, int n
#ifdef MESH_3D
	,o				//z dimension
#endif
#ifdef FIX_EDGE
	,int fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,int fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,float h			 //external magnetic field
#endif
#ifdef J_FIX
	,float JF			 //fix value of J
#endif
){
//sublattice 1
subsweep_ising_metropolis(spin1,spin2, J, kT, 0, m, n
#ifdef MESH_3D
	,o				//z dimension
#endif
#ifdef FIX_EDGE
	,fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,h			 //external magnetic field
#endif
#ifdef J_FIX
	,JF			 //fix value of J
#endif
);
//sublattice 2
subsweep_ising_metropolis(spin2,spin3, J, kT, 1, m, n
#ifdef MESH_3D
	,o				//z dimension
#endif
#ifdef FIX_EDGE
	,fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,h			 //external magnetic field
#endif
#ifdef J_FIX
	,JF			 //fix value of J
#endif
);
//sublattice 3
subsweep_ising_metropolis(spin3, spin1, J, kT, 2, m, n
#ifdef MESH_3D
	,o				//z dimension
#endif
#ifdef FIX_EDGE
	,fixedEdge		//number of side layers to fix
#endif
#ifdef FIX_SPIN
	,fixSpin		//spin to fix
#endif
#ifdef GLOBAL_FIELD
	,h			 //external magnetic field
#endif
#ifdef J_FIX
	,JF			 //fix value of J
#endif
);
}