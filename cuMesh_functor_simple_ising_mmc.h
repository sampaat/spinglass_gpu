	// Executes metropolis monte carlo sweep on defined sublattice
	// Returns updated spin vector
struct functor_ising_mmc
	{
			int m,n;			//dimensions
			int iter;			//iteration parameter for chessboard selection
			float kT;			//boltzmann factor times temperature, or reciprocal Beta parameter
		#ifdef MESH_3D
			int o;				//z dimension
		#endif
		#ifdef FIX_EDGE
			int fixedEdge;		//number of side layers to fix
		#endif
		#ifdef FIX_SPIN
			int fixSpin;		//spin to fix
		#endif
		#ifdef GLOBAL_FIELD
			float h;			 //external magnetic field
		#endif
		#ifdef J_FIX
			float JF;			 //fix value of J
		#endif


			__host__ __device__
			functor_ising_mmc(int _m,int _n, float _kT, int _iter
		#ifdef MESH_3D
						,int _o
		#endif
		#ifdef FIX_EDGE
						,int _fixedEdge
		#endif
		#ifdef FIX_SPIN
						,int _fixSpin
		#endif
		#ifdef GLOBAL_FIELD
						,float _h
		#endif
		#ifdef J_FIX
						,float _JF
		#endif
						) : m(_m),n(_n),kT(_kT),iter(_iter)
		#ifdef MESH_3D
						,o(_o)
		#endif
		#ifdef FIX_EDGE
						,fixedEdge(_fixedEdge)
		#endif
		#ifdef FIX_SPIN
						,fixSpin(_fixSpin)
		#endif
		#ifdef GLOBAL_FIELD
						,h(_h)
		#endif
		#ifdef J_FIX
						,JF(_JF)
		#endif
						{}
		#ifdef RAND_OWN
				//if we use random numbers generated in the mmc process
			__host__ __device__
				// from numerical recipes TODO but exactly where from :)
			float randUniform(unsigned long long int *ullState) // Ranq1.doub() in NR(3e)
			{
				*ullState ^= *ullState >> 21;
				*ullState ^= *ullState << 35;
				*ullState ^= *ullState >> 4;
				return 5.42101086242752217e-20*((*ullState)*2685821657736338717LL);
			}
		#endif
		
			template <typename Tuple>
			__host__ __device__
				char operator()( Tuple &t)
			{
				//thread properties
				J_str PROP = thrust::get<1>(t); 
				
				//Current id
				int i = PROP.idx;
				//Tray for random number
				float r;
				//Seed define
		#ifdef RAND_OWN
				unsigned long long int s = PROP.seed;
		#else
				curandState s = PROP.seed;
		#endif

				//Current point coordinates
		#ifdef MESH_3D
				int z = i / m / n;
				int y = (i - z * n * m) / m;
				int x = i - z * m * n - y * m;
		#else
				int y = i / m;
				int x = i - y * m;
		#endif
				// chessboard update rule with 2 sub-meshes chosen on the base of "iter"
		#ifdef MESH_HEX
				//in this topology there are 3 separate sub-meshes
				if((x-y)%3 != iter%3)
		#elif defined MESH_3D
				if( ( ( ( (z%2)==0 &&( (x%2)==(y%2) ) ) || ( (z%2)!=0 &&( (x%2)!=(y%2) ) ) ) && (iter%2!=0) ) ||
					( ( ( (z%2)==0 &&( (x%2)!=(y%2) ) ) || ( (z%2)!=0 &&( (x%2)==(y%2) ) ) ) && (iter%2==0) )
				)
		#else
				if(  (((x%2)==(y%2)) && (iter%2==1))  || (((x%2)!=(y%2)) && (iter%2==0)))
		#endif
				{
					return(thrust::get<0>(t)); //if current coordinate is not in the current sub-mesh, then the thread returns with untuched spin value
				}
		#ifdef J_GAUSS
				float dE = 0;
		#else
				int dE = 0;
		#endif
				float expE = 1;

		#ifdef FIX_SPIN
				// do not update fixed spin
				if( i == fixSpin ){

						return(thrust::get<0>(t));
	
				}
		#endif
		#ifdef FIX_EDGE
				// do not count spins in the fixed area
			#ifdef MESH_3D
				if( z<fixedEdge || z >= o-fixedEdge || y<fixedEdge || y >= n-fixedEdge || x<fixedEdge || x >= n-fixedEdge)
			#else
				if( y<fixedEdge || y >= n-fixedEdge || x<fixedEdge || x >= m-fixedEdge)
			#endif
				 {

					return(thrust::get<0>(t));

				}
		#endif
		
				//calculate field of neighboring spins
				// E(i) = - Sum_j{ J(i,j)*S(i)*S(j)) } -h_i*S(i)
				// dE = E(flipSpin(i)) - E(i) = E(-S(i)) - E(S(i)) = 2*(Sum_j{ J(i,j)*S(i)*S(j))} + h_i*s_i)
				
				//rigth neighbor
				
		#if (defined FIX_EDGE && !(defined PER_BOUND))
				if (x < m-1-fixedEdge)
		#else
				if (x < m-1)
		#endif

		#ifdef J_FIX
					dE += JF * thrust::get<2>(t);  // J*s(i+1)
		#else
					dE += PROP.J[0] * thrust::get<2>(t);  // J(i,i+1)*s(i+1)
		#endif
		
				//upper neighbor
				
		#if (defined FIX_EDGE && !(defined PER_BOUND))
				if (y < n-1-fixedEdge)
		#else
				if (y < n-1-0)
		#endif
				{
		#if (defined MESH_TRI || defined MESH_DIA )
				if(x%2==y%2)
		#endif
		#ifdef J_FIX
					dE += JF * thrust::get<3>(t);  // J*s(i+m)
		#else

					dE += PROP.J[1] * thrust::get<3>(t);  // J(i,i+m)*s(i+m)
		#endif
				}

				//left neighbor
				
		#if (defined FIX_EDGE && !(defined PER_BOUND))
				if (x > fixedEdge)
		#else
				if (x > 0)
		#endif
		#ifdef J_FIX
					dE += JF * thrust::get<4>(t);  // J*s(i-1)
		#else
					dE += PROP.J[2] * thrust::get<4>(t);  // J(i,i-1)*s(i-1)
		#endif
		
				//lower neighbor
				
		#if (defined FIX_EDGE && !(defined PER_BOUND))
				if (y > fixedEdge)
		#else
				if (y > 0)
		#endif
				{
		#if (defined MESH_TRI || defined MESH_DIA )
				if(x%2!=y%2)
		#endif
		#ifdef J_FIX
					dE += JF * thrust::get<5>(t);  // J*s(i-m)
		#else
					dE += PROP.J[3] * thrust::get<5>(t);  // J(i,i-m)*s(i-m)
		#endif
				}
		
		
		
		#ifdef MESH_3D
			#ifdef FIX_EDGE
				if (z < o-1-fixedEdge)
			#else
				if (z < o-1)
			#endif
				{
			#ifdef J_FIX
					dE += JF * si(thrust::get<7>(t));  // J*s(i+(m*n))
			#else
					dE += PROP.J[4] * thrust::get<7>(t);  // J(i,i+(m*n))*s(i+(m*n))
				}
			#endif
			#ifdef FIX_EDGE
				if (z > fixedEdge)
			#else
				if (z > 0)
			#endif
				{
			#ifdef J_FIX
					dE += JF * thrust::get<6>(t);  // J*s(i-(m*n))
			#else
					dE += PROP.J[5],thrust::get<6>(t));  // J(i,i-(m*n))*s(i-(m*n))
				}
			#endif
		#endif
		#ifdef MESH_HEX
			#ifdef FIX_EDGE
				if (y > fixedEdge || x < m-1-fixedEdge )
			#else
				if (y > 0 || x < m-1 )
			#endif
				{
			#ifdef J_FIX
					dE += JF * thrust::get<6>(t);  // J*s(i-m+1)
			#else
					dE += PROP.J[4] * thrust::get<6>(t);  // J(i,i-m+1)*s(i-m+1)
				}
			#endif
			#ifdef FIX_EDGE
				if (x > fixedEdge || x < n-1-fixedEdge )
			#else
				if (x > 0 || y < n-1 )
			#endif
				{
			#ifdef J_FIX
					dE += JF * thrust::get<7>(t);  // J*s(i+m-1)
			#else
					dE += PROP.J[5] * thrust::get<7>(t);  // J(i,i+m-1)*s(i+m-1)
				}
		#endif
		#endif
				//calculate additional components from global/local magnetic fields
		#ifdef GLOBAL_FIELD
				dE += h * thrust::get<0>(t);
		#elif defined LOCAL_FIELD
				dE += PROP.h * thrust::get<0>(t);
		#endif
				//final calculation of dE
				dE *= 2 * thrust::get<0>(t);
				
		#if !( defined J_GAUSS || defined GLOBAL_FIELD || defined LOCAL_FIELD)
				if( dE > 0 ){
					expE = PROP.E[int(dE/2-1)];
				}else{
				if( dE < 0)
					expE = 1/PROP.E[int(-dE/2-1)];
				}
				
		#else
				expE = exp(dE / kT);
		#endif

		#ifdef METROPOLIS
			#ifdef RAND_OWN
				r = randUniform(&(s)); //get random number from function
			#else
				r = curand_uniform(&(s));//get precalculated random number
			#endif
				PROP.seed = s; // update the seed before returning

				if (1/expE > r) { //flip or do not flip the spin
					return(-thrust::get<0>(t));
				}else{
					return(thrust::get<0>(t));
				}
		#else
				// Glauber dynamics
			#ifdef RAND_OWN
				r = randUniform(&(s)); //get random number from function
			#else
				r = curand_uniform(&(s)); //get precalculated random number
			#endif
				PROP.seed = s; // update the seed before returning
				if( r<1/(1+expE) ) {
					return(-thrust::get<0>(t));
				}else{
					return(thrust::get<0>(t));
				}
		#endif

		}


};
