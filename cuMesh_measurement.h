/*
* Measure functions
*
* TODO autocorrelation, spin overlap, link overlap, Binder parameter, Skewness
*/
float meas_sample_avg(thrust::device_ptr<spint> in_begin, thrust::device_ptr<spint> in_end){ //returns average of the values of in array

return(thrust::reduce(in_begin, in_end)/(in_end-in_begin));

}

void meas_time_avg(thrust::device_ptr<spint> in_begin, thrust::device_ptr<spint> in_end, thrust::device_ptr<float> sum_begin){ //sum up -> adds in vector to sum vector, element to element

thrust::transform(in_begin, in_end, sum_begin, sum_begin, thrust::plus<int>());

}

struct saxpy_functor //functor to multiple and add
{
    const float a;

    saxpy_functor(float _a) : a(_a) {}

    __host__ __device__
        float operator()(const float& x, const float& y) const { 
            return a * x + y;
        }
};

void meas_cross_corr(thrust::device_ptr<spint> in_begin, thrust::device_ptr<spint> in_end, thrust::device_ptr<float> sum_begin, int size, int n, int m 
#ifdef MESH_3D
																															, int o
#endif
){ //sum up -> adds in vector to sum vector, element to element

			for( int i = 0; i < m; i++){
				for( int j = 0; j < n; j++){ 
#ifdef MESH_3D
					for(int k = 0; k < o; k++){ 
#endif
					int l = m*i + j
					#ifdef MESH_3D
					* o
#endif
					;
					int size = n * m
#ifdef MESH_3D
					* o
#endif
					;		
					for(int ref = 0; ref < size; ref++){
					
					spint refSpin = in_begin[ref];
					// CC_i[] = CC_i[] + s_i * s_[]
					thrust::transform(in_begin, in_end, sum_begin+l*size, sum_begin+l*size, saxpy_functor(refSpin));
					}
#ifdef MESH_3D		
					}
#endif				
				}
			}

}