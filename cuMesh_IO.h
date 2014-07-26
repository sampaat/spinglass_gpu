/*
 * IO functions
 */
 
//loads a single parameter (paramName) from pFile
string load_Parameter(string paramName, FILE * pFile){
	int paramSize_1 = paramName.size()-1;
	string readName = paramName; readName[0] = '0';
	string paramValue = "";
	char tmp [2];

	while (readName != paramName)
	{
		for (int i=0; i<paramSize_1; i++)
			readName[i] = readName[i+1];

		fgets(tmp,2,pFile);

		if (tmp[0] == '%')
		{
			cerr << "Error: Parameter named \"" + paramName + "\" not found." << endl;
			exit (EXIT_FAILURE);
		}

		readName[paramSize_1] = tmp[0];
	}

	fgets(tmp,2,pFile);

	if (tmp[0] == '\n')
	{
		return paramValue = "";
	}
	else
	{
		while ((tmp[0] == ' ') || (tmp[0] == '\t'))
			fgets(tmp,2,pFile);

		while ((tmp[0] != ' ') && (tmp[0] != '\t') && (tmp[0] != '\n') && (tmp[0] != ',') && (tmp[0] != ';'))
		{
			paramValue += tmp;
			fgets(tmp,2,pFile);
		}
	
		tmp[0] = ' ';

		return paramValue;
	}
}

//loads the file containing the initial spin states
void load_Svector_int(string Sfname, thrust::host_vector<spint> spin_host, int size){
		cout << "Loading initial spin configuration from file" << endl; 
		ifstream Sfile;
		Sfile.open(Sfname.c_str());
		string spin;
		if (Sfile.is_open())
			{
				getline (Sfile,spin);
				cout << spin << endl; //cout information line of sfile
				for(int i=0;i<size;i++){
					getline (Sfile,spin);
					if(spin == "-1"){
						spin_host[i] = char(-1);
					}else{
						if(spin == "1")
						cout << "Error: not +/- 1 spin state, autouse S[i]=1" << endl;
						spin_host[i] = 1;
					}
				}
			Sfile.close();
		}

		else{
			cout << "Unable to open file" << endl;
			exit (EXIT_FAILURE);
		}
}

#ifndef J_FIX
//loads the file containing the J values
void load_Jvector_int(string Jfname, thrust::host_vector<J_str> J_host, int size){
		cout << "Loading initial J configuration from file" << endl; 
		ifstream Jfile (Jfname.c_str());
		string line, each;
		int links = 
#if ( defined MESH_3D || defined MESH_HEX)
					6;
#else
					4;
#endif

		if (Jfile.is_open()){
				getline (Jfile,line);
				cout << line << endl; //cout information line of sfile
				for(int i=0;i<size;i++){
					getline (Jfile,line);
					istringstream split(line);
					for(int j=0; j<links; j++){
						getline(line, each, "\t");
							if(each == "1" || each == "-1"){
								CUDA_CALL(J_host[i].J[j] = atoi(each.c_str()));
							}else{
								cout << "Error: not +/- 1 link state, autouse J[i][j]=1" << endl;
								CUDA_CALL(J_host[i].J[j] = 1);
							}
					}
				}
			Jfile.close();
		}else{
			cout << "Unable to open file" << endl; 
			exit (EXIT_FAILURE);
		}
}
#endif

#ifdef LOCAL_FIELD
//loads the file containing the initial H values
void load_Hvector(string Hfname, thrust::host_vector<J_str> J_host, int size){
		cout << "Loading initial magnetic field configuration from file" << endl; 
		ifstream Hfile;
		Hfile.open(Hfname.c_str());
		string line;
		if (Hfile.is_open()){
				getline (Hfile,line);
				cout << line << endl; //cout information line of sfile
				for(int i=0;i<size;i++){
					getline (Hfile,line);
					CUDA_CALL(J_host[i].h = atof(line.c_str()));
				}
			Hfile.close();
		}else{
			cout << "Unable to open file" << endl;
			exit (EXIT_FAILURE);
		}
}
#endif

//saves the spin states
void saveout_spin(stringstream infoline, thrust::host_vector<spint> spin, int start, int t, int size, float T, int n, int m
#ifdef MESH_3D
			, int o
#endif
		){
	cout << "Saving spin configuration to file" << endl; 
	stringstream filename;
	filename << "spinstate_start_" << start << "_t_" << t + "_T_" << T << "_dim_" << m << "_" << n;
#ifdef MESH_3D
	filename << "_" << o;
#endif
	filename << ".dat";
		
	ofstream Sfile;
		
	Sfile.open(filename.str().c_str());
	
	if (Sfile.is_open()){
		Sfile << infoline << endl;
		
		for(int i = 0; i < size; i++){
		
			Sfile << spin[i] << endl;
			
		}
		Sfile.close();
	}
	else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}
		
		
#ifdef GNUPLOT_DATA
		
	stringstream filename << "spinstate_start_" << start << "_t_" << t << "_T_" << T << "_dim_" << m << "_" << n; 
	#ifdef MESH_3D
	filename  << "_" << o;
	#endif
	filename  << ".gpdat";
		
		ofstream Sfile;
		
		Sfile.open(filename);
		
	if (Sfile.is_open()){
		Sfile << infoline << endl;
		
	#ifdef MESH_3D
		for(int i = 0; i < m; i++ ){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < o; k++){
					Sfile << i << " " << j << " " << k << " " << spin[k*o + j*n + i] << endl;
				}
				Sfile << endl;
			}
			Sfile << endl;
		}
	#else
		for(int i = 0; i < m; i++ ){
			for(int j = 0; j < n; j++){
				Sfile << i << " " << j << " " << spin[k*o + j*n + i] << endl;
			}
			Sfile << endl;
		}		
	#endif	
		Sfile.close();
	}
	else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}		
#endif
}

#ifndef J_FIX
//saves the J values
void saveout_J(stringstream infoline, thrust::host_vector<J_str> J_host, int size, int n, int m 
#ifdef MESH_3D
			, int o
#endif
		){
		cout << "Saving J configuration to file" << endl; 
		
		stringstream filename << "Jfile_dim_" << m << "_" << n;
#ifdef MESH_3D
		filename  << "_" << o;
#endif
		filename  << ".dat";
		ifstream Jfile (filename);
		int links = 
#if ( defined MESH_3D || defined MESH_HEX)
					6;
#else
					4;
#endif

		if (Jfile.is_open())
			{
				Jfile << infoline << endl;
				for(int i = 0; i < size; i++){
					for(int j = 0; j < links; j++){
						Jfile << J_host[i].J[j] << "\t";
					}
					Jfile << endl;
				}
			Jfile.close();
		}
		else{
			cout << "Unable to open file" << endl; 
			exit (EXIT_FAILURE);
		}
}
#endif

#ifdef LOCAL_FIELD
//saves the H values
void saveout_H(stringstream infoline; thrust::host_vector<J_str> J_host, int size, int n, int m 
#ifdef MESH_3D
			, int o
#endif
		){
			cout << "Saving H configuration to file" << endl; 
		
		stringstream filename << "Hfile_dim_" << m << "_" << n;
#ifdef MESH_3D
		filename  << "_" << o;
#endif
		filename  << ".dat";
		
	ofstream Hfile
		
	Hfile.open(filename);
	
	if (Hfile.is_open()){
		Hfile << infoline << endl;
		
		for(int i = 0; i < size; i++){
		
			Hfile << J_host[i].J[j] << endl;
			
			Hfile.close;
		}
	}
	else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}
}
#endif

#ifdef MEAS_MAGN
//saves the magnetization and correlation matrices
void saveout_magn(stringstream infoline; thrust::host_vector<float>sumSpin_host, int start, int t, int f_meas, int size, float T, int n, int m 
			#ifdef MESH_3D
			, int o
			#endif
	){
	cout << "Saving magnetisation values to file" << endl; 
	stringstteam filename << "magnetisation_start_" << start << "_t_" << t << "_T_" << T << "_dim_" << m << "_" << n;
#ifdef MESH_3D
	filename  << "_" << o;
#endif
	filename  << ".dat";
		
	ofstream Mfile;
		
	Mfile.open(filename);
	
	if (Mfile.is_open()){
		Mfile << infoline << endl;
		
			for( int i = 0; i < m; i++){
				for( int j = 0; j < n; j++){
#ifdef MESH_3D
					for(int k = 0; k < o; k++){
#endif
						int l =
#ifdef MESH_3D
						k*m*n
#endif					
						j*m+i;

						Mfile << i << " " << sumSpin_host[k]/(t/f_meas) << endl;
#ifdef MESH_3D
						}
#endif
					}
				}

	}else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}
	
#ifdef GNUPLOT_DATA
	stringstream filename << "magnetism_start_" << start << "_t_" << t << "_T_" << T << "_dim_" << m << "_" << n;
#ifdef MESH_3D
	filename  << "_" << o;
#endif
	filename  << ".gpdat";
		
	Mfile.open(filename);
	
	if (Mfile.is_open()){
		Mfile << infoline << endl;
		
			for( int i = 0; i < m; i++){
				for( int j = 0; j < n; j++){
#ifdef MESH_3D
					for(int k = 0; k < o; k++){
#endif
						int l =
#ifdef MESH_3D
						k*m*n
#endif					
						j*m+i;

						Mfile << i << " " << j << " ";
#ifdef MESH_3D
						Mfile << k << " ";
#endif
						Mfile << sumSpin_host[k]/(t/f_meas) << endl;
#ifdef MESH_3D
						}
						Mfile <<  endl;
#endif
					}
				Mfile <<  endl;
				}

	}else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}
#endif
}
#endif

#ifdef MEAS_CORR
//saves the correlation values
void saveout_corr(stringstream infoline; thrust::host_vector<float>sumSpin_host, thrust::host_vector<float>corrSpin_host, int start, int t, int f_meas, int size, float T, int n, int m 
			#ifdef MESH_3D
			, int o
			#endif
	){
	cout << "Saving correlation values to file" << endl; 
	stringstream filename << "correlation_start_" << start << "_t_" << t << "_T_" << T << "_dim_" << m << "_" << n; 
#ifdef MESH_3D
	filename  << "_" << o;
#endif
	filename  << ".dat";
		
	ofstream Cfile;
		
	Cfile.open(filename);
	
	if (Cfile.is_open()){
		Cfile << infoline << endl;
		
			for( int i = 0; i < m; i++){
				for( int j = 0; j < n; j++){
#ifdef MESH_3D
					for(int k = 0; k < o; k++){
#endif
						int l =
#ifdef MESH_3D
						k*m*n
#endif					
						j*m+i;
						for( int ii = 0; ii < i; ii++){
							for( int jj = 0; jj < j; jj++){
#ifdef MESH_3D
								for(int kk = 0; kk < k; k++){
#endif						
									int ll =
#ifdef MESH_3D
									kk*m*n
#endif					
									jj*m+ii;
#ifdef MESH_3D
									float c0 = corrSpin_host[k*m*n*o+kk]/(t/f_meas);
#else
									float c0 = corrSpin_host[k*mn+kk]/(t/f_meas);
#endif	

									float c1 = (sumSpin_host[k]*sumSpin_host[kk])/(t/f_meas);
								Cfile << k << " " << kk << " " << c0 << " " <<  c0-c1 <<endl;
#ifdef MESH_3D
								}
#endif
							}
						}
#ifdef MESH_3D
					}
#endif
				}
			}
	}
	
	else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}
#ifdef GNUPLOT_DATA
	stringstream filename << "correlation_start_" << start << "_t_" << t << "_T_" << T << "_dim_" << m << "_" << n;
#ifdef MESH_3D
	filename  << "_" << o;
#endif
	filename  << ".gpdat";
		
	Cfile.open(filename);
	
	if (Cfile.is_open()){
		Cfile << infoline << endl;
		
			for( int i = 0; i < m; i++){
				for( int j = 0; j < n; j++){
#ifdef MESH_3D
					for(int k = 0; k < o; k++){
#endif
						int l =
#ifdef MESH_3D
						k*m*n
#endif					
						j*m+i;
						
						Cfile << i << " " << j;
#ifdef MESH_3D					
						Cfile << " " << k;
#endif						
						
						for( int ii = 0; ii < i; ii++){
							for( int jj = 0; jj < j; jj++){
#ifdef MESH_3D
								for(int kk = 0; kk < k; k++){
#endif						
									int ll =
#ifdef MESH_3D
									kk*m*n
#endif					
									jj*m+ii;
#ifdef MESH_3D
									float c0 = corrSpin_host[kk*m*n*o+k];
#else
									float c0 = corrSpin_host[kk*mn+k];
#endif	

									float c1 = sumSpin_host[k]*sumSpin_host[kk];
									Cfile << " " << c0 << " " <<  c0-c1;
#ifdef MESH_3D
								}
#endif
							}
						}
						<< endl;
#ifdef MESH_3D
					}
				<< endl;
#endif
				}
			<< endl;
			}
	}
	
	else{
		cout << "Unable to open file" << endl; 
		exit (EXIT_FAILURE);
	}
#endif
}
#endif

