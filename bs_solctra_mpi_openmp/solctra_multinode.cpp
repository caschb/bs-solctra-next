#include "solctra_multinode.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <mpi.h>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_num_threads() 1
#endif


void initializeGlobals(Coil* rmi, Coil* rmf)
{
    for(unsigned int i = 0 ; i < TOTAL_OF_COILS ; ++i)
    {
        rmi[i].x = static_cast<double*>(_mmm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmi[i].y = static_cast<double*>(_mmm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmi[i].z = static_cast<double*>(_mmm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmf[i].x = static_cast<double*>(_mmm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmf[i].y = static_cast<double*>(_mmm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
        rmf[i].z = static_cast<double*>(_mmm_malloc(sizeof(double) * (TOTAL_OF_GRADES + 1), ALIGNMENT_SIZE));
    }
}

void finishGlobal(Coil* rmi, Coil* rmf)
{
    for(unsigned int i = 0 ; i < TOTAL_OF_COILS ; ++i)
    {
        _mmm_free(rmi[i].x);
        _mmm_free(rmi[i].y);
        _mmm_free(rmi[i].z);
        _mmm_free(rmf[i].x);
        _mmm_free(rmf[i].y);
        _mmm_free(rmf[i].z);
    }
}
void load_coil_data(double* x, double* y, double* z, const std::string& path)
{
    for (int num = 0; num < TOTAL_OF_COILS; num++)
    {
	std::ostringstream convert;
	convert << num;
	std::string value = convert.str();
	std::string tmp =  path + "/Bobina"+value+"m.txt";
        loadCartesianFile(&(x[num * TOTAL_OF_GRADES_PADDED]), &(y[num * TOTAL_OF_GRADES_PADDED]), &(z[num * TOTAL_OF_GRADES_PADDED]), TOTAL_OF_GRADES + 1, tmp);
    }

}
void e_roof(GlobalData& data)
{
    cartesian segment;
    for (int j = 0; j < TOTAL_OF_COILS; j++)
    {
        const int base = j * TOTAL_OF_GRADES_PADDED;
        #pragma GCC ivdep
    	for (int i = 0; i < TOTAL_OF_GRADES; i++)
            {
                segment.x = ( data.coils.x[base + i + 1] ) - ( data.coils.x[base + i] );
                segment.y = ( data.coils.y[base + i + 1] ) - ( data.coils.y[base + i] );
                segment.z = ( data.coils.z[base + i + 1] ) - ( data.coils.z[base + i] );
                data.leng_segment[base + i] = norm_of(segment);
                const double leng_segment_inverted = 1.0 / data.leng_segment[base + i];
                data.e_roof.x[base + i] = segment.x * leng_segment_inverted;
                data.e_roof.y[base + i] = segment.y * leng_segment_inverted;
                data.e_roof.z[base + i] = segment.z * leng_segment_inverted;
            }
    }
}

void printeroof(GlobalData& data, const int subsetIndex, const std::string& outputPath){                                                                              
        std::ostringstream indexString;                                         
        indexString << subsetIndex;                                                
        std::string valueIndex = indexString.str();                                
        
        FILE* handler;                                                             
        std::string file_name = outputPath +  "/eroof" + valueIndex + ".txt";   
        handler = fopen(file_name.c_str(), "a");                                   
        if(handler == NULL)                                                        
        {                                                                          
            printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
            exit(0);                                                               
        }                                                                       
                                                                                   
                                                                                    
        for (int j = 0; j < TOTAL_OF_COILS; j++){                                        
            const int base = j * TOTAL_OF_GRADES_PADDED;                            
            for (int i = 0; i < TOTAL_OF_GRADES; i++)                                   
            {                                                                      
                fprintf(handler, "%.17g\t%.17g\t%.17g\n", data.e_roof.x[base+i], data.e_roof.y[base + i], data.e_roof.z[base + i]);    
            }                                                                      
        }                                                                            
                                                                                   
        fclose(handler);                                                           
                                                                                   
}      

cartesian magnetic_field(Coil* rmi, Coil* rmf, const GlobalData& data, const Particle& point)
{
    const double multiplier = ( miu * I ) / ( 4 * PI );

    double Bx = 0;
    double By = 0;
    double Bz = 0;

    Coil* rmiA = (Coil*)__builtin_assume_aligned(rmi, 64);
    Coil* rmfA = (Coil*)__builtin_assume_aligned(rmf, 64);
    
    for (unsigned i = 0; i < TOTAL_OF_COILS; i++)
    {
        for (unsigned jj = 0; jj < TOTAL_OF_GRADES; jj += GRADES_PER_PAGE)
        {
            const unsigned final = (TOTAL_OF_GRADES < jj + GRADES_PER_PAGE) ? TOTAL_OF_GRADES : jj + GRADES_PER_PAGE;
            const int base = i * TOTAL_OF_GRADES_PADDED;
            double* x = &data.coils.x[base];
            double* y = &data.coils.y[base];
            double* z = &data.coils.z[base];

 
            #pragma omp simd aligned(rmiA,rmfA:64)
            for (unsigned j = jj; j < final ; ++j)
            {
                rmiA[i].x[j] = point.x - x[j];
                rmiA[i].y[j] = point.y - y[j];
                rmiA[i].z[j] = point.z - z[j];
                rmfA[i].x[j] = point.x - x[j + 1];
                rmfA[i].y[j] = point.y - y[j + 1];
                rmfA[i].z[j] = point.z - z[j + 1];
            }

            #pragma omp simd aligned(rmiA,rmfA:64) reduction(+:Bx) reduction(-:By) reduction(+:Bz)
            for (unsigned j = jj; j < final ; ++j)
            {
               const double norm_Rmi = sqrt((( rmiA[i].x[j] * rmiA[i].x[j] ) + ( rmiA[i].y[j] * rmiA[i].y[j] ) +
                                              ( rmiA[i].z[j] * rmiA[i].z[j] )));
               const double norm_Rmf = sqrt((( rmfA[i].x[j] * rmfA[i].x[j] ) + ( rmfA[i].y[j] * rmfA[i].y[j] ) +
                                              ( rmfA[i].z[j] * rmfA[i].z[j] )));

                //firts vector of cross product in equation 8
                cartesian U;
                U.x = multiplier * data.e_roof.x[base + j];
                U.y = multiplier * data.e_roof.y[base + j];
                U.z = multiplier * data.e_roof.z[base + j];

                //second vector of cross product in equation 8
                const double C = (
                        (( 2 * ( data.leng_segment[base + j] ) * ( norm_Rmi + norm_Rmf )) /
                          ( norm_Rmi * norm_Rmf )) *
                         (( 1 ) / (( norm_Rmi + norm_Rmf ) * ( norm_Rmi + norm_Rmf ) -
                                  data.leng_segment[base + j] * data.leng_segment[base + j] )));

                //printf("C value: %.17g\n",C);

                cartesian V;
                V.x = rmiA[i].x[j] * C;
                V.y = rmiA[i].y[j] * C;
                V.z = rmiA[i].z[j] * C;

                //cross product in equation 8
                Bx = Bx + (( U.y * V.z ) - ( U.z * V.y ));
                By = By - (( U.x * V.z ) - ( U.z * V.x ));
                Bz = Bz + (( U.x * V.y ) - ( U.y * V.x ));
            }
        }
    }

    cartesian B = {0.0, 0.0, 0.0};
    B.x = Bx;
    B.y = By;
    B.z = Bz;

    return B;
}


void printParallelIterationFile(const Particle* particle_array, const int iteration, const std::string& output, const int rank_id, const int length, const int offset){

	std::ostringstream convert;
    convert << iteration;
    std::string value = convert.str();

    MPI_File handler;
	MPI_Status status;
	std::string file_name = output +  "/iteration" + value + ".bin";
    MPI_File_open(MPI_COMM_WORLD, file_name.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&handler);
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }

	MPI_File_write_at(handler, offset*sizeof(struct Particle), particle_array, length, mpi_particle, &status);
	MPI_File_close(&handler);
}


void printIterationFile(const Particle* particle_array, const int iteration, const std::string& output, const int rank_id, const int length){

    std::ostringstream convert;
    convert << iteration;
    std::string value = convert.str();

    std::ostringstream rankstring;
    rankstring << rank_id;
    std::string valueRank = rankstring.str();

    FILE* handler;
    std::string file_name = output +  "/rank" + valueRank + "iteration" + value + ".bin";
    handler = fopen(file_name.c_str(), "ab");
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }

    fwrite(particle_array, sizeof(struct Particle), length,handler);

    fclose(handler);
}


void printRankExecutionTimeFile(const double compTime, const std::string& output, const int rank_id){

    std::ostringstream rankstring;
    rankstring << rank_id;
    std::string valueRank = rankstring.str();

    FILE* handler;
    std::string file_name = output +  "/rank_" + valueRank + "_compTime.txt" ;
    handler = fopen(file_name.c_str(), "a");
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }

    fprintf(handler, "%f,",compTime);
    

    fclose(handler);
}


void printExecutionTimeFile(const double compTime, const std::string& output, const int progress){

    FILE* handler;
    std::string file_name = output +  "/exec_compTime.txt" ;
    handler = fopen(file_name.c_str(), "a");
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }

    if(progress==0){
        fprintf(handler, "Halfway execution time: %f\n",compTime);
    }

    if(progress==1){
        fprintf(handler, "Second half execution time: %f\n",compTime);
    }

    if(progress==2){
        fprintf(handler, "Total execution time: %f\n",compTime);
    }
    fclose(handler);
}

/**
 * Compute next particle position
 *
 * This function implements the RK4 algorithm to approximate next particle
 * position. Divergence is checked in this function.
 *
 * @param data: contains the different coil data
 * @param start_point: reference to the particle for which we compute the next step
 * @param step_size: the size of each simulation step (measured in meters)
 * @param mode: mode refers to whether divergence is checked or not. If not, all steps are performed for all particles
 * @return void function. Particle position is updated in place using the particle reference
 */
bool computeIteration(const GlobalData& data, Particle& start_point, const double& step_size, const int mode, Coil* rmi, Coil* rmf,int &divergenceCounter)
{
    bool diverged=false;
    Particle p1 = {0, 0, 0};
    Particle p2 = {0, 0, 0};
    Particle p3 = {0, 0, 0};
    cartesian K1;
    cartesian K2;
    cartesian K3;
    cartesian K4;

    cartesian Ovect = {0, 0, 0};
    Particle p = {0, 0, 0};
    cartesian r_vector;
    double norm_temp;
    double r_radius;


    const double half = 1.0 / 2.0;

    K1 = magnetic_field(rmi, rmf, data, start_point);
    norm_temp = 1.0 / norm_of(K1);
    K1.x = ( K1.x * norm_temp ) * step_size;
    K1.y = ( K1.y * norm_temp ) * step_size;
    K1.z = ( K1.z * norm_temp ) * step_size;
    p1.x = ( K1.x * half ) + start_point.x;
    p1.y = ( K1.y * half ) + start_point.y;
    p1.z = ( K1.z * half ) + start_point.z;

    K2 = magnetic_field(rmi, rmf, data, p1);
    norm_temp = 1.0 / norm_of(K2);
    K2.x = ( K2.x * norm_temp ) * step_size;
    K2.y = ( K2.y * norm_temp ) * step_size;
    K2.z = ( K2.z * norm_temp ) * step_size;
    p2.x = ( K2.x * half ) + start_point.x;
    p2.y = ( K2.y * half ) + start_point.y;
    p2.z = ( K2.z * half ) + start_point.z;

    K3 = magnetic_field(rmi, rmf, data, p2);
    norm_temp = 1.0 / norm_of(K3);
    K3.x = ( K3.x * norm_temp ) * step_size;
    K3.y = ( K3.y * norm_temp ) * step_size;
    K3.z = ( K3.z * norm_temp ) * step_size;
    p3.x = K3.x + start_point.x;
    p3.y = K3.y + start_point.y;
    p3.z = K3.z + start_point.z;

    K4 = magnetic_field(rmi, rmf, data, p3);
    norm_temp = 1.0 / norm_of(K4);
    K4.x = ( K4.x * norm_temp ) * step_size;
    K4.y = ( K4.y * norm_temp ) * step_size;
    K4.z = ( K4.z * norm_temp ) * step_size;
    start_point.x = start_point.x + (( K1.x + 2 * K2.x + 2 * K3.x + K4.x ) / 6 );
    start_point.y = start_point.y + (( K1.y + 2 * K2.y + 2 * K3.y + K4.y ) / 6 );
    start_point.z = start_point.z + (( K1.z + 2 * K2.z + 2 * K3.z + K4.z ) / 6 );

    if (mode == 1)
    {
        p.x = start_point.x;
        p.y = start_point.y;
        Ovect.x = ( p.x / norm_of(p)) * 0.2381; //// Origen vector
        Ovect.y = ( p.y / norm_of(p)) * 0.2381;
        Ovect.z = 0;
        r_vector.x = start_point.x - Ovect.x;
        r_vector.y = start_point.y - Ovect.y;
        r_vector.z = start_point.z - Ovect.z;
        r_radius = norm_of(r_vector);
        if (r_radius > 0.0944165)
        {
	        start_point.x = MINOR_RADIUS;
	        start_point.y = MINOR_RADIUS;
	        start_point.z = MINOR_RADIUS;
            divergenceCounter += 1;
            diverged=true;
        }
     }
    return diverged;
}


void getMagneticProfile(const GlobalData& data, const int num_points, const int phi_angle, const std::string& output, const int dimension){

 	//Prepare parameters for magnetic_field function: rmi, rmf
	Coil rmi[TOTAL_OF_COILS];
	Coil rmf[TOTAL_OF_COILS];
	Particle* observation_points;
	Particle point={0,0,0};
	cartesian B_point;
	const double major_R = 0.2381;
	const double minor_r = 0.0944165;
	double width;
 	double radians = phi_angle*PI/180.0;

    initializeGlobals(rmi, rmf);

    //Prepare output file
    FILE* handler;
    std::string file_name = output + "/magnetic_field.txt";
    std::cout << "Before Handler open" << std::endl;
    handler = fopen(file_name.c_str(), "w");
    std::cout << "After handler open" << std::endl;
    if(nullptr == handler)
    {
        printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
        exit(0);
    }

    fprintf(handler, "x,y,z,|B|\n");

    if(dimension == 1){
        width = (2*minor_r)/num_points;
        observation_points = static_cast<Particle*>(malloc(sizeof(struct Particle) * num_points));
        //Generate observation points at phi_angle plane
        for(int i=0; i<num_points; i++){
            observation_points[i].x = ((major_R-minor_r+(width*i))+minor_r*cos(PI/2))*cos(radians);
            observation_points[i].y = ((major_R-minor_r+(width*i))+minor_r*cos(PI/2))*sin(radians);
            observation_points[i].z = 0.0;
        }
        for(int i=0;i<num_points;i++)
        {
            point.x = observation_points[i].x;
            point.y = observation_points[i].y;
            point.z = observation_points[i].z;
            B_point = magnetic_field(rmi,rmf,data,point);
            fprintf(handler, "%e,%e,%e,%e\n", point.x, point.y, point.z,norm_of(B_point));
        }
    }else if(dimension == 2){
        std::cout << "Entrando a dimension 2" << std::endl;
		width = minor_r/num_points;
        observation_points = static_cast<Particle*>(malloc(sizeof(struct Particle) * (num_points*360)));
        std::cout << "Inicializando puntos de observacion" << std::endl;
		for(int i=0; i<360; i++){
            for(int j=0; j<num_points; j++){
            //std::cout << "it:  " << ((num_points*i)+j) << std::endl;
	    	observation_points[((num_points*i)+j)].x = (major_R+((width*j)*sin(i*(PI/180))))*cos(radians);
            observation_points[((num_points*i)+j)].y = ((width*j)*cos(i*PI/180));
            observation_points[((num_points*i)+j)].z= (major_R+(width*j)* sin(i*PI/180))*sin(radians);
            }
        }
	std::cout << "Inicializacion finalizada" << std::endl;

        for(int i=0;i<num_points*360;i++)
        {
            point.x = observation_points[i].x;
            point.y = observation_points[i].y;
            point.z = observation_points[i].z;
            B_point = magnetic_field(rmi,rmf,data,point);
            fprintf(handler, "%e,%e,%e,%e\n", point.x, point.y, point.z,norm_of(B_point));
        }
	std::cout << "Campo calculado" << std::endl;

    }
     //For each observation point call magnetic_field
	fclose(handler);
	free(observation_points);
	finishGlobal(rmi, rmf);
}

/**
 * Run the actual simulation steps
 *
 * This function coordinates execution of simulation steps.
 * Each MPI_Rank iterate over its particles calling the computeIteration function for
 * each of them.
 *
 * @param data: contains the different coil data
 * @param output: string containing the output directory path ("results_<jobId>")
 * @param particles: Coil structure containing x,y and z arrays of particle positions
 * @param length: length of particle arrays (number of particles to simulate)
 * @param steps: amount of simulation steps to perform
 * @param step_size: the size of each simulation step (measured in meters)
 * @param mode: mode refers to whether divergence is checked or not. If not, all steps are performed for all particles
 * @param print_type: decides if printing is done with tabs or using commas as separation tokens
 * @param startPosition: used to assign names to output files when printing them by particle
 * @return void function. However, an output file per iteration is generated.
 */
void runParticles(const GlobalData& data, const std::string& output, Particle* particles, const int length, const int steps, const double& step_size,
	 const int mode,const int debugFlag)
{
	int myRank, prefixSize, offset;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int number_of_threads = omp_get_max_threads();
    int divergenceCounter=0;
    bool diverged = false;
	int particleIndex;
    double startIOTime = 0;
    double endIOTime = 0;
	double totalIOTime = 0;
    double compStartTime = 0;
    double compHalfTime = 0;
    double compEndTime = 0;
    double rankCompTime = 0;
    double balancedTotalCompTime = 0;
    double secondCompTime = 0; 
    double balancedTime = 0;
    double totalCompTime = 0;

	Coil rmi[TOTAL_OF_COILS];
	Coil rmf[TOTAL_OF_COILS];

	/*********Computing offset for parallel I/O*********************/
	MPI_Scan(&length,&prefixSize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	offset = prefixSize - length;

	/***************************************************************/

    if(debugFlag && (myRank==0)){
        printf("Number of threads being used: %d\n",number_of_threads);
        printf("Rank: %d, prefixSize: %d, offset: %d\n", myRank, prefixSize, offset*sizeof(struct Particle));
        printf("Running timesteps computations\n");
    }


	if(myRank==0 && debugFlag){ startIOTime=MPI_Wtime(); }
	
    //Write starting particle positions to file named: iteration0.bin
    /*printParallelIterationFile(particles, 0, output, myRank, length, offset);
    
    //printIterationFile(particles, 0, output, myRank,length);
	
    /*if(myRank==0 && debugFlag){ 
        endIOTime = MPI_Wtime();
        totalIOTime+=endIOTime-startIOTime;
    }*/

    MPI_Barrier(MPI_COMM_WORLD);
    compStartTime = MPI_Wtime();

    //This loop is where computeIterations is called and whre IO is also called.
	#pragma omp parallel shared(particles) private(rmi,rmf,particleIndex)
	 {

        initializeGlobals(rmi, rmf);

	 	for (int i = 1; i <= steps; i++){
		    #pragma omp for schedule(runtime) reduction(+:divergenceCounter)
	  	    for(int j=0; j < length ; j++)
	  	    {
	  			if((particles[j].x == MINOR_RADIUS) && (particles[j].y == MINOR_RADIUS) && (particles[j].z == MINOR_RADIUS)){
                    continue;
		 		}else{
	                diverged = computeIteration(data,particles[j],step_size,mode,rmi,rmf, divergenceCounter);
				}

                /*if(diverged){
                    particleIndex = (myRank*length)+j;
                    printf("Diverged: %d\n", particleIndex);    
                    diverged = false; 
                }*/
	  	    }

            
            
            /*if (i%(steps/2)==0 && i!=steps){
                #pragma omp single
                {
                    MPI_Barrier(MPI_COMM_WORLD);
                    compHalfTime = MPI_Wtime();
                    if (myRank==0){
                        printExecutionTimeFile(compHalfTime-compStartTime, output, 0);
                    }    
                }        
                    
            }*/    
                
            
		    /*#pragma omp single
		    {
		     	if(i%100 == 0){

					if(myRank==0 && debugFlag){ startIOTime=MPI_Wtime(); }

		     		printParallelIterationFile(particles, i, output, myRank, length, offset);

					if(myRank==0 && debugFlag){ 
                        endIOTime = MPI_Wtime();
                        totalIOTime+=endIOTime-startIOTime;
                    }
				}
	 	    }*/
 
	  	}

		 finishGlobal(rmi,rmf);

	}

    // Each rank computes its total execution time (rankCompTime)
    // Each rank computes its execution time for the second half of the iterations (secondCompTime)
    compEndTime = MPI_Wtime();
    rankCompTime = compEndTime - compStartTime;
    //secondCompTime = compEndTime-compHalfTime;
    
    // Synchronizing all ranks to compute total execution time and balanced iterations time
    
    MPI_Barrier(MPI_COMM_WORLD);
    totalCompTime = MPI_Wtime() - compStartTime;
    if(myRank==0){
        //balancedTotalCompTime = MPI_Wtime();
        //balancedTime = balancedTotalCompTime - compHalfTime;
        //printExecutionTimeFile(balancedTime, output, 1);
        printExecutionTimeFile(totalCompTime, output, 2);
    }
    
    if(debugFlag){
        //printRankExecutionTimeFile(secondCompTime, output, myRank);
        printf("Rank %d, computation time: %f\n",myRank, rankCompTime);
        printf("Rank %d, divergence counter: %d\n",myRank,divergenceCounter);
        int totalDiverged;
        MPI_Reduce(&divergenceCounter,&totalDiverged, 1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);

        if(myRank == 0){
            printf("Number of diverging particles: %d\n", totalDiverged);
        	printf("Total time in IO: %f\n", totalIOTime);
        }
    }
}
