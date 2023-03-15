
#include "solctra_multinode.h"
#include <fstream>
#include <mpi.h>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <math.h>
#include <ctime>


const unsigned DEFAULT_STRING_BUFFER = 100;
const unsigned DEFAULT_STEPS = 100000; //100000
const double DEFAULT_STEP_SIZE = 0.001;
const unsigned DEFAULT_PRECISION = 5;
const unsigned DEFAULT_PARTICLES= 1;
const unsigned DEFAULT_MODE= 1;
const std::string DEFAULT_OUTPUT = "results";
const std::string DEFAULT_RESOURCES = "resources";
const unsigned DEFAULT_MAGPROF = 0;
const unsigned DEFAULT_NUM_POINTS = 10000;
const unsigned DEFAULT_PHI_ANGLE = 0;
const unsigned DEFAULT_DIMENSION = 1;
const unsigned DEFAULT_DEBUG = 0;


MPI_Datatype mpi_particle;

unsigned getPrintPrecisionFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-precision")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_PRECISION;
}
unsigned getStepsFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-steps")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_STEPS;
}
double getStepSizeFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-stepSize")
        {
            return strtod(argv[i+1], nullptr);
        }
    }
    return DEFAULT_STEP_SIZE;
}
void LoadParticles(const int& argc, char** argv, Particle* particles, const int length, const int seedValue)
{
    bool found = false;
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-particles")
        {
            loadParticleFile(particles, length, argv[i+1]);
            found = true;
            break;
        }
    }
    if(!found)
    {
        printf("No file given. Initializing random particles\n");
        initializeParticles(particles, length,seedValue);
    }
}

std::string getResourcePath(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string param = argv[i];
        if("-resource" == param)
        {
            return std::string(argv[i+1]);
        }
    }
    return DEFAULT_RESOURCES;
}

unsigned getParticlesLengthFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-length")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    printf("ERROR: You must specify number of particles to simulate\n");
    exit(1);
}

unsigned getDebugFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-d")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_DEBUG;
}

unsigned getModeFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-mode")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_MODE;
}

std::string getJobId(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
        std::string tmp = argv[i];
        if(tmp == "-id")
        {
            return std::string(argv[i+1]);
        }
    }
    printf("ERROR: job id must be given!!\n");
    exit(1);
}

unsigned getMagneticProfileFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
            return static_cast<unsigned>(atoi(argv[i+1]));
        }
    }
    return DEFAULT_MAGPROF;
}


unsigned getNumPointsFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
            return static_cast<unsigned>(atoi(argv[i+2]));
        }
    }
    return DEFAULT_NUM_POINTS;
}

unsigned getAngleFromArgs(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
       	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
     	    return static_cast<unsigned>(atoi(argv[i+3]));
        }
    }
    return DEFAULT_PHI_ANGLE;
}

unsigned getDimension(const int& argc, char** argv)
{
    for(int i = 1 ; i < argc - 1 ; ++i)
    {
     	std::string tmp = argv[i];
        if(tmp == "-magnetic_prof")
        {
            return static_cast<unsigned>(atoi(argv[i+4]));
        }
    }
    return DEFAULT_DIMENSION;
}


int main(int argc, char** argv)
{
    /*****MPI variable declarations and initializations**********/
    int provided;
    MPI_Init_thread(&argc, &argv,MPI_THREAD_FUNNELED,&provided);
    int myRank;
    int commSize;
    int nameLen;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Get_processor_name(processorName, &nameLen);

    /********Create MPI particle type*****************/
    setupParticleType();
    /*************************************************/

    /*******Declaring program and runtime parameters*************/
    std::string resourcePath; //Coil directory path
    unsigned steps; //Amount of simulation steps
    double stepSize; //Size of each simulation step

    /*Variables for magnetic profile diagnostic*/
    unsigned magprof; //Flag to control whether magnetic profile is computed
    unsigned num_points; //Number of sampling points for magnetic profile
    unsigned phi_angle; //Angle at which the magnetic profile will be computed
    /******************************************/

    unsigned precision; //TBD
    unsigned int length; //Amount of particles to simulate
    unsigned int debugFlag;
    unsigned int mode; //Check divergence of simulation or not
    unsigned int dimension;
    std::string output; //Path of results directory
    std::string jobId; //JobID in the cluster
    std::ofstream handler;
    /*******Declaring program and runtime parameters*************/



    //Rank 0 reads input parameters from the command line
    //A log file is created to document the runtime parameters
    if(myRank == 0)
    {

	    resourcePath = getResourcePath(argc, argv);
    	steps = getStepsFromArgs(argc, argv);
    	stepSize = getStepSizeFromArgs(argc, argv);
    	precision = getPrintPrecisionFromArgs(argc, argv);
    	length = getParticlesLengthFromArgs(argc, argv);
    	mode = getModeFromArgs(argc, argv);
        debugFlag = getDebugFromArgs(argc, argv);
    	magprof = getMagneticProfileFromArgs(argc, argv);
    	num_points = getNumPointsFromArgs(argc, argv);
    	phi_angle = getAngleFromArgs(argc, argv);
    	jobId = getJobId(argc, argv);
    	dimension = getDimension(argc,argv);
    	output = "results_" + jobId;
    	createDirectoryIfNotExists(output);

	    std::cout.precision(precision);
    	std::cout << "Communicator Size=[" << commSize << "]." << std::endl;
	    std::cout << "Running with:" << std::endl;
    	std::cout << "Resource Path=[" << resourcePath << "]." << std::endl;
    	std::cout << "JobId=[" << jobId << "]." << std::endl;
    	std::cout << "Steps=[" << steps << "]." << std::endl;
    	std::cout << "Steps size=[" << stepSize << "]." << std::endl;
    	std::cout << "Particles=[" << length << "]." << std::endl;
    	std::cout << "Input Current=[" << I << "] A." << std::endl;
    	std::cout << "Mode=[" << mode << "]." << std::endl;
    	std::cout << "Output path=[" << output << "]." << std::endl;
    	std::string file_name = "stdout_"+jobId+".log";

    	handler.open(file_name.c_str());
    	if(!handler.is_open()){
        	std::cerr << "Unable to open stdout.log for appending. Nothing to do." << std::endl;
        	exit(0);
    	}

    	handler << "Running with:" << std::endl;
    	handler << "Steps=[" << steps << "]." << std::endl;
    	handler << "Steps size=[" << stepSize << "]." << std::endl;
    	handler << "Particles=[" << length << "]." << std::endl;
    	handler << "Mode=[" << mode << "]." << std::endl;
    	handler << "Output path=[" << output << "]." << std::endl;
	    handler << "MPI size=[" << commSize << "]." << std::endl;
	    handler << "Rank=[" << myRank << "] => Processor Name=[" << processorName << "]." << std::endl;
    }

    /*if(debugFlag){
        printf("Rank=[%d] in host=[%s]\n", myRank, processorName);
    }*/

    /*********** Rank 0 distributes runtime parameters amongst ranks********/
    MPI_Bcast(&steps, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stepSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&precision, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mode, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&debugFlag, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    unsigned int outputSize = static_cast<unsigned int>(output.size());
    MPI_Bcast(&outputSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    char* tmp = new char[outputSize + 1];
    if(0 == myRank)
    {
        std::strcpy(tmp, output.c_str());
    }
    MPI_Bcast(tmp, outputSize + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if(0 != myRank)
    {
       output = std::string(tmp);
    }
    delete[] tmp;
    /*********** Rank 0 distributes runtime parameters amongst ranks********/


    /*********** Rank 0 reads in all particles ******/
    Particle* particles = static_cast<Particle*>(_mmm_malloc(sizeof(struct Particle)*length, ALIGNMENT_SIZE));

    double startInitializationTime = 0;
    double endInitializationTime = 0;

    //Only rank 0 reads the information from the input file
    if(myRank==0)
    {
    	if(debugFlag){
            startInitializationTime = MPI_Wtime();
        }
        LoadParticles(argc, argv, particles, length,myRank);
        if(debugFlag){
            endInitializationTime = MPI_Wtime();
            std::cout << "Total initialization time=[" << (endInitializationTime - startInitializationTime) << "]." << std::endl;
        }
    	printf("Particles initialized\n");
    }
    /*********** All ranks declare a memory space for their particle ******/


    /*********** Calculating distribution parameters for scatterv of particles********/
    int myShare = floor(length/commSize); //Each rank computes how many particles it will take care of
    int startPosition; //This parameter was used to print out particle file names.

    /*In case the number of particles is not fully divisible between the number of ranks
      then the excess is distributed among ranks */
    if(myRank < length%commSize)
    {
	    myShare = myShare + 1;

    }

    int *groupmyShare = new int[commSize] ();
    int *displacements = new int[commSize] ();

    /*Rank 0 gathers share amounts from all ranks to then compute the displacement to
      to be used when distributing particles*/
    MPI_Gather(&myShare, 1, MPI_INT, groupmyShare, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(myRank == 0)
    {
	    /*if(debugFlag){
            std::cout << "Rank 0: Group Sizes are: " ;
            for(int i=0; i < commSize; i++)
            {
		        std::cout << groupmyShare[i] << " ";

	        } std::cout << std::endl;
	    }*/
        for(int i=1; i < commSize; i++)
        {
	        displacements[i] = displacements[i-1]+groupmyShare[i-1];
		    /*if(debugFlag){
                std::cout << "Rank " << i << " displacement: " << displacements[i] << std::endl;
            }*/
        }

    }

    MPI_Bcast(displacements,commSize, MPI_INT, 0, MPI_COMM_WORLD);
    startPosition = displacements[myRank];

    /*if(debugFlag){
        std::cout << "Rank " << myRank << " startPosition: " << startPosition << std::endl;
    }*/

    /*This is the rank-local particle array memory declaration and distribution*/
    Particle* localParticles = static_cast<Particle*>(_mmm_malloc(sizeof(struct Particle)*myShare, ALIGNMENT_SIZE));
    MPI_Scatterv(particles, groupmyShare, displacements, mpi_particle, localParticles, myShare, mpi_particle, 0, MPI_COMM_WORLD);

   /*********** Calculating distribution parameters for scatterv of particles********/


    /*********** All ranks declare a memory space for Coil data ******/
    const size_t sizeToAllocate = sizeof(double) * TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS;
    GlobalData data;
    data.coils.x = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.y = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.z = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));

    //Only Rank 0 reads the information from the coil files
    if(myRank == 0)
    {
    	load_coil_data(data.coils.x, data.coils.y, data.coils.z, resourcePath);
    	std::cout << "Coil data loaded" << std::endl;
    }
    //Rank 0 distributes coil data among ranks
    MPI_Bcast(data.coils.x, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data.coils.y, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data.coils.z, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*********** All ranks declare a memory space for Coil data ******/




    /*********** All ranks declare a memory space for e_roof computations ******/
    data.e_roof.x = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.y = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.z = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.leng_segment = static_cast<double*>(_mmm_malloc(sizeToAllocate, ALIGNMENT_SIZE));

    if(myRank == 0){
    	e_roof(data);
    	//std::cout << "e_roof data loaded" << std::endl;
    }
    MPI_Bcast(data.e_roof.x, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data.e_roof.y, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data.e_roof.z, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data.leng_segment, TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*********** All ranks declare a memory space for e_roof computations ******/


    /*********If the magneticProfile is requested, only rank 0 computes it.*****/
    if(myRank == 0 && magprof != 0)
    {
	std::cout << "Computing magnetic profiles" << std::endl;
	getMagneticProfile(data, num_points, phi_angle, output, dimension);
    }
    /******************End magnetic profile computation*************************/



    /******************Call the trajectory simulation********************/
    MPI_Barrier(MPI_COMM_WORLD);
    
    double startTime = 0;
    double endTime = 0;
    if(myRank==0)
    {
    	startTime= MPI_Wtime();
        std::cout << "Executing simulation" << std::endl;
    }
    runParticles(data, output, localParticles, myShare, steps, stepSize, mode, debugFlag);

    MPI_Barrier(MPI_COMM_WORLD);

    /******************End of the trajectory simulation********************/


    _mmm_free(data.coils.x);
    _mmm_free(data.coils.y);
    _mmm_free(data.coils.z);
    _mmm_free(data.e_roof.x);
    _mmm_free(data.e_roof.y);
    _mmm_free(data.e_roof.z);
    _mmm_free(data.leng_segment);
    _mmm_free(particles);
    _mmm_free(localParticles);




    if(myRank == 0){
        endTime= MPI_Wtime();
        std::cout << "Simulation finished" << std::endl;
        std::cout << "Total execution time=[" << (endTime - startTime) << "]." << std::endl;
        handler << "Total execution time=[" << (endTime - startTime) << "]." << std::endl;
    	handler.close();
    	handler.open("stats.csv", std::ofstream::out | std::ofstream::app);
    	if(!handler.is_open())
    	{
        	std::cerr << "Unable to open stats.csv for appending. Nothing to do." << std::endl;
        	exit(0);
    	}
    	handler << jobId << "," << length << "," << steps << "," <<  stepSize << "," << output << "," << (endTime - startTime) << std::endl;
    	handler.close();
    	time_t now = time(0);
    	char* dt = ctime(&now);

    	std::cout << "Timestamp: " << dt << std::endl;
    }


    MPI_Finalize();
    return 0;
}
