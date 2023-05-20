#include "solctra_multinode.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <utils.h>
#include <vector>

constexpr auto DEFAULT_STEPS = 100000;
constexpr auto DEFAULT_STEP_SIZE = 0.001;
constexpr auto DEFAULT_PRECISION = 5;
constexpr auto DEFAULT_MODE = 1;
constexpr auto DEFAULT_RESOURCES = "resources";
constexpr auto DEFAULT_MAGPROF = 0;
constexpr auto DEFAULT_NUM_POINTS = 10000;
constexpr auto DEFAULT_PHI_ANGLE = 0;
constexpr auto DEFAULT_DIMENSION = 1;
constexpr auto DEFAULT_DEBUG = 0;

unsigned getPrintPrecisionFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-precision") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  return DEFAULT_PRECISION;
}
unsigned getStepsFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-steps") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  return DEFAULT_STEPS;
}
double getStepSizeFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-stepSize") {
      return strtod(argv[i + 1], nullptr);
    }
  }
  return DEFAULT_STEP_SIZE;
}

void LoadParticles(const int &argc, char **argv, Particles &particles,
                   const int length, const int seedValue) {
  bool found = false;
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-particles") {
      std::cout << argv[i + 1] << '\n';
      loadParticleFile(particles, length, argv[i + 1]);
      found = true;
      break;
    }
  }
  if (!found) {
    std::cout << "No file given. Initializing random particles\n";
    initializeParticles(particles, seedValue);
  }
}

std::string getResourcePath(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string param = argv[i];
    if ("-resource" == param) {
      return std::string(argv[i + 1]);
    }
  }
  return DEFAULT_RESOURCES;
}

unsigned getParticlesLengthFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-length") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  std::cerr << "ERROR: You must specify number of particles to simulate\n";
  exit(1);
}

unsigned getDebugFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-d") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  return DEFAULT_DEBUG;
}

unsigned getModeFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-mode") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  return DEFAULT_MODE;
}

std::string getJobId(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-id") {
      return std::string(argv[i + 1]);
    }
  }
  std::cerr << "ERROR: job id must be given!!\n";
  exit(1);
}

unsigned getMagneticProfileFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-magnetic_prof") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  return DEFAULT_MAGPROF;
}

unsigned getNumPointsFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-magnetic_prof") {
      return static_cast<unsigned>(atoi(argv[i + 2]));
    }
  }
  return DEFAULT_NUM_POINTS;
}

unsigned getAngleFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-magnetic_prof") {
      return static_cast<unsigned>(atoi(argv[i + 3]));
    }
  }
  return DEFAULT_PHI_ANGLE;
}

unsigned getDimension(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-magnetic_prof") {
      return static_cast<unsigned>(atoi(argv[i + 4]));
    }
  }
  return DEFAULT_DIMENSION;
}

int main(int argc, char **argv) {
  /*****MPI variable declarations and initializations**********/
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

  unsigned int myRank = 0;
  unsigned int commSize = 0;
  int nameLen = 0;
  char processorName[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, (int *)&commSize);
  MPI_Comm_rank(MPI_COMM_WORLD, (int *)&myRank);
  MPI_Get_processor_name(processorName, &nameLen);

  /********Create MPI particle type*****************/
  auto MPI_Cartesian = setupMPICartesianType();
  auto MPI_Coil = setupMPIArray(MPI_Cartesian, TOTAL_OF_GRADES);
  auto MPI_LengthSegment = setupMPIArray(MPI_DOUBLE, TOTAL_OF_GRADES);
  /*************************************************/

  /*******Declaring program and runtime parameters*************/
  std::string resourcePath; // Coil directory path
  unsigned steps = 0;       // Amount of simulation steps
  double stepSize = NAN;    // Size of each simulation step

  /*Variables for magnetic profile diagnostic*/
  unsigned magprof = 0; // Flag to control whether magnetic profile is computed
  unsigned num_points = 0; // Number of sampling points for magnetic profile
  unsigned phi_angle =
      0; // Angle at which the magnetic profile will be computed
  /******************************************/

  unsigned int length; // Amount of particles to simulate
  unsigned precision;  // TBD
  unsigned int debugFlag;
  unsigned int mode; // Check divergence of simulation or not
  unsigned int dimension;
  std::string output; // Path of results directory
  std::string jobId;  // JobID in the cluster
  std::ofstream handler;
  /*******Declaring program and runtime parameters*************/

  // Rank 0 reads input parameters from the command line
  // A log file is created to document the runtime parameters
  if (myRank == 0) {

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
    dimension = getDimension(argc, argv);
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
    std::string file_name = "stdout_" + jobId + ".log";

    handler.open(file_name.c_str());
    if (!handler.is_open()) {
      std::cerr << "Unable to open stdout.log for appending. Nothing to do."
                << std::endl;
      exit(0);
    }

    handler << "Running with:" << std::endl;
    handler << "Steps=[" << steps << "]." << std::endl;
    handler << "Steps size=[" << stepSize << "]." << std::endl;
    handler << "Particles=[" << length << "]." << std::endl;
    handler << "Mode=[" << mode << "]." << std::endl;
    handler << "Output path=[" << output << "]." << std::endl;
    handler << "MPI size=[" << commSize << "]." << std::endl;
    handler << "Rank=[" << myRank << "] => Processor Name=[" << processorName
            << "]." << std::endl;
  }

  /*********** Rank 0 distributes runtime parameters amongst ranks********/
  MPI_Bcast(&steps, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&stepSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&precision, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&length, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mode, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&debugFlag, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  unsigned int outputSize = static_cast<unsigned int>(output.size());
  MPI_Bcast(&outputSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  char *tmp = new char[outputSize + 1];
  if (0 == myRank) {
    std::strcpy(tmp, output.c_str());
  }
  MPI_Bcast(tmp, outputSize + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (0 != myRank) {
    output = std::string(tmp);
  }
  delete[] tmp;
  /*********** Rank 0 distributes runtime parameters amongst ranks********/

  /*********** Rank 0 reads in all particles ******/
  std::vector<Particle> particles(length);

  double startInitializationTime = 0.0;
  double endInitializationTime = 0.0;

  // Only rank 0 reads the information from the input file
  if (myRank == 0) {
    if (debugFlag) {
      startInitializationTime = MPI_Wtime();
    }
    LoadParticles(argc, argv, particles, length, myRank);

    if (debugFlag) {
      endInitializationTime = MPI_Wtime();
      std::cout << "Total initialization time=["
                << (endInitializationTime - startInitializationTime) << "]."
                << std::endl;
    }
    std::cout << "Particles initialized\n";
  }

  int myShare = length / commSize;

  if (myRank < length % commSize) {
    myShare = myShare + 1;
  }

  std::vector<int> groupMyShare(commSize);
  std::vector<int> displacements(commSize);

  MPI_Gather(&myShare, 1, MPI_INT, &groupMyShare.front(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  if (myRank == 0) {
    for (unsigned int i = 1; i < commSize; i++) {
      displacements[i] = displacements[i - 1] + groupMyShare[i - 1];
    }
  }

  MPI_Bcast(&displacements.front(), commSize, MPI_INT, 0, MPI_COMM_WORLD);

  Particles local_particles(myShare);

  MPI_Scatterv(particles.data(), groupMyShare.data(), displacements.data(),
               MPI_Cartesian, local_particles.data(), myShare, MPI_Cartesian, 0,
               MPI_COMM_WORLD);

  Coils coils;
  Coils e_roof;
  LengthSegments length_segments;
  if (myRank == 0) {
    loadCoilData(coils, resourcePath);
    computeERoof(coils, e_roof, length_segments);
  }

  MPI_Bcast(&coils.front(), TOTAL_OF_COILS, MPI_Coil, 0, MPI_COMM_WORLD);
  MPI_Bcast(&e_roof.front(), TOTAL_OF_COILS, MPI_Coil, 0, MPI_COMM_WORLD);
  MPI_Bcast(&length_segments.front(), TOTAL_OF_COILS, MPI_LengthSegment, 0,
            MPI_COMM_WORLD);

  if (myRank == 0 && magprof != 0) {
    std::cout << "Computing magnetic profiles" << std::endl;
    computeMagneticProfile(coils, e_roof, length_segments, num_points,
                           phi_angle, dimension);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  double startTime = 0;
  double endTime = 0;
  if (myRank == 0) {
    startTime = MPI_Wtime();
    std::cout << "Executing simulation" << std::endl;
  }
  runParticles(coils, e_roof, length_segments, output, local_particles, length,
               steps, stepSize, mode, debugFlag);
  MPI_Barrier(MPI_COMM_WORLD);

  if (myRank == 0) {
    endTime = MPI_Wtime();
    std::cout << "Simulation finished" << std::endl;
    std::cout << "Total execution time=[" << (endTime - startTime) << "]."
              << std::endl;
    handler << "Total execution time=[" << (endTime - startTime) << "]."
            << std::endl;
    handler.close();
    handler.open("stats.csv", std::ofstream::out | std::ofstream::app);
    if (!handler.is_open()) {
      std::cerr << "Unable to open stats.csv for appending. Nothing to do."
                << std::endl;
      exit(0);
    }
    handler << jobId << "," << length << "," << steps << "," << stepSize << ","
            << output << "," << (endTime - startTime) << std::endl;
    handler.close();
    time_t now = time(0);
    char *dt = ctime(&now);

    std::cout << "Timestamp: " << dt << std::endl;
  }

  MPI_Finalize();
  return 0;
}
