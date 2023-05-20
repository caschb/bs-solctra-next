#include "solctra_multinode.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <ostream>
#include <sstream>
#include <string>
#include <utils.h>
#include <vector>

constexpr auto DEFAULT_STEPS = 100000u;
constexpr auto DEFAULT_STEP_SIZE = 0.001;
constexpr auto DEFAULT_MODE = 1u;
constexpr auto DEFAULT_RESOURCES = std::string("resources");
constexpr auto DEFAULT_MAGPROF = 0u;
constexpr auto DEFAULT_NUM_POINTS = 10000u;
constexpr auto DEFAULT_PHI_ANGLE = 0;
constexpr auto DEFAULT_DIMENSION = 1u;
constexpr auto DEFAULT_DEBUG = 0u;

auto getStepsFromArgs(const int &argc, char **argv) {
  for (int i = 1; i < argc - 1; ++i) {
    std::string tmp = argv[i];
    if (tmp == "-steps") {
      return static_cast<unsigned>(atoi(argv[i + 1]));
    }
  }
  return DEFAULT_STEPS;
}
auto getStepSizeFromArgs(const int &argc, char **argv) {
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

auto getResourcePath(const int &argc, char **argv) {
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
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, nullptr);

  auto my_rank = 0u;
  auto comm_size = 0u;
  auto name_len = 0u;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, (int *)&comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, (int *)&my_rank);
  MPI_Get_processor_name(processor_name, (int *)&name_len);

  /********Create MPI particle type*****************/
  auto MPI_Cartesian = setupMPICartesianType();
  auto MPI_Coil = setupMPIArray(MPI_Cartesian, TOTAL_OF_GRADES);
  auto MPI_LengthSegment = setupMPIArray(MPI_DOUBLE, TOTAL_OF_GRADES);
  /*************************************************/

  /*******Declaring program and runtime parameters*************/
  auto resource_path = DEFAULT_RESOURCES; // Coil directory path
  auto steps = DEFAULT_STEPS;       // Amount of simulation steps
  auto step_size = DEFAULT_STEP_SIZE;    // Size of each simulation step

  /*Variables for magnetic profile diagnostic*/
  auto magprof = DEFAULT_MAGPROF; // Flag to control whether magnetic profile is computed
  auto num_points = DEFAULT_NUM_POINTS; // Number of sampling points for magnetic profile
  auto phi_angle = DEFAULT_PHI_ANGLE; // Angle at which the magnetic profile will be computed
  /******************************************/

  auto length = 0u; // Amount of particles to simulate
  auto debug_flag = DEFAULT_DEBUG;
  auto mode = DEFAULT_MODE; // Check divergence of simulation or not
  auto dimension = DEFAULT_DIMENSION;
  std::string output; // Path of results directory
  std::string jobId;  // JobID in the cluster
  std::ofstream handler;
  /*******Declaring program and runtime parameters*************/

  // Rank 0 reads input parameters from the command line
  // A log file is created to document the runtime parameters
  if (my_rank == 0) {

    resource_path = getResourcePath(argc, argv);
    steps = getStepsFromArgs(argc, argv);
    step_size = getStepSizeFromArgs(argc, argv);
    length = getParticlesLengthFromArgs(argc, argv);
    mode = getModeFromArgs(argc, argv);
    debug_flag = getDebugFromArgs(argc, argv);
    magprof = getMagneticProfileFromArgs(argc, argv);
    num_points = getNumPointsFromArgs(argc, argv);
    phi_angle = getAngleFromArgs(argc, argv);
    jobId = getJobId(argc, argv);
    dimension = getDimension(argc, argv);
    output = "results_" + jobId;
    createDirectoryIfNotExists(output);

    std::cout << "Communicator Size=[" << comm_size << "]." << std::endl;
    std::cout << "Running with:" << std::endl;
    std::cout << "Resource Path=[" << resource_path << "]." << std::endl;
    std::cout << "JobId=[" << jobId << "]." << std::endl;
    std::cout << "Steps=[" << steps << "]." << std::endl;
    std::cout << "Steps size=[" << step_size << "]." << std::endl;
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
    handler << "Steps size=[" << step_size << "]." << std::endl;
    handler << "Particles=[" << length << "]." << std::endl;
    handler << "Mode=[" << mode << "]." << std::endl;
    handler << "Output path=[" << output << "]." << std::endl;
    handler << "MPI size=[" << comm_size << "]." << std::endl;
    handler << "Rank=[" << my_rank << "] => Processor Name=[" << processor_name
            << "]." << std::endl;
  }

  /*********** Rank 0 distributes runtime parameters amongst ranks********/
  MPI_Bcast(&steps, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&step_size, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&length, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mode, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&debug_flag, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  int output_size = output.size();
  MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (0 != my_rank) {
    output.resize(output_size);
  }
  MPI_Bcast(const_cast<char *>(output.data()), output_size, MPI_CHAR, 0, MPI_COMM_WORLD);
  /*********** Rank 0 distributes runtime parameters amongst ranks********/

  /*********** Rank 0 reads in all particles ******/
  std::vector<Particle> particles(length);

  double startInitializationTime = 0.0;
  double endInitializationTime = 0.0;

  // Only rank 0 reads the information from the input file
  if (my_rank == 0) {
    if (debug_flag) {
      startInitializationTime = MPI_Wtime();
    }
    LoadParticles(argc, argv, particles, length, my_rank);

    if (debug_flag) {
      endInitializationTime = MPI_Wtime();
      std::cout << "Total initialization time=["
                << (endInitializationTime - startInitializationTime) << "]."
                << std::endl;
    }
    std::cout << "Particles initialized\n";
  }

  int myShare = length / comm_size;

  if (my_rank < length % comm_size) {
    myShare = myShare + 1;
  }

  std::vector<int> groupMyShare(comm_size);
  std::vector<int> displacements(comm_size);

  MPI_Gather(&myShare, 1, MPI_INT, &groupMyShare.front(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  if (my_rank == 0) {
    for (unsigned int i = 1; i < comm_size; i++) {
      displacements[i] = displacements[i - 1] + groupMyShare[i - 1];
    }
  }

  MPI_Bcast(&displacements.front(), comm_size, MPI_INT, 0, MPI_COMM_WORLD);

  Particles local_particles(myShare);

  MPI_Scatterv(particles.data(), groupMyShare.data(), displacements.data(),
               MPI_Cartesian, local_particles.data(), myShare, MPI_Cartesian, 0,
               MPI_COMM_WORLD);

  Coils coils;
  Coils e_roof;
  LengthSegments length_segments;
  if (my_rank == 0) {
    loadCoilData(coils, resource_path);
    computeERoof(coils, e_roof, length_segments);
  }

  MPI_Bcast(&coils.front(), TOTAL_OF_COILS, MPI_Coil, 0, MPI_COMM_WORLD);
  MPI_Bcast(&e_roof.front(), TOTAL_OF_COILS, MPI_Coil, 0, MPI_COMM_WORLD);
  MPI_Bcast(&length_segments.front(), TOTAL_OF_COILS, MPI_LengthSegment, 0,
            MPI_COMM_WORLD);

  if (my_rank == 0 && magprof != 0) {
    std::cout << "Computing magnetic profiles" << std::endl;
    computeMagneticProfile(coils, e_roof, length_segments, num_points,
                           phi_angle, dimension);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  double startTime = 0;
  double endTime = 0;
  if (my_rank == 0) {
    startTime = MPI_Wtime();
    std::cout << "Executing simulation" << std::endl;
  }
  runParticles(coils, e_roof, length_segments, output, local_particles, length,
               steps, step_size, mode, debug_flag);
  MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0) {
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
    handler << jobId << "," << length << "," << steps << "," << step_size << ","
            << output << "," << (endTime - startTime) << std::endl;
    handler.close();
    time_t now = time(0);
    char *dt = ctime(&now);

    std::cout << "Timestamp: " << dt << std::endl;
  }

  MPI_Finalize();
  return 0;
}
