#include "solctra_multinode.h"
#include <CLI/App.hpp>
#include <CLI/CLI.hpp>
#include <CLI/Validators.hpp>
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

using std::string;

constexpr auto DEFAULT_STEPS = 100000U;
constexpr auto DEFAULT_STEP_SIZE = 0.001;
constexpr auto DEFAULT_MODE = 1U;
constexpr auto DEFAULT_RESOURCES = std::string("resources");
constexpr auto DEFAULT_MAGPROF = 0U;
constexpr auto DEFAULT_NUM_POINTS = 10000U;
constexpr auto DEFAULT_PHI_ANGLE = 0;
constexpr auto DEFAULT_DIMENSION = 1U;
constexpr auto DEFAULT_DEBUG = 0U;

void load_particles(const std::string &particles_file, Particles &particles, const int length, const int seed_value)
{
  std::cout << particles_file << '\n';
  if (!particles_file.empty()) {
    loadParticleFile(particles, length, particles_file);
  } else {
    initializeParticles(particles, seed_value);
  }
}

int main(int argc, char **argv)
{
  CLI::App app("BS-Solctra");

  /*****MPI variable declarations and initializations**********/
  auto provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  auto my_rank = 0U;
  auto comm_size = 0U;
  auto name_len = 0U;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, reinterpret_cast<int *>(&comm_size));
  MPI_Comm_rank(MPI_COMM_WORLD, reinterpret_cast<int *>(&my_rank));
  MPI_Get_processor_name(processor_name, reinterpret_cast<int *>(&name_len));

  /********Create MPI particle type*****************/
  auto *MPI_Cartesian = setupMPICartesianType();
  auto *MPI_Coil = setupMPIArray(MPI_Cartesian, TOTAL_OF_GRADES);
  auto *MPI_LengthSegment = setupMPIArray(MPI_DOUBLE, TOTAL_OF_GRADES);
  /*************************************************/

  /*******Declaring program and runtime parameters*************/
  auto resource_path = DEFAULT_RESOURCES;// Coil directory path
  auto steps = DEFAULT_STEPS;// Amount of simulation steps
  auto step_size = DEFAULT_STEP_SIZE;// Size of each simulation step

  /*Variables for magnetic profile diagnostic*/
  auto magprof = DEFAULT_MAGPROF;// Flag to control whether magnetic profile is computed
  auto num_points = DEFAULT_NUM_POINTS;// Number of sampling points for magnetic profile
  auto phi_angle = DEFAULT_PHI_ANGLE;// Angle at which the magnetic profile will be computed
  /******************************************/

  auto length{ 0U };// Amount of particles to simulate
  auto debug_flag{ DEFAULT_DEBUG };
  auto mode{ DEFAULT_MODE };// Check divergence of simulation or not
  auto dimension{ DEFAULT_DIMENSION };
  std::string particles_file;
  std::string output;// Path of results directory
  std::string jobId;// JobID in the cluster
  std::ofstream handler;
  /*******Declaring program and runtime parameters*************/

  // Rank 0 reads input parameters from the command line
  // A log file is created to document the runtime parameters
  if (my_rank == 0) {
    app.add_option("--particles", particles_file, "Particles file")->check(CLI::ExistingFile);
    app.add_option("--resource", resource_path, "Resource path")->check(CLI::ExistingPath)->capture_default_str();
    app.add_option("--steps", steps, "Number of simulation steps")->check(CLI::PositiveNumber);
    app.add_option("--step-size", step_size, "Step Size")->check(CLI::PositiveNumber);
    app.add_option("--length", length, "Length of the particle array")->check(CLI::PositiveNumber);
    app.add_option("--mode", mode, "Mode of the simulation")->check(CLI::PositiveNumber);
    app.add_option("--debug", debug_flag, "Enable debug messages")->check(CLI::PositiveNumber);
    app.add_option("--magnetic-profile", magprof, "Magnetic profile");
    app.add_option("--num-points", num_points, "Number of points");
    app.add_option("--phi-angle", phi_angle, "Angle phi");
    app.add_option("--dimension", "Dimension");
    app.add_option("--job-id", jobId, "Job ID")->required()->check(CLI::PositiveNumber);
    CLI11_PARSE(app, argc, argv);
    output = "results_" + jobId;
    createDirectoryIfNotExists(output);
    std::string file_name = "stdout_" + jobId + ".log";

    handler.open(file_name.c_str());
    if (!handler.is_open()) {
      std::cerr << "Unable to open stdout.log for appending. Nothing to do." << std::endl;
      exit(0);
    }

    handler << "Running with:" << std::endl;
    handler << "Steps=[" << steps << "]." << std::endl;
    handler << "Steps size=[" << step_size << "]." << std::endl;
    handler << "Particles=[" << length << "]." << std::endl;
    handler << "Mode=[" << mode << "]." << std::endl;
    handler << "Output path=[" << output << "]." << std::endl;
    handler << "MPI size=[" << comm_size << "]." << std::endl;
    handler << "Rank=[" << my_rank << "] => Processor Name=[" << processor_name << "]." << std::endl;
  }

  /*********** Rank 0 distributes runtime parameters amongst ranks********/
  MPI_Bcast(&steps, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&step_size, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&length, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mode, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(&debug_flag, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  int output_size = output.size();
  MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (0 != my_rank) { output.resize(output_size); }
  MPI_Bcast(const_cast<char *>(output.data()), output_size, MPI_CHAR, 0, MPI_COMM_WORLD);
  /*********** Rank 0 distributes runtime parameters amongst ranks********/

  /*********** Rank 0 reads in all particles ******/
  std::vector<Particle> particles(length);

  double startInitializationTime = 0.0;
  double endInitializationTime = 0.0;

  // Only rank 0 reads the information from the input file
  if (my_rank == 0) {
    if (debug_flag != 0U) { startInitializationTime = MPI_Wtime(); }
    load_particles(particles_file, particles, length, my_rank);

    if (debug_flag != 0U) {
      endInitializationTime = MPI_Wtime();
      std::cout << "Total initialization time=[" << (endInitializationTime - startInitializationTime) << "]."
                << std::endl;
    }
    std::cout << "Particles initialized\n";
  }

  int myShare = length / comm_size;

  if (my_rank < length % comm_size) { myShare = myShare + 1; }

  std::vector<int> groupMyShare(comm_size);
  std::vector<int> displacements(comm_size);

  MPI_Gather(&myShare, 1, MPI_INT, &groupMyShare.front(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    for (unsigned int i = 1; i < comm_size; i++) { displacements[i] = displacements[i - 1] + groupMyShare[i - 1]; }
    std::cout << "Sending displacements\n";
  }

  MPI_Bcast(&displacements.front(), comm_size, MPI_INT, 0, MPI_COMM_WORLD);

  Particles local_particles(myShare);

  if (my_rank == 0) { std::cout << "Sending particles\n"; }

  MPI_Scatterv(particles.data(),
    groupMyShare.data(),
    displacements.data(),
    MPI_Cartesian,
    local_particles.data(),
    myShare,
    MPI_Cartesian,
    0,
    MPI_COMM_WORLD);

  Coils coils;
  Coils e_roof;
  LengthSegments length_segments;
  if (my_rank == 0) {
    std::cout << "Loading coils\n";
    loadCoilData(coils, resource_path);
    computeERoof(coils, e_roof, length_segments);
    std::cout << "Broadcasting coil data\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&coils.front(), TOTAL_OF_COILS, MPI_Coil, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    std::cout << "Broadcasting e_roof\n";
  }

  MPI_Bcast(&e_roof.front(), TOTAL_OF_COILS, MPI_Coil, 0, MPI_COMM_WORLD);
  if (my_rank == 0) {
    std::cout << "Broadcasting length segments\n";
  }
  MPI_Bcast(&length_segments.front(), TOTAL_OF_COILS, MPI_LengthSegment, 0, MPI_COMM_WORLD);

  if (my_rank == 0 && magprof != 0) {
    std::cout << "Computing magnetic profiles" << std::endl;
    computeMagneticProfile(coils, e_roof, length_segments, num_points, phi_angle, dimension);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  double startTime = 0;
  double endTime = 0;
  if (my_rank == 0) {
    startTime = MPI_Wtime();
    std::cout << "Executing simulation" << std::endl;
  }
  Timings timings(comm_size);
  runParticles(coils,
    e_roof,
    length_segments,
    output,
    local_particles,
    length,
    steps,
    step_size,
    mode,
    debug_flag,
    timings[my_rank]);

  MPI_Gatherv(local_particles.data(),
    myShare,
    MPI_Cartesian,
    particles.data(),
    groupMyShare.data(),
    displacements.data(),
    MPI_Cartesian,
    0,
    MPI_COMM_WORLD);
  MPI_Gather(&timings[my_rank], 1, MPI_Cartesian, timings.data(), 1, MPI_Cartesian, 0, MPI_COMM_WORLD);

  if (my_rank == 0) {
    for (const auto &timing : timings) {
      std::cout << timing << '\n';
      handler << timing << '\n';
    }
    endTime = MPI_Wtime();
    std::cout << "Simulation finished" << std::endl;
    std::cout << "Total execution time=[" << (endTime - startTime) << "]." << std::endl;
    handler << "Total execution time=[" << (endTime - startTime) << "]." << std::endl;
    handler.close();
    handler.open("stats.csv", std::ofstream::out | std::ofstream::app);
    if (!handler.is_open()) {
      std::cerr << "Unable to open stats.csv for appending. Nothing to do." << std::endl;
      std::quick_exit(0);
    }
    handler << jobId << "," << length << "," << steps << "," << step_size << "," << output << ","
            << (endTime - startTime) << std::endl;
    handler.close();
    time_t now = time(nullptr);
    char *deltaTime = ctime(&now);

    std::cout << "Timestamp: " << deltaTime << std::endl;
  }

  MPI_Finalize();
  return 0;
}
