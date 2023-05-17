#include "solctra_multinode.h"
#include "utils.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <mpi_proto.h>
#include <sstream>
#include <string>

void load_coil_data(double *x, double *y, double *z, const std::string &path) {
  for (int num = 0; num < TOTAL_OF_COILS; num++) {
    std::ostringstream convert;
    convert << num;
    std::string value = convert.str();
    std::string tmp = path + "/Bobina" + value + "m.txt";
    loadCartesianFile(
        &(x[num * TOTAL_OF_GRADES_PADDED]), &(y[num * TOTAL_OF_GRADES_PADDED]),
        &(z[num * TOTAL_OF_GRADES_PADDED]), TOTAL_OF_GRADES + 1, tmp);
  }
}

void printIterationFileTxt(Particles &particles, const int iteration,
                           const int rank, const std::string &output) {
  std::ostringstream convert;
  convert << iteration;
  convert << "_";
  convert << rank;
  std::string value = convert.str();
  FILE *handler;
  std::string file_name = output + "/iteration_" + value + ".txt";
  handler = fopen(file_name.c_str(), "a");
  if (nullptr == handler) {
    std::cout << "Unable to open file=[" << file_name.c_str()
              << "]. Nothing to do\n";
    exit(0);
  }
  fprintf(handler, "x,y,z\n");
  for (auto &particle : particles) {
    fprintf(handler, "%f,%f,%f\n", particle.x, particle.y, particle.z);
  }
  fclose(handler);
}

void printIterationFile(const Particle *particle_array, const int iteration,
                        const std::string &output, const int rank_id,
                        const int length) {

  std::ostringstream convert;
  convert << iteration;
  std::string value = convert.str();

  std::ostringstream rankstring;
  rankstring << rank_id;
  std::string valueRank = rankstring.str();

  FILE *handler;
  std::string file_name =
      output + "/rank" + valueRank + "iteration" + value + ".bin";
  handler = fopen(file_name.c_str(), "ab");
  if (nullptr == handler) {
    printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
    exit(0);
  }

  fwrite(particle_array, sizeof(Particle), length, handler);

  fclose(handler);
}

void printRankExecutionTimeFile(const double compTime,
                                const std::string &output, const int rank_id) {

  std::ostringstream rankstring;
  rankstring << rank_id;
  std::string valueRank = rankstring.str();

  FILE *handler;
  std::string file_name = output + "/rank_" + valueRank + "_compTime.txt";
  handler = fopen(file_name.c_str(), "a");
  if (nullptr == handler) {
    printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
    exit(0);
  }

  fprintf(handler, "%f,", compTime);

  fclose(handler);
}

void printExecutionTimeFile(const double compTime, const std::string &output,
                            const int progress) {

  FILE *handler;
  std::string file_name = output + "/exec_compTime.txt";
  handler = fopen(file_name.c_str(), "a");
  if (nullptr == handler) {
    printf("Unable to open file=[%s]. Nothing to do\n", file_name.c_str());
    exit(0);
  }

  if (progress == 0) {
    fprintf(handler, "Halfway execution time: %f\n", compTime);
  }

  if (progress == 1) {
    fprintf(handler, "Second half execution time: %f\n", compTime);
  }

  if (progress == 2) {
    fprintf(handler, "Total execution time: %f\n", compTime);
  }
  fclose(handler);
}

bool computeIteration(const Coils &coils, const Coils &e_roof,
                      const LengthSegments &length_segments,
                      Particle &start_point, const double &step_size,
                      const int mode, Coils &rmi, Coils &rmf,
                      int &divergenceCounter) {
  bool diverged = false;
  Particle p1;
  Particle p2;
  Particle p3;

  Cartesian k1;
  Cartesian k2;
  Cartesian k3;
  Cartesian k4;

  Cartesian zero_vect;
  Particle p;
  Cartesian r_vector;
  double norm_temp;
  double r_radius;

  constexpr double half = 1.0 / 2.0;
  k1 = computeMagneticField(coils, e_roof, rmi, rmf, length_segments,
                            start_point);
  norm_temp = 1.0 / norm_of(k1);
  k1.x = (k1.x * norm_temp) * step_size;
  k1.y = (k1.y * norm_temp) * step_size;
  k1.z = (k1.z * norm_temp) * step_size;
  p1.x = (k1.x * half) + start_point.x;
  p1.y = (k1.y * half) + start_point.y;
  p1.z = (k1.z * half) + start_point.z;

  k2 = computeMagneticField(coils, e_roof, rmi, rmf, length_segments, p1);
  norm_temp = 1.0 / norm_of(k2);
  k2.x = (k2.x * norm_temp) * step_size;
  k2.y = (k2.y * norm_temp) * step_size;
  k2.z = (k2.z * norm_temp) * step_size;
  p2.x = (k2.x * half) + start_point.x;
  p2.y = (k2.y * half) + start_point.y;
  p2.z = (k2.z * half) + start_point.z;

  k3 = computeMagneticField(coils, e_roof, rmi, rmf, length_segments, p2);
  norm_temp = 1.0 / norm_of(k3);
  k3.x = (k3.x * norm_temp) * step_size;
  k3.y = (k3.y * norm_temp) * step_size;
  k3.z = (k3.z * norm_temp) * step_size;
  p3.x = k3.x + start_point.x;
  p3.y = k3.y + start_point.y;
  p3.z = k3.z + start_point.z;

  k4 = computeMagneticField(coils, e_roof, rmi, rmf, length_segments, p3);
  norm_temp = 1.0 / norm_of(k4);
  k4.x = (k4.x * norm_temp) * step_size;
  k4.y = (k4.y * norm_temp) * step_size;
  k4.z = (k4.z * norm_temp) * step_size;
  start_point.x = start_point.x + ((k1.x + 2 * k2.x + 2 * k3.x + k4.x) / 6);
  start_point.y = start_point.y + ((k1.y + 2 * k2.y + 2 * k3.y + k4.y) / 6);
  start_point.z = start_point.z + ((k1.z + 2 * k2.z + 2 * k3.z + k4.z) / 6);

  if (mode == 1) {
    p.x = start_point.x;
    p.y = start_point.y;
    zero_vect.x = (p.x / norm_of(p)) * 0.2381; //// Origen vector
    zero_vect.y = (p.y / norm_of(p)) * 0.2381;
    zero_vect.z = 0;
    r_vector.x = start_point.x - zero_vect.x;
    r_vector.y = start_point.y - zero_vect.y;
    r_vector.z = start_point.z - zero_vect.z;
    r_radius = norm_of(r_vector);
    if (r_radius > 0.0944165) {
      start_point.x = MINOR_RADIUS;
      start_point.y = MINOR_RADIUS;
      start_point.z = MINOR_RADIUS;
      divergenceCounter += 1;
      diverged = true;
    }
  }
  return diverged;
}

void runParticles(Coils &coils, Coils &e_roof, LengthSegments &length_segments,
                  const std::string &output, Particles &particles,
                  const int length, const int steps, const double &step_size,
                  const int mode, const int debugFlag) {
  int myRank, prefixSize, offset;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  int divergenceCounter = 0;
  double totalIOTime = NAN;
  double compStartTime = NAN;
  double compEndTime = NAN;
  double rankCompTime = NAN;
  double totalCompTime = NAN;

  Coils rmi;
  Coils rmf;

  MPI_Scan(&length, &prefixSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  offset = prefixSize - length;

  if (myRank == 0 && debugFlag) {
    std::cout << "Running timestep computations\n";
    std::cout << "Rank: " << myRank << ", prefixSize: " << prefixSize
              << ", offset: " << offset << '\n';
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printIterationFileTxt(particles, 0, myRank, output);
  compStartTime = MPI_Wtime();

  for (int step = 1; step <= steps; ++step) {
    for (auto &particle : particles) {
      if ((particle.x == MINOR_RADIUS) && (particle.y == MINOR_RADIUS) &&
          (particle.z == MINOR_RADIUS)) {
        continue;
      } else {
        computeIteration(coils, e_roof, length_segments, particle, step_size,
                         mode, rmi, rmf, divergenceCounter);
      }
    }
    if (step % 10 == 0) {
      printIterationFileTxt(particles, step, myRank, output);
    }
  }

  compEndTime = MPI_Wtime();
  rankCompTime = compEndTime - compStartTime;

  MPI_Barrier(MPI_COMM_WORLD);
  totalCompTime = MPI_Wtime() - compStartTime;

  if (myRank == 0) {
    printExecutionTimeFile(totalCompTime, output, 2);
  }

  if (debugFlag) {
    std::cout << "Rank " << myRank << ", computation time: " << rankCompTime
              << '\n';
    std::cout << "Rank " << myRank
              << ", divergence counter: " << divergenceCounter << '\n';
    int totalDiverged;
    MPI_Reduce(&divergenceCounter, &totalDiverged, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (myRank == 0) {
      std::cout << "Number of diverging particles: " << totalDiverged << '\n';
      std::cout << "Total time in IO: " << totalIOTime << '\n';
    }
  }
}