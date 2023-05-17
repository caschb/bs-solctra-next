//
//
// Created by lchavarr on 4/19/16.
//

#include "utils.h"
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <mpi.h>
#include <random>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>

void loadCartesianFile(double *x, double *y, double *z, const int length,
                       const std::string &path) {
  FILE *file_buff;
  // Open file
  file_buff = fopen(path.c_str(), "r");
  if (file_buff == nullptr) {
    printf("Error al abrir archivo \n");
    exit(0);
  } else {
    double localX, localY, localZ;
    printf("Loading %s with length=%d\n", path.c_str(), length);
    for (int point = 0; point < length; point++) {
      fscanf(file_buff, "%le %le %le", &localX, &localY, &localZ);
      x[point] = localX;
      y[point] = localY;
      z[point] = localZ;
    }
    fclose(file_buff);
  }
}

void loadParticleFile(Particles &particles, const int numberOfParticles,
                      const std::string_view path) {
  static const auto delimeter = "\t";
  std::ifstream particles_file(path.data());
  std::string line;
  auto line_number = 0;
  while (std::getline(particles_file, line) &&
         line_number < numberOfParticles) {
    size_t position = 0;
    std::array<double, 3> data;
    auto idx = 0;
    while (position != std::string::npos && position < line.size()) {
      auto current_pos = position;
      position = line.find(delimeter, position);
      auto tok = line.substr(current_pos, position - current_pos);
      data[idx] = strtod(tok.c_str(), nullptr);
      idx += 1;
      if (position == std::string::npos || position >= line.size()) {
        break;
      }
      position += 1;
    }
    particles[line_number] = Particle(data[0], data[1], data[2]);
    line_number += 1;
  }
  particles_file.close();
}

void loadCoilData(Coils &coils, const std::string_view path) {

  static const auto delimeter = "\t";

  for (int coil_number = 0; coil_number < TOTAL_OF_COILS; coil_number++) {
    auto coil_number_token = std::to_string(coil_number);
    std::ostringstream filename_oss;
    filename_oss << path << "/Bobina" << coil_number_token << "m.txt";
    std::string filename = filename_oss.str();
    std::ifstream coil_file(filename);
    std::string line;
    auto line_number = 0;
    while (std::getline(coil_file, line)) {
      size_t position = 0;
      std::array<double, 3> data;
      auto idx = 0;
      while (position != std::string::npos && position < line.size()) {
        auto current_pos = position;
        position = line.find(delimeter, position);
        auto tok = line.substr(current_pos, position - current_pos);
        data[idx] = strtod(tok.c_str(), nullptr);
        idx += 1;
        if (position == std::string::npos || position >= line.size()) {
          break;
        }
        position += 1;
      }
      coils[coil_number][line_number] = Cartesian(data[0], data[1], data[2]);
      line_number += 1;
    }
    coil_file.close();
  }
}

Cartesian computeMagneticField(const Coils &coils, const Coils &e_roof,
                               Coils &rmi, Coils &rmf,
                               const LengthSegments &length_segments,
                               const Particle &point) {
  static const auto multiplier = (miu * I) / (4.0 * PI);
  Cartesian B;

  for (auto i = 0; i < TOTAL_OF_COILS; ++i) {
    for (auto j = 0; j < TOTAL_OF_GRADES; ++j) {
      rmi[i][j].x = point.x - coils[i][j].x;
      rmi[i][j].y = point.y - coils[i][j].y;
      rmi[i][j].z = point.z - coils[i][j].z;
      rmf[i][j].x = point.x - coils[i][j + 1].x;
      rmf[i][j].y = point.y - coils[i][j + 1].y;
      rmf[i][j].z = point.z - coils[i][j + 1].z;

      const auto norm_rmi = norm_of(rmi[i][j]);
      const auto norm_rmf = norm_of(rmf[i][j]);

      Cartesian U;
      U.x = multiplier * e_roof[i][j].x;
      U.y = multiplier * e_roof[i][j].y;
      U.z = multiplier * e_roof[i][j].z;

      const auto C = (((2 * (length_segments[i][j]) * (norm_rmi + norm_rmf)) /
                       (norm_rmi * norm_rmf)) *
                      ((1) / ((norm_rmi + norm_rmf) * (norm_rmi + norm_rmf) -
                              length_segments[i][j] * length_segments[i][j])));

      Cartesian V;
      V.x = C * rmi[i][j].x;
      V.y = C * rmi[i][j].y;
      V.z = C * rmi[i][j].z;

      B.x = B.x + ((U.y * V.z) - (U.z * V.y));
      B.y = B.y - ((U.x * V.z) - (U.z * V.x));
      B.z = B.z + ((U.x * V.y) - (U.y * V.x));
    }
  }
  return B;
}

void computeMagneticProfile(
    Coils &coils, Coils &e_roof, LengthSegments &length_segments,
    const int num_points, const int phi_angle,
    /*const std::string &output,*/ const int dimension) {
  Coils rmi;
  Coils rmf;
  Cartesian point;
  Cartesian b_point;
  static const auto major_R = 0.2381;
  static const auto minor_r = 0.0944165;
  auto width = 0.0;
  auto radians = phi_angle * PI / 180.0;

  // TODO: prepare output file

  if (dimension == 1) {
    width = (2 * minor_r) / num_points;
    std::vector<Particle> observation_particles(num_points);
    for (auto i = 0; auto &particle : observation_particles) {
      particle.x = ((major_R - minor_r + (width * i)) + minor_r * cos(PI / 2)) *
                   cos(radians);
      particle.y = ((major_R - minor_r + (width * i)) + minor_r * cos(PI / 2)) *
                   sin(radians);
      particle.z = 0.0;
      i += 1;
    }
    for (auto &particle : observation_particles) {
      point.x = particle.x;
      point.y = particle.y;
      point.z = particle.z;
      b_point =
          computeMagneticField(coils, e_roof, rmi, rmf, length_segments, point);
      // TODO: write points
    }
  } else if (dimension == 2) {
    width = minor_r / num_points;
    std::vector<Particle> observation_particles(num_points * TOTAL_OF_GRADES);
    for (auto count = 0; auto &particle : observation_particles) {
      auto i = count / TOTAL_OF_GRADES;
      auto j = count % num_points;
      count += 1;
      particle.x =
          (major_R + ((width * j) * sin(i * (PI / 180)))) * cos(radians);
      particle.y = ((width * j) * cos(i * PI / 180));
      particle.z = (major_R + (width * j) * sin(i * PI / 180)) * sin(radians);
    }

    for (auto &particle : observation_particles) {
      point.x = particle.x;
      point.y = particle.y;
      point.z = particle.z;
      b_point =
          computeMagneticField(coils, e_roof, rmi, rmf, length_segments, point);
      // TODO: write points
    }
  }
}

void computeERoof(Coils &coils, Coils &e_roof,
                  LengthSegments &length_segments) {
  Cartesian segment;
  for (int i = 0; i < TOTAL_OF_COILS; i++) {
    // #pragma GCC ivdep
    for (int j = 0; j < TOTAL_OF_GRADES; j++) {

      segment.x = (coils[i][j + 1].x) - (coils[i][j].x);
      segment.y = (coils[i][j + 1].y) - (coils[i][j].y);
      segment.z = (coils[i][j + 1].z) - (coils[i][j].z);

      length_segments[i][j] = norm_of(segment);

      const double length_segment_inverted = 1.0 / length_segments[i][j];
      e_roof[i][j].x = segment.x * length_segment_inverted;
      e_roof[i][j].y = segment.y * length_segment_inverted;
      e_roof[i][j].z = segment.z * length_segment_inverted;
    }
  }
}

double randomGenerator(const double min, const double max,
                       const int seedValue) {
  // std::random_device rd; //Used to obtain a seed for the random number
  std::mt19937 gen(seedValue);
  std::uniform_real_distribution<double> distribution(min, max);
  double result = distribution(gen);
  return result;
}

/*double randomGeneratorExponential(const double maxNumber){
    std::random_device rd; //Used to obtain a seed for the random number
    std::mt19937 gen(rd());
    std::exponential_distribution<double> distribution(1.0);
    double result = distribution(gen);
    if(result>1.0){
        double factor = gen.max()/maxNumber;
        result = result/factor;
    }else{ result = result*maxNumber;}
    return result;
}*/

void initializeParticles(Particles &particles, const int seedValue) {
  double radius;
  double minorRadius;
  double toroidalAngle;
  double poloidalAngle;
  double toroidalRadians;
  double poloidalRadians;
  auto length = particles.size();
  printf("Initializing with seed %d\n", seedValue);
  for (size_t point = 0; point < length; point++) {
    radius = randomGenerator(0.0, 0.4, seedValue);
    toroidalAngle = randomGenerator(0.0, 360.0, seedValue);
    poloidalAngle = randomGenerator(0.0, 360.0, seedValue);

    toroidalRadians = toroidalAngle * PI / 180.0;
    poloidalRadians = poloidalAngle * PI / 180.0;
    minorRadius = MINOR_RADIUS * radius;
    particles[point].x = (MAJOR_RADIUS + (minorRadius)*cos(poloidalRadians)) *
                         cos(toroidalRadians);
    particles[point].y = (MAJOR_RADIUS + (minorRadius)*cos(poloidalRadians)) *
                         sin(toroidalRadians);
    particles[point].z = (minorRadius)*sin(poloidalRadians);
  }
}

double getCurrentTime() {
  struct timeval tod;
  gettimeofday(&tod, nullptr);
  return static_cast<double>(tod.tv_sec) +
         static_cast<double>(tod.tv_usec) * 1.0e-6;
}

void createDirectoryIfNotExists(const std::string &path) {
  if (!directoryExists(path)) {
    const int error =
        mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (0 > error) {
      printf("Error creating directory!n");
      exit(1);
    }
  }
}

bool directoryExists(const std::string &path) {
  struct stat info;
  if (stat(path.c_str(), &info) != 0) {
    return false;
  } else if (info.st_mode & S_IFDIR) {
    return true;
  } else {
    return false;
  }
}

std::string getZeroPadded(const int num) {
  std::ostringstream convert;
  convert << num;
  std::string value = convert.str();
  const size_t numSize = value.size();
  int numberZeros = 5; // Number of 0s wanted in file name:
                       // particle00000...X.txt
  for (size_t i = 0; i < numberZeros - numSize; ++i) {
    value = "0" + value;
  }
  return value;
}

MPI_Datatype setupMPICartesianType() {
  MPI_Datatype MPI_Cartesian;
  static const int NDIMS = 3;
  int blockLengths[NDIMS] = {1, 1, 1};
  MPI_Datatype types[NDIMS] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets[NDIMS];

  offsets[0] = offsetof(Cartesian, x);
  offsets[1] = offsetof(Cartesian, y);
  offsets[2] = offsetof(Cartesian, z);

  MPI_Type_create_struct(NDIMS, blockLengths, offsets, types, &MPI_Cartesian);
  MPI_Type_commit(&MPI_Cartesian);
  return MPI_Cartesian;
}

MPI_Datatype setupMPIArray(MPI_Datatype base_type, int count) {
  MPI_Datatype MPI_Arr;
  MPI_Type_contiguous(count, base_type, &MPI_Arr);
  MPI_Type_commit(&MPI_Arr);
  return MPI_Arr;
}
