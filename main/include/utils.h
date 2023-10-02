//
// utils.h
// Created by Diego Jimenez
// Basic simulation structures and function prototypes declarations
//

#ifndef SOLCTRA_UTILS_H
#define SOLCTRA_UTILS_H

#include <array>
#include <cmath>
#include <mpi.h>
#include <spdlog/spdlog.h>
#include <string>
#include <string_view>
#include <vector>

struct Cartesian
{
  double x{ 0.0 }, y{ 0.0 }, z{ 0.0 };
  void print() { spdlog::info("{}, {}, {}", x, y, z); }
  Cartesian(double x_e = 0.0, double y_e = 0.0, double z_e = 0.0) : x(x_e), y(y_e), z(z_e) {}
  friend std::ostream &operator<<(std::ostream &os, const Cartesian &cartesian)
  {
    os << cartesian.x << ',' << cartesian.y << ',' << cartesian.z;
    return os;
  }
};

constexpr auto PI = std::numbers::pi;
constexpr auto MIU = 1.2566e-06;
constexpr auto I = -4350;
constexpr auto MINOR_RADIUS = 0.0944165;
constexpr auto MAJOR_RADIUS = 0.2381;
constexpr auto TOTAL_OF_GRADES = 360;
constexpr auto TOTAL_OF_COILS = 12;

using Coil = std::array<Cartesian, TOTAL_OF_GRADES + 1>;
using Coils = std::array<Coil, TOTAL_OF_COILS>;
using LengthSegments = std::array<std::array<double, TOTAL_OF_GRADES>, TOTAL_OF_COILS>;
using Particle = Cartesian;
using Particles = std::vector<Particle>;
using Timing = Cartesian;
using Timings = std::vector<Timing>;

MPI_Datatype setupMPICartesianType();
MPI_Datatype setupMPIArray(MPI_Datatype base_type, int count);

void initializeParticles(Particles &particles, const int seedValue);

void loadParticleFile(Particles &particles, const int numberOfParticles, const std::string_view path);

void loadCoilData(Coils &coil, const std::string_view path);

void computeERoof(Coils &coils, Coils &e_roof, LengthSegments &length_segments);

void computeMagneticProfile(Coils &coils,
  Coils &e_roof,
  LengthSegments &length_segments,
  const int num_points,
  const int phi_angle,
  /*const std::string &output,*/ const int dimension);

Cartesian computeMagneticField(const Coils &coils,
  const Coils &e_roof,
  Coils &rmi,
  Coils &rmf,
  const LengthSegments &length_segments,
  const Particle &point);

auto getCurrentTime();
void createDirectoryIfNotExists(const std::string &path);
auto directoryExists(const std::string &path);
auto getZeroPadded(const int num);
double randomGenerator(const double min, const double max, const int seedValue);
inline auto norm_of(const Cartesian &vec) { return std::sqrt((vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z)); }

#endif
