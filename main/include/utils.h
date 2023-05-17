//
// utils.h
// Created by Diego Jimenez
// Basic simulation structures and function prototypes declarations
//

#ifndef SOLCTRA_UTILS_H
#define SOLCTRA_UTILS_H

#define PI 3.141592654
#define MIU 1.2566e-06
#define I -4350
#define MINOR_RADIUS 0.0944165
#define MAJOR_RADIUS 0.2381

#define TOTAL_OF_GRADES 360
#define TOTAL_OF_COILS 12

#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

struct Cartesian {
  double x, y, z;
  void print() { std::cout << x << ',' << y << ',' << z << '\n'; }
  Cartesian(double x = 0.0, double y = 0.0, double z = 0.0)
      : x(x), y(y), z(z) {}
};

typedef std::array<Cartesian, TOTAL_OF_GRADES + 1> Coil;
typedef std::array<Coil, TOTAL_OF_COILS> Coils;
typedef std::array<std::array<double, TOTAL_OF_GRADES>, TOTAL_OF_COILS>
    LengthSegments;
typedef Cartesian Particle;
typedef std::vector<Particle> Particles;

MPI_Datatype setupMPICartesianType();
MPI_Datatype setupMPIArray(MPI_Datatype base_type, int count);

void initializeParticles(Particles &particles, const int seedValue);

void loadParticleFile(Particles &particles, const int numberOfParticles,
                      const std::string_view path);

void loadCoilData(Coils &coil, const std::string_view path);

void computeERoof(Coils &coils, Coils &e_roof, LengthSegments &length_segments);

void computeMagneticProfile(Coils &coils, Coils &e_roof,
                            LengthSegments &length_segments,
                            const int num_points, const int phi_angle,
                            /*const std::string &output,*/ const int dimension);

Cartesian computeMagneticField(const Coils &coils, const Coils &e_roof,
                               Coils &rmi, Coils &rmf,
                               const LengthSegments &length_segments,
                               const Particle &point);

double getCurrentTime();
void createDirectoryIfNotExists(const std::string &path);
bool directoryExists(const std::string &path);
std::string getZeroPadded(const int num);
double randomGenerator(const double min, const double max, const int seedValue);
inline double norm_of(const Cartesian &vec) {
  return sqrt((vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z));
}

#endif
