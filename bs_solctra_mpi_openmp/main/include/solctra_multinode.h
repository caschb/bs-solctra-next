#ifndef SOLCTRA_SOLCTRA_H
#define SOLCTRA_SOLCTRA_H

#include <string>
#include <cmath>
#include "utils.h"
#include <sstream>
void load_coil_data(double* x, double* y, double* z, const std::string& path);
void e_roof(GlobalData& data);
void printeroof(GlobalData& data, const int subsetIndex, const std::string& outputPath);
cartesian magnetic_field(Coil* rmi, Coil* rmf, const GlobalData& data, const Particle& point);
void initializeGlobals(Coil* rmi, Coil* rmf);
void finishGlobals(Coil* rmi, Coil* rmf);
void getMagneticProfile(const GlobalData& data, const int num_points, const int phi_angle, const std::string& output, const int dimension);
bool computeIteration(const GlobalData& data, Particle& start_point, const double& step_size, const int mode, Coil* rmi, Coil* rmf,int &divergenceCounter);
inline double norm_of(const cartesian& vec)
{
    return sqrt(( vec.x * vec.x ) + ( vec.y * vec.y ) + ( vec.z * vec.z ));
}
inline double norm_of(const Particle& vec)
{
    return sqrt(( vec.x * vec.x ) + ( vec.y * vec.y ) + ( vec.z * vec.z ));
}
void printIterationFileTxt(const Particle *particle_array, const int iteration, const int rank, const std::string &output, const int length);
void printParallelIterationFile(const Particle* particle_array, const int iteration, const std::string& output, const int rank_id, const int length, const int offset);
void printIterationFile(const Particle* particle_array, const int iteration, const std::string& output, const int rank_id, const int length);
void printRankExecutionTimeFile(const double compTime, const std::string& output, const int rank_id);
void printExecutionTimeFile(const double compTime, const std::string& output, const int progress);
void runParticles(const GlobalData& data, const std::string& output, Particle* particles, const int length, const int steps, const double& step_size, const int mode, const int debugFlag);

#endif //SOLCTRA_SOLCTRA_H
