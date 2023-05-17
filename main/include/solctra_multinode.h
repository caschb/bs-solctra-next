#ifndef SOLCTRA_SOLCTRA_H
#define SOLCTRA_SOLCTRA_H

#include "utils.h"
#include <cmath>
#include <sstream>
#include <string>

void printIterationFile(const Particle *particle_array, const int iteration,
                        const std::string &output, const int rank_id,
                        const int length);

void printRankExecutionTimeFile(const double compTime,
                                const std::string &output, const int rank_id);

void printExecutionTimeFile(const double compTime, const std::string &output,
                            const int progress);

void runParticles(Coils &coils, Coils &e_roof, LengthSegments &length_segments,
                  const std::string &output, Particles &particles,
                  const int length, const int steps, const double &step_size,
                  const int mode, const int debugFlag);
#endif
