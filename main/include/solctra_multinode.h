#ifndef SOLCTRA_SOLCTRA_H
#define SOLCTRA_SOLCTRA_H

#include "utils.h"
#include <string_view>

void runParticles(Coils &coils, Coils &e_roof, LengthSegments &length_segments,
                  const std::string_view output, Particles &particles,
                  const int length, const int steps, const double &step_size,
                  const int mode, const int debugFlag, Timing& timing);
#endif
