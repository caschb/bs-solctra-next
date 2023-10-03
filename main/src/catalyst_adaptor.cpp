#include "catalyst_adaptor.h"
#include <catalyst.hpp>
#include <spdlog/spdlog.h>

void insitu::initialize()
{
  spdlog::debug("Initializing catalyst");
  conduit_cpp::Node node;
  node["catalyst_load/implementation"] = "stub";
}

void insitu::execute()
{}

void insitu::finalize()
{}