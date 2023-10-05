#include "catalyst_adaptor.h"
#include <catalyst.hpp>
#include <catalyst_api.h>
#include <catalyst_conduit.hpp>
#include <spdlog/spdlog.h>

void insitu::initialize()
{
  spdlog::debug("Initializing catalyst");
  conduit_cpp::Node node;
  node["catalyst_load/implementation"] = "stub";
  catalyst_status error = catalyst_initialize(conduit_cpp::c_node(&node));
  if (error != catalyst_status_ok)
  {
    spdlog::error("ERROR: Failed to initialize Catalyst {}", error);
  }
}

void insitu::execute()
{}

void insitu::finalize()
{}