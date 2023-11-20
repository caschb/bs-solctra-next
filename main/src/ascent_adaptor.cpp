#include <ascent_adaptor.h>
#include <ascent.hpp>
#include <conduit_blueprint.hpp>

#include <spdlog/spdlog.h>

void insitu::initialize()
{
  spdlog::debug("Initializing ascent");
  ascent::Ascent ascent_client;
  ascent_client.open();
  ascent_client.close();
}

void insitu::execute()
{}

void insitu::finalize()
{

}