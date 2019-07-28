#include "simulation.h"
#include "compute.h"
#include "output.h"

int main()
{
	simulation::setup();
	simulation::run();
	mean_r_sq::output();
}
