#include "simulation.h"
#include "compute.h"
#include "output.h"
#include <iostream>
#include <conio.h>

int main()
{
	std::cout << "*****MDSim -- Molecular Dynamics Simulation by by Bharath Raghavan*****\n";
	simulation::setup();
	simulation::run();
	if(global::input["mean_r_sq_dump"]) mean_r_sq::output();
	std::cout << "Simulation Done.";
}
