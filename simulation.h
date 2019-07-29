##ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include "global.h"
#include "output.h"
#include "compute.h"
#include "mean_r_sq.h"

namespace simulation
{
	void setup()
	{
		std::cout << "Start simulation setup...\n";
		global::input["xlen"] = 100;
		global::input["ylen"] = 100;
		global::input["zlen"] = 100;
		global::input["dT"] = 0.003;
		global::input["temp"] = 1;
		global::input["T_f"] = 1000;
		global::input["ens"] = 0;
		global::input["mean_r_sq_dump"] = 0;
		global::input["cut"] = 2.5;
		global::input["dump_time"] = 1;
		global::input["skin"] = 0.3;
		global::input["mean_r_sq_time_length"] = 2;
		
		real v_x, v_y, v_z, a, b, dis = 1.0;
		global::sumv_sq = 0;
		
		vec lat_init;
		std::vector<vec> lat;

		std::ifstream ip;
		ip.open("input.xyz");
		
		char c;
		
		ip >> global::n_prt;
		
		std::ifstream in;
		in.open("input.in");
		
		std::string param;
		real param_val;
		
		while(!in.eof())
		{
			in >> param >> param_val;
			global::input[param] = param_val;
		}
		
		global::T_f = (int)global::input["T_f"];
		global::dump_time = (int)global::input["dump_time"];
		
		global::ecut = 4 * pow(global::input["cut"], -6) * (pow(global::input["cut"], -6) - 1);
				
		if(global::input["ens"] == 1) global::set_temp = global::input["temp"];
		for(int i=0; i<global::n_prt; i++)
		{
			ip >> c >> lat_init.x >> lat_init.y >> lat_init.z;
			lat.push_back(lat_init);
		}

		for(int i = 0; i < global::n_prt; i++)
		{
			srand (i);
			v_x = ((rand() % 10 + 0)-5)/10.0;

			srand (i+1);
			v_y = ((rand() % 10 + 0)-5)/10.0;
			
			srand (i+2);
			v_z = ((rand() % 10 + 0)-5)/10.0;

			global::prt_list.push_back(atom( lat[i], vec(v_x, v_y, v_z) ) );
			
			global::comv.x += global::prt_list[i].vel.x;
			global::comv.y += global::prt_list[i].vel.y;
			global::comv.z += global::prt_list[i].vel.z;
			
			global::sumv_sq = global::sumv_sq + global::prt_list[i].vel.magnitude();
		}
		 
		global::comv.x = global::comv.x/global::n_prt;
		global::comv.y = global::comv.y/global::n_prt; 
		global::comv.z = global::comv.z/global::n_prt;
		
		real fs = pow((3*global::input["temp"]*global::n_prt)/global::sumv_sq, 0.5);
		
		std::vector<atom>::iterator i;
		
		for(i = global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			i->vel.x = (i->vel.x - global::comv.x) * fs;
			i->vel.y = (i->vel.y - global::comv.y) * fs;
			i->vel.z = (i->vel.z - global::comv.z) * fs;

			i->posp.x = i->pos.x - (i->vel.x*global::input["dT"]);
			i->posp.y = i->pos.y - (i->vel.y*global::input["dT"]);
			i->posp.z = i->pos.z - (i->vel.z*global::input["dT"]);
		}

		compute::initList();
		
		compute::forces();
		
		if(global::input["mean_r_sq_dump"]) mean_r_sq::init();
		
		output::setDumpFile("out.xyz", "log.txt");
	}

	void run()
	{
		std::cout << "Running simulation...\nTime Step\n";
		for(int time = 0; time <= global::T_f; ++time)
		{
			compute::work();
			if(global::input["mean_r_sq_dump"])
				mean_r_sq::calculate(time);
			if(time%global::dump_time == 0)
			{	
				output::trajectory();
				output::log(time*global::input["dT"]);
				std::cout << time << '\n';
			}	
		}
	}
}

#endif
