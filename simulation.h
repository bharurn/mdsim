#ifndef SIMULATION_H
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
		global::temp = 2.483;

		real sumv_sq = 0, v_x, v_y, v_z, a, b, dis = 1.0;
		
		vec lat_init;
		std::vector<vec> lat;

		std::ifstream ip;
		ip.open("input.xyz");
		
		ip >> global::n_prt;
		
		ip >> global::xlen >> global::ylen >> global::zlen >> global::cut >> global::dT >> global::T_f >> global::dump_time >> global::temp;
				
		for(int i=0; i<global::n_prt; i++)
		{
			ip >> lat_init.x >> lat_init.y >> lat_init.z;
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
			
			sumv_sq = sumv_sq + global::prt_list[i].vel.magnitude();
		}
		 
		global::comv.x = global::comv.x/global::n_prt;
		global::comv.y = global::comv.y/global::n_prt; 
		global::comv.z = global::comv.z/global::n_prt;
		sumv_sq = sumv_sq/global::n_prt;
		real fs = pow((3*global::temp)/sumv_sq, 0.5);
		
		std::vector<atom>::iterator i;
		
		for(i = global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			i->vel.x = (i->vel.x - global::comv.x) * fs;
			i->vel.y = (i->vel.y - global::comv.y) * fs;
			i->vel.z = (i->vel.z - global::comv.z) * fs;

			i->posp.x = i->pos.x - (i->vel.x*global::dT);
			i->posp.y = i->pos.y - (i->vel.y*global::dT);
			i->posp.z = i->pos.z - (i->vel.z*global::dT);
		}

		compute::initList();
		
		mean_r_sq::init();
		
		output::setDumpFile("out.xyz", "log.txt");
	}

	void run()
	{
		for(int time = 0; time <= global::T_f; ++time)
		{
			compute::work();
			mean_r_sq::calculate(time);
			if(time%global::dump_time == 0)
			{	
				output::trajectory();
				output::log(time*global::dT);
				std::cout << time << '\n';
			}	
		}
	}
}

#endif
