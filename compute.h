#ifndef COMPUTE_H
#define COMPUTE_H

#include "global.h"

namespace compute
{
	std::vector< std::vector<int> > neigh_list;
	
	real fric;
	const real Q = 500; 
	
	void initList()
	{
		std::vector<int> newRow;
		
		for(int i = 0; i < global::n_prt; i++)
		{
			neigh_list.push_back(newRow);
		}
				
	}
	
	void refreshList()
	{
		std::vector<atom>::iterator i,j;
		vec dr;
		
		for(i= global::prt_list.begin(); i != global::prt_list.end()-1; i++)
		{
			neigh_list[i - global::prt_list.begin()].clear();
			for(j= i+1; j != global::prt_list.end(); j++)
			{
				dr.x = i->pos.x - j->pos.x;
				if (dr.x >= 0.5 * global::input["xlen"]) dr.x -= global::input["xlen"];
				else if (dr.x < -0.5 * global::input["xlen"]) dr.x += global::input["xlen"];
				
				dr.y = i->pos.y - j->pos.y;
				if (dr.y >= 0.5 * global::input["ylen"]) dr.y -= global::input["ylen"];
				else if (dr.y < -0.5 * global::input["ylen"]) dr.y += global::input["ylen"];

				dr.z = i->pos.z - j->pos.z;
				if (dr.z >= 0.5 * global::input["zlen"]) dr.z -= global::input["zlen"];
				else if (dr.z < -0.5 * global::input["zlen"]) dr.z += global::input["zlen"];
				
				if(abs(dr.magnitude()) < (global::input["cut"]+global::input["skin"])*(global::input["cut"]+global::input["skin"]))
					neigh_list[i - global::prt_list.begin()].push_back(j - global::prt_list.begin());
			}
		}
	}
	
	void forces()
	{
		std::vector<atom>::iterator i,k;
		std::vector<int>::iterator l;
		real f_by_dr, r_sq, r_i_sq, r_i_six;
		vec dr;
		global::pe = 0;

		for(k= global::prt_list.begin(); k != global::prt_list.end(); k++)
		{
			k->force.x = 0;
			k->force.y = 0;
			k->force.z = 0;
		}

		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			for(l= neigh_list[i - global::prt_list.begin()].begin(); l != neigh_list[i - global::prt_list.begin()].end(); l++)
			{				
				dr.x = i->pos.x - global::prt_list[*l].pos.x;
				if (dr.x >= 0.5 * global::input["xlen"]) dr.x -= global::input["xlen"];
				else if (dr.x < -0.5 * global::input["xlen"]) dr.x += global::input["xlen"];
				
				dr.y = i->pos.y - global::prt_list[*l].pos.y;
				if (dr.y >= 0.5 * global::input["ylen"]) dr.y -= global::input["ylen"];
				else if (dr.y < -0.5 * global::input["ylen"]) dr.y += global::input["ylen"];

				dr.z = i->pos.z - global::prt_list[*l].pos.z;
				if (dr.z >= 0.5 * global::input["zlen"]) dr.z -= global::input["zlen"];
				else if (dr.z < -0.5 * global::input["zlen"]) dr.z += global::input["zlen"];

				r_sq = dr.magnitude();
				
				if(r_sq < (global::input["cut"]*global::input["cut"]))
				{
					r_i_sq = 1/r_sq;
					r_i_six =  pow(r_i_sq, 3);
					
					f_by_dr = 48 * r_i_sq * r_i_six * ( r_i_six - 0.5); 
					
					i->force.x += f_by_dr * dr.x;
					global::prt_list[*l].force.x -= f_by_dr *dr.x;

					i->force.y += f_by_dr *dr.y;
					global::prt_list[*l].force.y -= f_by_dr *dr.y;

					i->force.z += f_by_dr * dr.z;
					global::prt_list[*l].force.z -= f_by_dr * dr.z;

					global::pe += 4 * r_i_six * (r_i_six - 1) - global::ecut;
				}
			}
		}
	}

	void pbc()
	{
		std::vector<atom>::iterator i;
		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			if (i->pos.x >= 0.5 * global::input["xlen"])
			{
				i->pos.x -= global::input["xlen"];
				i->ix++;
			} 
			else if (i->pos.x < -0.5 * global::input["xlen"])
			{
				i->pos.x += global::input["xlen"];
				i->ix--;
			}

			else if (i->pos.y >= 0.5 * global::input["ylen"])
			{
				i->pos.y -= global::input["ylen"];
				i->iy++;
			}
			else if (i->pos.y < -0.5 * global::input["ylen"])
			{
				i->pos.y += global::input["ylen"];
				i->iy--;
			}						
			else if (i->pos.z >= 0.5 * global::input["zlen"])
			{
				i->pos.z -= global::input["zlen"];
				i->iz++;
			}
			else if (i->pos.z < -0.5 * global::input["zlen"])
			{
				i->pos.z += global::input["zlen"];
				i->iz--;
			}
		}
	}

	void verlet_nve()
	{
		real dT_sq = pow(global::input["dT"], 2);
		global::sumv_sq = 0;
		vec v;
		std::vector<atom>::iterator i,j;
		global::comv.x = 0;
		global::comv.y = 0;
		global::comv.z = 0;
		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			v.x = (2*i->pos.x) - (i->posp.x) + (i->force.x*dT_sq);
			i->vel.x = (v.x - (i->posp.x))/(2*global::input["dT"]);
			global::comv.x += i->vel.x;
			i->posp.x = i->pos.x;
			i->pos.x = v.x;
			i->dist.x += i->pos.x - i->posp.x;

			v.y = (2*i->pos.y) - (i->posp.y) + (i->force.y*dT_sq);
			i->vel.y = (v.y - (i->posp.y))/(2*global::input["dT"]);
			global::comv.y += i->vel.y;
			i->posp.y = i->pos.y;
			i->pos.y = v.y;
			i->dist.y += i->pos.y - i->posp.y;
			
			v.z = (2*i->pos.z) - (i->posp.z) + (i->force.z*dT_sq);
			i->vel.z = (v.z - (i->posp.z))/(2*global::input["dT"]);
			global::comv.z += i->vel.z;
			i->posp.z = i->pos.z;
			i->pos.z = v.z;
			i->dist.z += i->pos.z - i->posp.z;

			global::sumv_sq += i->vel.magnitude();
		}
		
	}
	
	void verlet_nvt()
	{
		std::vector<atom>::iterator i;
		global::comv.x = 0;
		global::comv.y = 0;
		global::comv.z = 0;
		
		fric += (global::sumv_sq - (3*global::n_prt + 1)*global::set_temp)*global::input["dT"]*0.25/Q; 
		
		global::sumv_sq = 0;
		
		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			i->v.x = i->vel.x + .5*global::input["dT"]*(i->force.x - fric*i->vel.x);
			i->posp.x = i->pos.x;
			i->pos.x += i->v.x*global::input["dT"];
			i->dist.x += i->pos.x - i->posp.x;
			
			i->v.y = i->vel.y + .5*global::input["dT"]*(i->force.y - fric*i->vel.y);
			i->posp.y = i->pos.y;
			i->pos.y += i->v.y*global::input["dT"];
			i->dist.y += i->pos.y - i->posp.y;
			
			i->v.z = i->vel.z + .5*global::input["dT"]*(i->force.z - fric*i->vel.z);
			i->posp.z = i->pos.z;
			i->pos.z += i->v.z*global::input["dT"];
			i->dist.z += i->pos.z - i->posp.z;
			
			global::sumv_sq += i->v.magnitude();
		
		}
		
		fric += (global::sumv_sq - (3*global::n_prt + 1)*global::set_temp)*global::input["dT"]*0.25/Q; 
		
		compute::forces();
		
		global::sumv_sq = 0;
		
		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			i->vel.x = (i->v.x + 0.5*global::input["dT"]*i->force.x)/(1 + (0.5*global::input["dT"]*fric));
			global::comv.x += i->vel.x;
			
			i->vel.y = (i->v.y + 0.5*global::input["dT"]*i->force.y)/(1 + (0.5*global::input["dT"]*fric));
			global::comv.y += i->vel.y;
			
			i->vel.z = (i->v.z + 0.5*global::input["dT"]*i->force.z)/(1 + (0.5*global::input["dT"]*fric));
			global::comv.z += i->vel.z;
			
			global::sumv_sq += i->vel.magnitude(); 
		}
		
	}

	void work()
	{
		std::vector<atom>::iterator i;
		real dist_max=0;
		
		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			if(dist_max < i->dist.magnitude())
			{
				dist_max = i->dist.magnitude();
			}
		}
		
		if(dist_max > global::input["skin"]*global::input["skin"])
		{
			compute::refreshList();
			for(i= global::prt_list.begin(); i != global::prt_list.end(); i++) i->dist = 0;
		}
		compute::forces();
		if(global::input["ens"] == 0)	
			compute::verlet_nve();
		else if (global::input["ens"] == 1)
			compute::verlet_nvt();
		
		global::comv.x = global::comv.x/global::n_prt;
		global::comv.y = global::comv.y/global::n_prt; 
		global::comv.z = global::comv.z/global::n_prt;
		
		global::input["temp"] = global::sumv_sq/(3*global::n_prt + 1); 
		global::etot = (global::pe + 0.5*global::sumv_sq)/global::n_prt;
		
		compute::pbc();
	}
}
#endif
