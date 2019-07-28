#ifndef COMPUTE_H
#define COMPUTE_H

#include "global.h"

namespace compute
{
	std::vector< std::vector<int> > neigh_list; 
	
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
				if (dr.x >= 0.5 * global::xlen) dr.x -= global::xlen;
				else if (dr.x < -0.5 * global::xlen) dr.x += global::xlen;
				
				dr.y = i->pos.y - j->pos.y;
				if (dr.y >= 0.5 * global::ylen) dr.y -= global::ylen;
				else if (dr.y < -0.5 * global::ylen) dr.y += global::ylen;

				dr.z = i->pos.z - j->pos.z;
				if (dr.z >= 0.5 * global::zlen) dr.z -= global::zlen;
				else if (dr.z < -0.5 * global::zlen) dr.z += global::zlen;
				
				if(abs(dr.magnitude()) < (global::cut+global::skin)*(global::cut+global::skin))
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
			for(l= neigh_list[i - global::prt_list.begin()].begin();
			 l != neigh_list[i - global::prt_list.begin()].end(); l++)
			{				
				dr.x = i->pos.x - global::prt_list[*l].pos.x;
				if (dr.x >= 0.5 * global::xlen) dr.x -= global::xlen;
				else if (dr.x < -0.5 * global::xlen) dr.x += global::xlen;
				
				dr.y = i->pos.y - global::prt_list[*l].pos.y;
				if (dr.y >= 0.5 * global::ylen) dr.y -= global::ylen;
				else if (dr.y < -0.5 * global::ylen) dr.y += global::ylen;

				dr.z = i->pos.z - global::prt_list[*l].pos.z;
				if (dr.z >= 0.5 * global::zlen) dr.z -= global::zlen;
				else if (dr.z < -0.5 * global::zlen) dr.z += global::zlen;

				r_sq = dr.magnitude();
				
				if(r_sq < (global::cut*global::cut))
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
			if (i->pos.x >= 0.5 * global::xlen)
			{
				i->pos.x -= global::xlen;
				i->ix++;
			} 
			else if (i->pos.x < -0.5 * global::xlen)
			{
				i->pos.x += global::xlen;
				i->ix--;
			}

			else if (i->pos.y >= 0.5 * global::ylen)
			{
				i->pos.y -= global::ylen;
				i->iy++;
			}
			else if (i->pos.y < -0.5 * global::ylen)
			{
				i->pos.y += global::ylen;
				i->iy--;
			}						
			else if (i->pos.z >= 0.5 * global::zlen)
			{
				i->pos.z -= global::zlen;
				i->iz++;
			}
			else if (i->pos.z < -0.5 * global::zlen)
			{
				i->pos.z += global::zlen;
				i->iz--;
			}
		}
	}

	void verlet()
	{
		real sumv_sq=0, dT_sq = pow(global::dT, 2);
		vec v;
		std::vector<atom>::iterator i,j;
		global::comv.x = 0;
		global::comv.y = 0;
		global::comv.z = 0;
		for(i= global::prt_list.begin(); i != global::prt_list.end(); i++)
		{
			v.x = (2*i->pos.x) - (i->posp.x) + (i->force.x*dT_sq);
			i->vel.x = (v.x - (i->posp.x))/(2*global::dT);
			global::comv.x += i->vel.x;
			i->posp.x = i->pos.x;
			i->pos.x = v.x;
			i->dist.x += i->pos.x - i->posp.x;

			v.y = (2*i->pos.y) - (i->posp.y) + (i->force.y*dT_sq);
			i->vel.y = (v.y - (i->posp.y))/(2*global::dT);
			global::comv.y += i->vel.y;
			i->posp.y = i->pos.y;
			i->pos.y = v.y;
			i->dist.y += i->pos.y - i->posp.y;
			
			v.z = (2*i->pos.z) - (i->posp.z) + (i->force.z*dT_sq);
			i->vel.z = (v.z - (i->posp.z))/(2*global::dT);
			global::comv.z += i->vel.z;
			i->posp.z = i->pos.z;
			i->pos.z = v.z;
			i->dist.z += i->pos.z - i->posp.z;

			sumv_sq += i->vel.magnitude();
		}
		global::temp = sumv_sq/(3*global::n_prt); 
		global::etot = (global::pe + 0.5*sumv_sq)/global::n_prt;
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
		
		if(dist_max > global::skin*global::skin)
		{
			compute::refreshList();
			for(i= global::prt_list.begin(); i != global::prt_list.end(); i++) i->dist = 0;
		}
		compute::forces();
		compute::verlet();
		compute::pbc();
	}
}
#endif


