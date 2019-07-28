#ifndef MEAN_R_SQ_H
#define MEAN_R_SQ_H

#include "global.h"

namespace mean_r_sq
{
	struct diff_value_datatype
	{
		std::vector<vec> ori;
		real diff;
	};
	
	int time_length = 2;
	int max_time = global::T_f/time_length;
	
	std::vector<diff_value_datatype> diff_value;
	
	std::ofstream f;
	
	void init()
	{
		std::vector<atom>::iterator i;
		diff_value_datatype d;		
		
		for(int k=0;k<max_time; k++)
		{
			for(i = global::prt_list.begin(); i != global::prt_list.end(); i++)
			{
				d.ori.push_back(i->pos);
			}
			d.diff = 0;
			diff_value.push_back(d);
		}
		
		f.open("mean_r_sq_diff.txt");
		f << "All quantities in LJ units.\nt		D\n"; 
	}
	
	void calculate(int time)
	{
		std::vector<vec>::iterator origin;
		std::vector<atom>::iterator j;
		vec disp, dr;
		
		for(int i=0;i<max_time; i++)
		{
			if(time%((i+1)*time_length) == 0)
			{
				for(int j=0; j<global::n_prt; j++)	
				{	
					disp.x = global::prt_list[j].pos.x + global::prt_list[j].ix*global::xlen;
					disp.y = global::prt_list[j].pos.y + global::prt_list[j].iy*global::ylen;
					disp.z = global::prt_list[j].pos.z + global::prt_list[j].iz*global::zlen;
					
					dr.x = disp.x - diff_value[i].ori[j].x;
					dr.y = disp.y - diff_value[i].ori[j].y;
					dr.z = disp.z - diff_value[i].ori[j].z;
					
					diff_value[i].diff += dr.magnitude();
					diff_value[i].ori[j].x = disp.x;
					diff_value[i].ori[j].y = disp.y;
					diff_value[i].ori[j].z = disp.z;
				}
			}
		}
	}
	
	void output()
	{
		int no;
		for(int i=0; i<max_time; i++)
		{
			no = (int)(global::T_f/((i+1)*time_length));
			diff_value[i].diff /= (global::n_prt*no);
		}
		
		int val = global::T_f/global::dump_time;
		
		for(int i=1; i<max_time; i+=val)
			f << i*global::dT*time_length << "		" << diff_value[i].diff << '\n';
	}
}

#endif
