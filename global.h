#ifndef GLOBAL_H
#define GLOBAL_H

#include "types.h"
#include <cmath>
#include <vector>
#include <map> 

namespace global
{
	real pe, etot, sumv_sq=0, set_temp;
	int n_prt = 1000, T_f, dump_time;
	
	std::map<std::string, real> input;
	
	vec comv(0.0, 0.0, 0.0);
	real ecut;
	std::vector<atom> prt_list(0);
}

#endif
