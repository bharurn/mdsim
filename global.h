#ifndef GLOBAL_H
#define GLOBAL_H

#include "types.h"
#include <cmath>
#include <vector>

namespace global
{
	real xlen = 100, ylen = 100, zlen = 100, dT = 0.003;
	real temp, pe, etot;
	int T_f = 1000, n_prt = 1000;
	vec comv(0.0, 0.0, 0.0);
	real cut = 2.5, ecut = 4 * pow(cut, -6) * (pow(cut, -6) - 1), skin = 0.3;
	std::vector<atom> prt_list(0);
	int dump_time = 1;
}

#endif
