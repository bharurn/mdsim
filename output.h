#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <string>

namespace output
{
	std::ofstream f, d;
	real a,b,c;
	int n;

	void trajectory()
	{
		static int t = 0;
		f << global::n_prt << "\n\n";
		for(int i=0;i<global::n_prt;i++)
			f << "H " << global::prt_list[i].pos.x << " " << global::prt_list[i].pos.y << " " << global::prt_list[i].pos.z << " " 
			<< global::prt_list[i].vel.x << " " << global::prt_list[i].vel.y << " " << global::prt_list[i].vel.z << '\n';
		t++;
	}

	void log(real t)
	{
		d << t << "		" << global::etot << "		" << global::temp << '\n'; 
	}
	
	void setDumpFile(std::string str1, std::string str2)
	{
		f.open(str1);
		d.open(str2);
		d << "All quantities in LJ units.\nt		Etot		T\n"; 
	}
}

#endif
