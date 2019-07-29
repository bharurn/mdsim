#ifndef OUTPUT_H
#define OUTPUT_H

#define BOUNDBOX 0

#include <fstream>
#include <string>

namespace output
{
	std::ofstream f, d;
	real a,b,c;
	int n;

	void boundBox()
	{
			n = global::n_prt+1;

			for(int i=0; i<2; i++)
			{
				if(i%2 == 0) a = 0.5*global::input["xlen"];
				else a = -0.5*global::input["xlen"];
				for(int j=0; j<2; j++)
				{
					if(j%2 == 0) b = 0.5*global::input["ylen"];
					else b = -0.5*global::input["ylen"];
					for(int k=0; k<2; k++)
					{
						if(k%2 == 0) c = 0.5*global::input["zlen"];
						else c = -0.5*global::input["zlen"];
						f << "C " << a <<  " " << b << " " << c << '\n';
						++n;
					}
				}
			}
			n = global::n_prt+1;
	}

	void trajectory()
	{
		static int t = 0;
		if(BOUNDBOX) f << global::n_prt+8 << "\n\n";
		else f << global::n_prt << "\n\n";
		for(int i=0;i<global::n_prt;i++)
			f << "H " << global::prt_list[i].pos.x << " " << global::prt_list[i].pos.y << " " << global::prt_list[i].pos.z << " " 
			<< global::prt_list[i].vel.x << " " << global::prt_list[i].vel.y << " " << global::prt_list[i].vel.z << '\n';
		if(BOUNDBOX) boundBox();
		t++;
	}

	void log(real t)
	{
		d << t << "		" << global::etot << "		" << global::input["temp"] << '\n'; 
	}
	
	void setDumpFile(std::string str1, std::string str2)
	{
		f.open(str1);
		d.open(str2);
		d << "All quantities in LJ units.\nt		Etot		T\n"; 
	}
}

#endif
