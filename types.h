#ifndef TYPES_H
#define TYPES_H

typedef double real;

class vec
{
public:
	real x,y,z;
	vec(){}
	vec(real a, real b, real c)
	{
		x = a; y = b; z = c;
	}
	
	real magnitude()
	{
		return ((x*x) + (y*y) + (z*z));
	}
	
	void operator=(const int &v )
    { 
    	this->x = v;
        this->y = v;
        this->z = v;
    }
    
    void operator=(vec &v )
    { 
    	this->x = v.x;
        this->y = v.y;
        this->z = v.z;
	}
};

class atom
{
public:
	vec pos, vel, force, posp, dist;
	int ix=0, iy=0, iz=0;
	atom() {}
	atom(vec v, vec ve)
	{
		pos = v;
		vel = ve;
		force = 0.0;
		dist = 0.0;
		ix = 0;
		iy = 0;
		iz = 0;
	}
};

#endif
