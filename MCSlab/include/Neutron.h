#pragma once

class Neutron
{
public:
	Neutron(Point position, double angle);
	
	
private:
	Point _pos;
	double _ang;

	double randomIsoAngle();
};	
