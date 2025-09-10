#pragma once

class Point
{
public:
	Point(double x, double y);

	// define getters
	double getX() { return _x};
	double getY() {return _y};

private:
	double _x;
	double _y;
