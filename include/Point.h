#pragma once

class Point
{
public:
  Point(double x, double y)
  {
    _x = x;
    _y = y;
  };

  // define getters
  double getX() { return _x; }
  double getY() { return _y; }

  // define setter functions
  void moveX(double x) { _x = x; }
  void moveY(double y) { _y = y; }

private:
  double _x;
  double _y;
};
