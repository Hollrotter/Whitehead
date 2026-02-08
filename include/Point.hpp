#pragma once

class Point
{
    double x;
    double y;
public:
    Point() = default;
    Point(double _x, double _y) : x(_x), y(_y){};
    double X() const
    {
        return x;
    }
    double Y() const
    {
        return y;
    }
};