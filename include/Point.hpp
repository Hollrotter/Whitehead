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
    void perpendicular()
    {
        std::swap(x, y);
        x *=-1;
    }
    inline Point &operator+=(const Point &P)
    {
        x += P.x;
        y += P.y;
        return *this;
    }
    inline Point &operator-=(const Point &P)
    {
        x -= P.x;
        y -= P.y;
        return *this;
    }
    template <Number C> Point &operator*=(const C c)
    {
        x *= c;
        y *= c;
        return *this;
    }
    template <Number C> Point &operator/=(const C c)
    {
        x /= c;
        y /= c;
        return *this;
    }

};

inline double distance(Point P1, Point P2)
{
    return sqrt(pow(P1.X()-P2.X(), 2) + pow(P1.Y()-P2.Y(), 2));
}

inline Point operator+(Point P1, const Point &P2)
{
    return P1 += P2;
}

inline Point operator-(Point P1, const Point &P2)
{
    return P1 -= P2;
}

template <Number C> Point operator*(Point P, C c)
{
    return P *= c;
}

template <Number C> Point operator*(C c, Point P)
{
    return P *= c;
}

template <Number C> Point operator/(Point P, C c)
{
    return P /= c;
}