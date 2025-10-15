#include<bits/stdc++.h>

struct Point {
    double x, y;
    bool operator < (const Point &a) const {if (x != a.x) return x < a.x;return y < a.y;}
    bool operator == (const Point &a) const {return x == a.x && y == a.y;}
    double mod() { return sqrt(x * x + y * y); }
    Point operator - (const Point &a) const {return {x - a.x, y - a.y};}
    double cross(const Point &a) const {return x * a.y - y * a.x;}
};
