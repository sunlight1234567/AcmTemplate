/*
Computing Geometry
计算几何需要充分灵活运用数学解析几何知识
*/
#include <bits/stdc++.h>
using namespace std;
using bd = double;
using ldb = long double;

class point
{
public:
    ldb x, y;
    point(){}
    point(ldb _x, ldb _y){x=_x; y=_y;}
    ldb dis(){
        return sqrt(x*x+y*y);
    }
    ldb dis(const point& a){
        return sqrt((x-a.x)*(x-a.x) + (y-a.y)*(y-a.y));
    }
    point operator+ (const point& a){
        return {x+a.x, y+a.y};
    }
    point operator- (const point& a){
        return {x-a.x, y-a.y};
    }
    point operator/(const ldb& k){
        return {x/k, y/k};
    }
};

int main()
{
    return 0;
}