#include <bits/stdc++.h>
using namespace std;

struct Point {
    double x, y;
    bool operator < (const Point &a) const {if (x != a.x) return x < a.x;return y < a.y;}
    bool operator == (const Point &a) const {return x == a.x && y == a.y;}
    double mod() { return sqrt(x * x + y * y); }
    Point operator - (const Point &a) const {return {x - a.x, y - a.y};}
    double cross(const Point &a) const {return x * a.y - y * a.x;}
};


void convexHull() {
    int n;
    cin >> n;
    vector<Point> a(n);
    for (int i = 0; i < n; ++i) cin >> a[i].x >> a[i].y;

    if (n == 1) {
        cout << "0.00\n";
        return;
    }

    // 排序并去重
    sort(a.begin(), a.end());
    a.erase(unique(a.begin(), a.end()), a.end());
    n = a.size();

    vector<Point> hull;
    // 构建下凸包
    for (int i = 0; i < n; ++i) {
        while (hull.size() >= 2 && 
               (hull.back() - hull[hull.size()-2]).cross(a[i] - hull.back()) <= 0) {
            hull.pop_back();
        }
        hull.push_back(a[i]);
    }
    // 构建上凸包
    for (int i = n - 2, t = hull.size(); i >= 0; --i) {
        while (hull.size() > t && 
               (hull.back() - hull[hull.size()-2]).cross(a[i] - hull.back()) <= 0) {
            hull.pop_back();
        }
        hull.push_back(a[i]);
    }

    // 计算周长（凸包是闭合的，hull[0] 和 hull.back() 是同一个点）
    double ans = 0;
    for (int i = 1; i < hull.size(); ++i) {
        ans += (hull[i] - hull[i-1]).mod();
    }
    cout << fixed << setprecision(2) << ans << '\n';
}

// int main() {
//     ios::sync_with_stdio(0);
//     cin.tie(0), cout.tie(0);
//     convexHull();
//     return 0;
// }