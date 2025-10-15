#include <iostream>
using namespace std;
using ll= long long;
//带权并查集

const int N = 10 + 2e5;
ll n, m;
ll d[N];
ll s[N];
ll ans = 0;
ll find(ll x) {
    if(s[x] != x){
        int t = s[x];
        s[x] = find(s[x]);
        d[x] += d[t];
    }
    return s[x];
}
void merge(ll x, ll y, ll _d) {// x <= y       x->y
    ll fx = find(x), fy = find(y);
    if(fx != fy) {
        s[fx] = fy;
        d[fx] = d[y] - d[x] + _d;
    }
    else {
        if(d[x] - d[y] != _d) ++ans;
    }
}
int main() {
    std :: cin >> n >> m;
    //初始化？
    for(int i=0; i<=n; ++i) {s[i] = i; d[i] = 0;}
    while(m--) {
        ll a, b, v;
        cin >> a >> b >> v;//a <= b
        a--;
        merge(a, b, v);
    }
    cout << ans; 
}