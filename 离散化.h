#include<bits/stdc++.h>
using namespace std;
void fun(){
    int n;
    vector<int> a(n);
    for(auto &x : a) cin >> x;
    vector<int> unic = a;
    sort(unic.begin(), unic.end());
    unic.resize(unique(unic.begin(), unic.end()) - unic.begin());
    for(int i=0; i<n; ++i) {
        a[i] = lower_bound(unic.begin(), unic.end(), a[i]) - unic.begin();
    }
}
