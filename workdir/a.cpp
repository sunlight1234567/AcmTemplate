#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const int MAXN = 2e5 + 10;
int n, m;

void check(int n) {
    vector<int> a(n);
    for(int i=0; i<n; ++i) a[i] = i;
    while(next_permutation(a.begin(), a.end())){
        for(int i=0; i<n; ++i) {

        }
    }
    
}
void sol() {
    cin >> n;
    if(n <= 6) {
        if(n == 2) {cout << -1 << '\n';}
        else if(n == 3) {
            cout << "1 3\n2 3\n";
        }
        else if(n==4) {
            cout << "1 2\n3 1\n4 1\n";
        }
        else if(n == 5) {
            cout << "2 5\n3 5\n4 5\n1 4\n";
        }
        else if(n == 6) {
            cout << "1 6\n5 6\n 2 5\n3 5\n4 5\n";
        }
        // else if(n == 7) {
        //     cout << "1 7\n2 7\n3 6\n4 7\n5 7\n6 7\n";
        // }
        // else if(n == 8) {
        //     cout << "1 8\n2 8\n3 8\n4 1\n5 8\n6 8\n7 8\n";
        // }
        // else if(n == 9) {
        //     cout << "1 9\n2 9\n3 9\n4 9\n 5 9\n6 9\n7 9\n8 9\n";
        // }
        return;

    }
    
    for(int i=1; i<=n; ++i) {
        if(i==2) continue;
        if(i==n-4) cout << i << " 1\n";
        else cout << i << " 2\n";
    }
}
signed main() {
    ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
    int T=1;
    cin>>T;
    while(T--){sol();}
    return 0;
}