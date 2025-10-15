/*
经典莫队算法
https://www.luogu.com.cn/problem/SP3267
https://www.luogu.com.cn/problem/P1494
https://www.luogu.com.cn/problem/CF617E
https://www.luogu.com.cn/problem/P3245
*/
#include <bits/stdc++.h>
using namespace std;
int MAXN = 50010;
int MAXM = 50010;
int unitsize = sqrt(MAXN);
struct Query{
    int l, r, id;
    bool operator < (const Query &x) const {
        if( l/unitsize != x.l/unitsize ) return l/unitsize < x.l/unitsize;
        if((l/unitsize) & 1) return r < x.r;
        else return r > x.r;
    }
};

vector<int> arr(MAXN, 0);
vector<Query> qs(MAXM);
vector<int> freq(1e6+10, 0);
int kind = 0;
long long A=0, B = 1;
void move(int pos, int type) {
    if(type == 1) { // add   
        //if(freq[arr[pos]]++ == 0) ++kind;
        if(freq[arr[pos]] == 0) ++kind;
        freq[arr[pos]] ++;
        long long t = freq[arr[pos]];
        A -= 1LL*(t-1)*(t-2);
        A += 1LL*t*(t-1);
    }
    else { // delete
        if(freq[arr[pos]] == 1) --kind;
        freq[arr[pos]] --;
        long long t = freq[arr[pos]];
        A -= 1LL*(t+1)*(t);
        A += 1LL*t*(t-1);
    }
} 
void MoAlg() {
    int n, m;
    cin >> n >> m;
    for(int i=0; i<n; ++i) cin >> arr[i];
    for(int i=0; i<m; ++i) {
        int l,r;cin >> l >> r;
        qs[i] = {l-1, r-1, i};
    }
    sort(qs.begin(), qs.begin()+m);
    int l = 0, r = -1;
    
    vector<pair<long long, long long>> ans(m); 
    for(int i=0; i<m; ++i) {
        auto [ql, qr, qid] = qs[i];
        if(ql == qr) {ans[qid] = {0, 1}; continue;}
        while(l > ql) move(--l, 1);
        while(r < qr) move(++r, 1);
        while(l < ql) move(l++, -1);
        while(r > qr) move(r--, -1);
        
        B = 1LL*(r-l+1)*(r-l);
        long long g = __gcd(A,B);
        ans[qid] = {A/g, B/g};
    }
    for(int i=0; i<m; ++i) cout << ans[i].first << '/' << ans [i].second << '\n';
    
}
int main() {
    MoAlg();
    return 0;
}