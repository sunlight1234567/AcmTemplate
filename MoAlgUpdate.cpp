//带修改莫队
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
int MAXN = 140000;
int MAXM = 140000;
int unitsize = pow(MAXN, 2.0/3.0);
ll n, m;
struct Query{
    int l, r, id, t;
    bool operator < (const Query &x) const {
        if( l/unitsize != x.l/unitsize ) return l/unitsize < x.l/unitsize;
        if( r/unitsize != x.r/unitsize ) return r/unitsize < x.r/unitsize;
        return t < x.t;
    }
};
struct Update{
    int pos, upvalue;
};
vector<int> arr(MAXN, 0);
vector<Query> qs(MAXM);
vector<Update> upds(MAXM);
vector<ll> freq(1e6+10, 0);
vector<int> ans(MAXM, 0);
int kind = 0;
void add(int x) {
    if(freq[x]++ == 0) ++kind;
}
void del(int x) {
    if(freq[x]-- == 1) --kind;
}
void move(int pos, int type) {    
    int x = arr[pos];
    if(type == 1) { // add   
        if(freq[x]++ == 0) ++kind;
    }
    else if(type == -1){ // delete  
        if(freq[x]-- == 1) --kind;
    }
} 
void MoAlg() {
    cin >> n >> m;
    unitsize = pow(n, 2.0/3.0);
    for(int i=1; i<=n; ++i) cin >> arr[i];
    int s1 = 0, s2 = 0;//querysize,updatesize
    int tick = 0;
    for(int i=0; i<m; ++i) {
        char c;int l,r;cin >> c >> l >> r;
        if(c == 'Q'){
            qs[s1] = {l,r,s1,tick};
            s1++;
        }
        else {
            upds[++s2] = {l,r};++tick;
        }
    }
    sort(qs.begin(), qs.begin()+s1);
    int l=1, r=0, t=0;
    for(int i=0; i<s1; ++i) {
        auto [ql, qr, qid, qt] = qs[i];
        while(l > ql) move(--l, 1);
        while(r < qr) move(++r, 1);
        while(l < ql) move(l++, -1);
        while(r > qr) move(r--, -1);
        while(t < qt) {
            ++t;
            auto [pos, val] = upds[t];
            if(pos>=l && pos<=r) {
                del(arr[pos]);
                add(val);
            }
            swap(arr[pos], upds[t].upvalue);
        }
        while(t > qt) {
            auto [pos, val] = upds[t];
            if(pos>=l && pos<=r) {
                del(arr[pos]);
                add(val);
            }
            swap(arr[pos], upds[t].upvalue);
            --t;//这里应该先撤销修改再移动t
        }
        ans[qid] = kind;
    }
    for(int i=0; i<s1; ++i) cout << ans[i] << '\n';
    
    
    
}
int main() {
    MoAlg();
    return 0;
}
