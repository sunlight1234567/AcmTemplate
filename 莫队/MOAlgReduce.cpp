/*
回滚莫队
https://www.luogu.com.cn/problem/P1997
https://loj.ac/p/2874 oiwiki例题
区间增加时可以维护答案信息，减少时不可以维护
或者反之
*/
#include <bits/stdc++.h>
#define int long long
using namespace std;
using ll = long long;
int MAXN = 1e5 + 10;
int MAXM = 2e5 + 10;
int unitsize = pow(MAXN, 1.0/2.0);
int n, m;
struct Query{
    int l, r, id;
    bool operator < (const Query &x) const {
        //下标从0开始的块号应该为(ind-1)/unitsize
        //普通莫队将块号当作(ind)/unitsize影响不大，但是回滚莫队对分块要求高
        if( (l-1)/unitsize == (x.l-1)/unitsize ) return r < x.r;
        return (l-1)/unitsize < (x.l-1)/unitsize;
    }
};

vector<int> bi(MAXN), br(MAXN);//分块信息:所在块号，右边界
int blen, bnum;
vector<int> arr(MAXN, 0);
vector<Query> qs(MAXM);
vector<int> freq(MAXN, 0);
vector<int> ans(MAXM, 0);
vector<int> decode(MAXN, 0);//离散化复原
int curans = 0;
int force(int l, int r) {
    vector<int> freq(MAXN, 0);
    int res = 0 ;
    for(int i=l; i<=r; ++i) {
        res = max(res, ++freq[arr[i]] * decode[arr[i]]);
    }
    return res;
}
void add(int x){
    //curans = max(curans, ++freq[x]);
    curans = max(curans, ++freq[x] * decode[x]);
}
void del(int x){
    --freq[x];
}
void MoAlg() {
    cin >> n >> m;
    unitsize = sqrt(n);
    for(int i=1; i<=n; ++i) {
        cin >> arr[i];//按升序排列
    }
    for(int i=0; i<m; ++i) {
        int l,r;
        cin >> l >> r;
        qs[i] = {l,r,i};
    }
    //离散化a{}
    vector<int> unic(arr.begin()+1, arr.begin()+1+n);
    sort(unic.begin(), unic.end());
    unic.resize( unique(unic.begin(), unic.end()) - unic.begin());
    for(int i=1; i<=n; ++i) {
        int t = arr[i];
        arr[i] = lower_bound(unic.begin(), unic.end(), arr[i])-unic.begin();
        decode[arr[i]] = t;
    }
    blen = unitsize;
    bnum =  (n + blen - 1) / blen;
    for(int i=1; i<=n; ++i) {
        bi[arr[i]] = (i - 1) / blen;
    }
    for(int i=0 ;i<bnum; ++i) br[i] = min(n, (i+1)*blen);
    
    sort(qs.begin(), qs.begin()+m);
    int l=1, r=0;
    for(int bid=0, i=0; bid<bnum && i<m; ++bid) {
        int en = min((n),(bid+1)*blen);
        while(i<m && qs[i].r<=en) {
            auto [ql, qr, qid] = qs[i];
            ans[qid] = force(ql, qr);
            ++i;
        }
        l=en+1, r=en;
        freq = vector<int>(MAXN, 0);//注意对上一块信息的清空
        curans = 0;
        int backup = 0;
        while(i<m && qs[i].l<=en){
            auto [ql, qr, qid] = qs[i];
            while(r < qr) add(arr[++r]);
            backup = curans;
            while(l > ql) add(arr[--l]);
            ans[qid] = curans;
            while(l < en+1) del(arr[l++]);
            curans = backup;
            ++i;
        }
    }
    for(int i=0; i<m; ++i) cout << ans[i] << '\n';
}
signed main() {
    MoAlg();
    return 0;
}
