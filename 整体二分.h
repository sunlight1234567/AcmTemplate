//RangeKth
#include<bits/stdc++.h>
using namespace std;

struct node {
    int pos, num;
    bool operator <(const node &x)const{
        if(num != x.num)
            return num < x.num;
        return pos < x.pos;
    }
};
const int N = 10 + 2e5;
int n, m;
node arr[N] = {0};
//查询
int qid[N];
int l[N];
int r[N];
int k[N];

int ans[N];

int tree[N];//树状数组
void add(int i, int v) {
    while (i <= n) {
        tree[i] += v;
        i += (i & (-i));
    }
}
int sum(int i) {
    int res = 0;
    while (i > 0) {
        res += tree[i];
        i -= (i & (-i));
    }
    return res;
}
int query(int l, int r) {
    return sum(r) - sum(l-1);
}

void compute(int ql, int qr, int vl, int vr) {
    if(ql > qr) return;
    if(vl == vr) {
        for(int i=ql; i<=qr; ++i) {
            ans[qid[i]] = arr[vl].num; 
        }
    }
    else {
        int mid = (vl + vr) / 2;
        for(int i=vl; i<=mid; ++i) {
            add(arr[i].pos, 1);
        }
        vector<int> lset, rset;
        for(int i=ql; i<=qr; ++i) {
            int id = qid[i];
            int cnt = query(l[id], r[id]);
            
            if(cnt >= k[id]) {
                lset.push_back(id);
            }
            else if(cnt < k[id]) {
                rset.push_back(id);
                k[id] -= cnt;
            }
        }
        for(int i=0; i<lset.size(); ++i) qid[ql + i] = lset[i];
        for(int i=0; i<rset.size(); ++i) qid[ql + lset.size() + i] = rset[i];

        for(int i=vl; i<=mid; ++i) {
            add(arr[i].pos, -1);
        }

        compute(ql, ql + lset.size() - 1, vl, mid);
        compute(ql + lset.size(), qr, mid + 1, vr);
    }
}
int main() {
    cin >> n >> m;
    for(int i=1; i<=n; ++i) {
        arr[i].pos = i;
        cin >> arr[i].num;
    }
    for(int i=1; i<=m; ++i) {
        cin >> l[i] >> r[i] >> k[i];
        qid[i] = i;
    }
    sort(arr+1, arr+1+n, less<node>());
    compute(1, m, 1, n);
    for(int i=1; i<=m; ++i)cout << ans[i] << endl;
    return 0;
}