//倍增https://www.luogu.com.cn/problem/P4155
#include<bits/stdc++.h>
using namespace std;

struct node{
    int x, y, id;
    bool operator<(const node &_n) {return pair{x, y} < pair{_n.x, _n.y};}
};
void sol() {
    int n,m;
    cin >> n >> m;
    vector<node> a(n);
    for(int i=0; i<n; ++i) {
        cin >> a[i].x >> a[i].y;
        a[i].id = i;
        if(a[i].y < a[i].x) a[i].y += m;
    }
    sort(a.begin(), a.end());
    vector<vector<int>> st(2*n, vector<int>(log2(2*n)+1, 0));
    for(int i=0; i<n; ++i) {
        a.push_back({a[i].x + m, a[i].y + m});
    }
    

    for(int i=0, j=0; i<2*n; ++i) {
        while( j+ 1 < 2*n && a[j+1].x <= a[i].y) ++j;
        if(a[j].x <= a[i].y && j > i)st[i][0] = j;
    }
    for(int j=1; j<=log2(2*n); ++j) {
        for(int i=0; i<2*n; ++i) {
            int mid = st[i][j-1];
            if(mid > i)
                st[i][j] = st[mid][j-1];
        }
    }
    // for(int i=0; i<2*n; ++i) {
    //     for(int j=0; j<=log2(2*n); ++j)cout <<st[i][j] << ' ';
    //     cout << endl;
    // }
    vector<int> ans(n, 0);
    for(int i=0; i<n; ++i) {
        int p = log2(n);//最多跳n步
        int en = i + n;
        int now = i;
        while(now < en) {
            while(p && st[now][p] == 0) --p;
            while(p && st[now][p] > en) {
                --p;
            }
            now = st[now][p];
            ans[a[i].id] += 1 << p;
        }
    }
    for(int i=0; i<n; ++i)cout << ans[i] << ' ';
}

//ST表维护区间最值https://www.luogu.com.cn/problem/P2880
const int inf = 10 + 1e6;
void sol() {
    int n, q;
    cin >> n >> q;
    vector<int> a(n), l(q), r(q);
    for(int i=0; i<n; ++i) cin >> a[i];
    for(int i=0; i<q; ++i) cin >> l[i] >> r[i];
    vector<vector<int>> stmax(n, vector<int>(log2(n) + 1, 0)), stmin(n, vector<int>(log2(n) + 1, inf));
    //init
    for(int i=0; i<n; ++i) {
        stmax[i][0] = a[i];
        stmin[i][0] = a[i];
    }
    for(int j=1; j<=log2(n); ++j) {
        for(int i=0; i<n; ++i) {
            if(i + (1 << j) >= n) continue;
            int ni = i + (1 << (j-1));
            stmax[i][j] = max(stmax[ni][j-1], stmax[i][j-1]);
            stmin[i][j] = min(stmin[ni][j-1], stmin[i][j-1]);
        }
    }  
    vector<int> ansmax(q, 0), ansmin(q, inf);
    for(int i=0; i<q; ++i) {
        //max
        int ql = l[i]-1, qr = r[i]-1;
        int now = ql;
        int p = 0;
        while(now <= qr) {
            //相较于从大到小推，小到大推P不断置0，会有重复的运算
            p = 0;
            while(now + (1 << (p + 1) ) <= qr) ++p;
            ansmax[i] = max(ansmax[i], stmax[now][p]);
            now += 1 << p;
        }
        //min
        p = 0;
        now = ql;
        while(now <= qr) {
            p = 0;
            while(now + (1 << (p + 1) ) <= qr) ++p;
            ansmin[i] = min(ansmin[i], stmin[now][p]);
            now += 1 << p;
        }
    }
    for(int i=0; i<q; ++i) cout << ansmax[i] - ansmin[i] << endl;

}

//维护gcd https://www.luogu.com.cn/problem/P1890
int gcd (int a, int b) {
    while(b) {
        a %= b;
        swap(a, b);
    }
    return a;
}
const int inf = 10 + 1e6;
void sol() {
    int n, q;
    cin >> n >> q;
    vector<int> a(n), l(q), r(q);
    for(int i=0; i<n; ++i) cin >> a[i];
    for(int i=0; i<q; ++i) cin >> l[i] >> r[i];
    vector<vector<int>> stgcd(n, vector<int>(log2(n) + 1, 0));
    //init
    for(int i=0; i<n; ++i) {
        stgcd[i][0] = a[i];
    }
    for(int j=1; j<=log2(n); ++j) {
        for(int i=0; i<n; ++i) {
            if(i + (1 << j) >= n) continue;
            int ni = i + (1 << (j-1));
            if(stgcd[i][j-1] && stgcd[ni][j-1]) 
                stgcd[i][j] = gcd(stgcd[i][j-1], stgcd[ni][j-1]);
        }
    }  
    vector<int> ans(q, 0);
    for(int i=0; i<q; ++i) {
        //max
        int ql = l[i]-1, qr = r[i]-1;
        int now = ql;
        int p = 0;
        ans[i] = stgcd[ql][0];
        while(now <= qr) {
            //相较于从大到小推，小到大推P不断置0，会有重复的运算
            p = 0;
            while(now + (1 << (p + 1) ) <= qr) ++p;
            ans[i] = gcd(ans[i], stgcd[now][p]);
            now += 1 << p;
        }
        
    }
    for(int i=0; i<q; ++i) cout << ans[i] << endl;

}

//倍增LCA https://www.luogu.com.cn/problem/P3339
void sol() {
    int n, m, root;
    cin >> n >> m >> root;
    vector<int> g[n+1], dep(n+1, -1);
    vector<int> qa(m), qb(m), ans(m);
    for(int i=0; i<n-1; ++i) {
        int u, v;
        cin >> u >> v;
        g[u].push_back(v);
        g[v].push_back(u);
    }
    for(int i=0; i<m; ++i) {
        cin >> qa[i] >> qb[i];
    }
    vector<vector<int>> stjump(n+1, vector<int>( log2(n) + 1, -1));
    function<void(int)> dfs = [&] (int cur){
        for(auto nx : g[cur]) {
            if(dep[nx] == -1) {
                dep[nx] = dep[cur] + 1;
                stjump[nx][0] = cur;
                dfs(nx);
            }
        }
    };
    dep[root] = 0;
    dfs(root);

    for(int j=1; j<=log2(n); ++j) {
        for(int i=1; i<=n; ++i) {
            int nx = stjump[i][j-1];
            if(nx >= 0) stjump[i][j] = stjump[nx][j-1];
        }
    }
    auto jump = [&](int _a, int k) -> int{
        while(k) {
            int p = 0;
            while(1 << p + 1 < k) ++ p;
            _a = stjump[_a][p];
            k -= 1 << p;
        }
        return _a;
    };
    for(int i=0; i<m; ++i) {
        int a = qa[i], b = qb[i];
        //int p = log2(max(dep[a], dep[b]));
        if(dep[a] > dep[b]) {
            a = jump(a, dep[a] - dep[b]);
        }
        else if(dep[a] < dep[b]) b = jump(b, dep[b] - dep[a]);
        if(a == b){ 
            ans[i] = a;
            continue;
        }

        while(stjump[a][0] != stjump[b][0]) {
            int p = 0;
            while(dep[a] >= 1<<p && stjump[a][p+1] != stjump[b][p+1]) ++p;
            a = stjump[a][p];
            b = stjump[b][p];
        }
        ans[i] = stjump[a][0];
    }
    for(int i=0; i<m; ++i) cout << ans[i] << '\n';
}