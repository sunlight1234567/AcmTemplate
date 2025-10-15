// 树上莫队入门题，C++版
// 一共有n个节点，每个节点给定颜色值，给定n-1条边，所有节点连成一棵树
// 一共有m条查询，格式 u v : 打印点u到点v的简单路径上，有几种不同的颜色
// 1 <= n <= 4 * 10^4
// 1 <= m <= 10^5
// 1 <= 颜色值 <= 2 * 10^9
// 测试链接 : https://www.luogu.com.cn/problem/SP10707
// 测试链接 : https://www.spoj.com/problems/COT2/
#include <bits/stdc++.h>

using namespace std;

struct Query {
   int l, r, lca, id;
};

struct Update {
   int pos, val;
};

const int MAXN = 100001;
const int MAXP = 20;
int n, m;
int color[MAXN];
int decode[MAXN];
int kind = 0;
vector<int> graph[MAXN];

Query query[MAXN];

int dep[MAXN];
int seg[MAXN << 1];
int st[MAXN];
int ed[MAXN];
int stjump[MAXN][MAXP];
int cntd;

int bi[MAXN << 1];
bool vis[MAXN];
int cnt[MAXN];
long long curans;
long long ans[MAXN];

void dfs(int u, int fa) {
    dep[u] = dep[fa] + 1;
    seg[++cntd] = u;
    st[u] = cntd;
    stjump[u][0] = fa;
    for (int p = 1; p < MAXP; p++) {
        stjump[u][p] = stjump[stjump[u][p - 1]][p - 1];
    }
    // 修改：遍历邻接表的方式
    for (int v : graph[u]) {
        if (v != fa) {
            dfs(v, u);
        }
    }
    seg[++cntd] = u;
    ed[u] = cntd;
}

int lca(int a, int b) {
    if (dep[a] < dep[b]) {
        swap(a, b);
    }
    for (int p = MAXP - 1; p >= 0; p--) {
        if (dep[stjump[a][p]] >= dep[b]) {
            a = stjump[a][p];
        }
    }
    if (a == b) {
        return a;
    }
    for (int p = MAXP - 1; p >= 0; p--) {
        if (stjump[a][p] != stjump[b][p]) {
            a = stjump[a][p];
            b = stjump[b][p];
        }
    }
    return stjump[a][0];
}

bool QueryCmp(Query &a, Query &b) {
    if (bi[a.l] != bi[b.l]) {
        return bi[a.l] < bi[b.l];
    }
    return (bi[a.r] & 1) ^ (a.r < b.r);
}

void invert(int node) {
    int col = color[node];
    if (vis[node]) {
        if(--cnt[col] == 0) --kind;
    } else {
        if(++cnt[col] == 1) ++kind; 
    }
    vis[node] = !vis[node];
}

void compute() {
    int winl = 1, winr = 0, wint = 0;
    for (int i = 1; i <= m; i++) {
        int jobl = query[i].l;
        int jobr = query[i].r;
        int lca = query[i].lca;
        int id = query[i].id;
        while (winl > jobl) {
            invert(seg[--winl]);
        }
        while (winr < jobr) {
            invert(seg[++winr]);
        }
        while (winl < jobl) {
            invert(seg[winl++]);
        }
        while (winr > jobr) {
            invert(seg[winr--]);
        }
        if (lca > 0) {
            invert(lca);
        }
        ans[id] = kind;
        if (lca > 0) {
            invert(lca);
        }
    }
}

void prapare() {
    int blen = sqrt(cntd);
    for (int i = 1; i <= cntd; i++) {
        bi[i] = (i - 1) / blen + 1;
    }
    sort(query + 1, query + m + 1, QueryCmp);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cin >> n >> m;
    for (int i = 1; i <= n; i++) {
        cin >> color[i];
    }
    vector<int> unic(n);
    for(int i=1; i<=n; ++i) unic[i-1] = color[i];
    sort(unic.begin(), unic.end());
    unic.resize(unique(unic.begin(),unic.end())-unic.begin());
    for(int i=1; i<=n; ++i) {
        int t = color[i];
        color[i] = lower_bound(unic.begin(), unic.end(), color[i])-unic.begin();
        decode[color[i]] = t;
    }
    for (int i = 1, u, v; i < n; i++) {
        cin >> u >> v;
        graph[u].push_back(v);
        graph[v].push_back(u);
    }
    dfs(1, 0);
    for (int i = 1, op, x, y; i <= m; i++) {
        cin >> x >> y;
        
        if (st[x] > st[y]) {
            swap(x, y);
        }
        int xylca = lca(x, y);
        if (x == xylca) {
            query[i] = {st[x], st[y], 0, i};
        } else {
            query[i] = {ed[x], st[y], xylca, i};
        }
        
    }
    prapare();
    compute();
    for (int i = 1; i <= m; i++) {
        cout << ans[i] << '\n';
    }
    return 0;
}
