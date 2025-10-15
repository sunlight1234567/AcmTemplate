// manacher
// 扩展字符串
// 考虑当前位置i与回文右边界r的位置；考虑i关于回文中心c的对称点j的回文区间对i的回文区间的贡献
// 对每一个位置求最长回文半径p[i]，原字符串回文长度len = p[i] - 1;
// 扩展串回文子串的末位置 l2 对应在原字符串的末位置 l1 : l1 = l2/2 - 1
#include<bits/stdc++.h>
using namespace std;

const int N = 1e7 + 1e6 + 10;
void sol() {
    string s;
    cin >> s;
    int n = s.length();
    vector<char> ss(2*n + 1);
    vector<int> p(2*n + 1);
    for(int i=0; i<2*n+1; ++i) {
        p[i] = 1;
        ss[i] = (i & 1 ? s[i/2] : '#');
    }
    //for(auto x : ss) cout << x;
    //cout << endl;
    n = 2*n + 1;
    int r = 0, c = -1;
    for(int i=0; i<n; ++i) {
        int j = 2*c - i;
        int len = (i < r? min(p[j], r-i) : 1);
        while(i + len < n && i-len >=0 && ss[i+len] == ss[i-len]) ++len;
        p[i] = len;
        if( i + len > r) {
            r = i + len;
            c = i;
        }
    }
    int ans = 0;
    for(int i=0; i<n; ++i) ans= max(ans, p[i]);
    --ans;
    cout << ans;
}
signed main() {
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int T = 1;
    //cin >> T;
    while(T--) sol();
    return 0;
}

std::vector<int> manacher(std::string s) {
    std::string t = "#";
    for (auto it : s) {
        t += it;
        t += '#';
    }
    int n = t.size();
    std::vector<int> r(n);
    for (int i = 0, j = 0; i < n; i++) {
        if (2 * j - i >= 0 && j + r[j] > i) {
            r[i] = std::min(r[2 * j - i], j + r[j] - i);
        }
        while (i - r[i] >= 0 && i + r[i] < n && t[i - r[i]] == t[i + r[i]]) {
            r[i] += 1;
        }
        if (i + r[i] > j + r[j]) {
            j = i;
        }
    }
    return r;
}