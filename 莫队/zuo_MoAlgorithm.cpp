// 普通莫队入门题，C++版
// 给定一个长度为n的数组arr，一共有q条查询，格式如下
// 查询 l r : 打印arr[l..r]范围上有几种不同的数字
// 1 <= n <= 3 * 10^4
// 1 <= arr[i] <= 10^6
// 1 <= q <= 2 * 10^5
// 测试链接 : https://www.luogu.com.cn/problem/SP3267
// 测试链接 : https://www.spoj.com/problems/DQUERY/
// 如下实现是C++的版本，C++版本和java版本逻辑完全一样
// 提交如下代码，可以通过所有测试用例

#include <bits/stdc++.h>

using namespace std;

struct Query {
    int l, r, id;
};

const int MAXV = 1000001;
int n, q;
vector<int> arr;
vector<Query> query;

vector<int> bi;
vector<int> cnt;
int kind = 0;

vector<int> ans;

bool QueryCmp1(const Query &a, const Query &b) {
    if (bi[a.l] != bi[b.l]) {
        return bi[a.l] < bi[b.l];
    }
    return a.r < b.r;
}

bool QueryCmp2(const Query &a, const Query &b) {
    if (bi[a.l] != bi[b.l]) {
        return bi[a.l] < bi[b.l];
    }
    if ((bi[a.l] & 1) == 1) {
        return a.r < b.r;
    } else {
        return a.r > b.r;
    }
}

void del(int num) {
    if (--cnt[num] == 0) {
        kind--;
    }
}

void add(int num) {
    if (++cnt[num] == 1) {
        kind++;
    }
}

void prepare() {
    int blen = (int)sqrt(n);
    bi.resize(n + 1);
    for (int i = 1; i <= n; i++) {
        bi[i] = (i - 1) / blen + 1;
    }
    sort(query.begin() + 1, query.end(), QueryCmp2);
}

void compute() {
    int winl = 1, winr = 0;
    for (int i = 1; i <= q; i++) {
        int jobl = query[i].l;
        int jobr = query[i].r;
        while (winl > jobl) {
            add(arr[--winl]);
        }
        while (winr < jobr) {
            add(arr[++winr]);
        }
        while (winl < jobl) {
            del(arr[winl++]);
        }
        while (winr > jobr) {
            del(arr[winr--]);
        }
        ans[query[i].id] = kind;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n;
    arr.resize(n + 1);
    for (int i = 1; i <= n; i++) {
        cin >> arr[i];
    }
    
    cin >> q;
    query.resize(q + 1);
    ans.resize(q + 1);
    for (int i = 1; i <= q; i++) {
        cin >> query[i].l >> query[i].r;
        query[i].id = i;
    }
    
    cnt.resize(MAXV, 0);
    prepare();
    compute();
    
    for (int i = 1; i <= q; i++) {
        cout << ans[i] << '\n';
    }
    return 0;
}
