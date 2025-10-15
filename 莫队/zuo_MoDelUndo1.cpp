// 只删回滚莫队入门题，C++版
// 本题最优解为主席树，讲解158，题目2，已经讲述
// 给定一个长度为n的数组arr，一共有m条查询，格式如下
// 查询 l r : 打印arr[l..r]内没有出现过的最小自然数，注意0是自然数
// 0 <= n、m、arr[i] <= 2 * 10^5
// 测试链接 : https://www.luogu.com.cn/problem/P4137
// 如下实现是C++的版本，C++版本和java版本逻辑完全一样
// 提交如下代码，可以通过所有测试用例

#include <bits/stdc++.h>

using namespace std;

struct Query {
    int l, r, id;
};

const int MAXN = 200001;
const int MAXB = 501;
int n, m;
vector<int> arr;
vector<Query> query;

int blen, bnum;
vector<int> bi;
vector<int> bl;

vector<int> cnt;
int mex;
vector<int> ans;

bool QueryCmp(const Query &a, const Query &b) {
    if (bi[a.l] != bi[b.l]) {
        return bi[a.l] < bi[b.l];
    }
    return b.r < a.r;
}

void del(int num) {
    if (--cnt[num] == 0) {
        mex = min(mex, num);
    }
}

void add(int num) {
    cnt[num]++;
}

void compute() {
    // 初始化计数数组
    for (int i = 1; i <= n; i++) {
        cnt[arr[i]]++;
    }
    
    // 计算初始mex
    mex = 0;
    while (cnt[mex] != 0) {
        mex++;
    }
    
    int winl = 1, winr = n;
    for (int block = 1, qi = 1; block <= bnum && qi <= m; block++) {
        // 移动左边界到当前块的起始位置
        while (winl < bl[block]) {
            del(arr[winl++]);
        }
        
        int beforeJob = mex;
        
        // 处理当前块内的所有查询
        for (; qi <= m && bi[query[qi].l] == block; qi++) {
            int jobl = query[qi].l;
            int jobr = query[qi].r;
            int id = query[qi].id;
            
            // 移动右边界到查询右边界
            while (winr > jobr) {
                del(arr[winr--]);
            }
            
            int backup = mex;
            
            // 移动左边界到查询左边界
            while (winl < jobl) {
                del(arr[winl++]);
            }
            
            ans[id] = mex;
            mex = backup;
            
            // 回滚左边界到当前块的起始位置
            while (winl > bl[block]) {
                add(arr[--winl]);
            }
        }
        
        // 恢复右边界到数组末尾
        while (winr < n) {
            add(arr[++winr]);
        }
        
        mex = beforeJob;
    }
}

void prepare() {
    blen = (int)sqrt(n);
    bnum = (n + blen - 1) / blen;
    
    // 初始化分块数组
    bi.resize(n + 1);
    bl.resize(bnum + 1);
    
    for (int i = 1; i <= n; i++) {
        bi[i] = (i - 1) / blen + 1;
    }
    for (int i = 1; i <= bnum; i++) {
        bl[i] = (i - 1) * blen + 1;
    }
    
    // 排序查询
    sort(query.begin() + 1, query.end(), QueryCmp);
}

int main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);

	cin >> n >> m;

	// 使用vector替代数组
	arr.resize(n + 1);
	query.resize(m + 1);
	cnt.resize(MAXN, 0);
	ans.resize(m + 1);

	for (int i = 1; i <= n; i++) {
		cin >> arr[i];
	}
	for (int i = 1; i <= m; i++) {
		cin >> query[i].l >> query[i].r;
		query[i].id = i;
	}

	prepare();
	compute();

	for (int i = 1; i <= m; i++) {
		cout << ans[i] << '\n';
	}
	return 0;
}
