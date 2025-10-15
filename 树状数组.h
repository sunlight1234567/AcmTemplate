/*
树形数组
支持单点修改，单点查询；单点修改，区间查询；区间修改，单点查询（不能区间查询）
进行区间修改时，只需将数组差分，用树状数组维护差分数组即可；
https://blog.csdn.net/TheWayForDream/article/details/118436732
*/
#include <bits/stdc++.h>
using namespace std;

int lowbit(int x) { return x & -x; }
const int N = 500;
struct TreeArr {
    int _size;
    vector<int> tree;
    TreeArr(){
        _size = N;
        tree.resize(N + 10, 0);
    }
    TreeArr(int size) {
        _size = size;
        tree.resize(_size + 10, 0);
    }
    void add(int i, int d) {
        while(i <= _size) {
            tree[i] += d;
            i += (i & -i);
        }
    }
    int sum(int i) {
        int res = 0;
        while(i > 0) {
            res += tree[i];
            i -= (i & -i);
        }
        return res;
    }
    int query(int l, int r) {
        return sum(r) - sum(l-1);
    }
};

struct TreeArrMatrix {
    int _size;
    vector<vector<int>> tree;
    TreeArrMatrix(){
        _size = N;
        tree.resize(N + 10, vector<int>(N + 10, 0));
    }
    TreeArrMatrix(int size) {
        _size = size;
        tree.resize(_size + 10, vector<int>(_size + 10, 0));
    }
    void resize(int size){
        _size = size;
        tree.resize(_size + 10, vector<int>(_size + 10, 0));
    }
    void add(int x, int y, int d) {
        for(int i = x; i <= _size; i += lowbit(i)) {
            for(int j=y; j<=_size; j += lowbit(j)) {
                tree[i][j] += d;
            }
        }
    }
    int sum(int x, int y) {
        int res = 0;
        for(int i=x; i>0; i-= lowbit(i)) {
            for(int j=y; j>0; j -= lowbit(j)) {
                res += tree[i][j];
            }
        }
        return res;
    }
    int query(int xl, int yl, int xr, int yr) {
        return sum(xr, yr) - sum(xr, yl-1) - sum(xl-1, yr) + sum(xl-1, yl-1);
    }
};