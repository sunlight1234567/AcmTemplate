#include<bits/stdc++.h>
using ll = long long;
using namespace std;

class segTree{
    //下标从1开始

    int N;//size of original array
    vector<ll> a;
    vector<ll> tree;//N<<2
    vector<ll> tag;//用于记录区间的修改量
public:
    segTree(int _N = 1e5){
        N = _N + 10;
        a = vector<ll>(N, 0);
        tree = vector<ll>(N*4, 0);
        tag = vector<ll>(N*4, 0);
    }
    segTree(vector<ll> _a) {
        a = _a;
        N = _a.size()+10;
        tree = vector<ll>(N*4, 0);
        tag = vector<ll>(N*4, 0);
        build_tree(1, 1, N);
    }
    void build_tree(int f, int l, int r)//f结点表示[l,r]的线段
    {
        if(l==r){tag[f]=0;tree[f]=a[l];return;}//叶子节点
        int mid=(l+r)/2;
        int lc=f*2, rc=f*2+1;
        build_tree(lc, l, mid);
        build_tree(rc, mid+1, r);
        tree[f]=tree[lc] + tree[rc];//push_up
    }

    void push_up(int f)
    {
        int lc=2*f, rc=2*f+1;
        tree[f]=tree[lc]+tree[rc];
    }

    void push_down(int f, int l, int r)
    {
        //f代表的[l,r]区间都加上k
        //将父亲结点的信息tag[]传给孩子节点，并更新孩子节点tree[]
        int mid=(l+r)/2;
        int lc=2*f, rc=2*f+1;
        tag[lc]+=tag[f];
        //tree[lc]+=(mid-l+1)*tag[lc];由于之前传递来的tag[]已经增加过了，现在只需增加新增的量tag[f]即可
        tree[lc]+=(mid-l+1)*tag[f];
        tag[rc]+=tag[f];
        //tree[rc]+=(r-mid)*tag[rc];
        tree[rc]+=(r-mid)*tag[f];
        tag[f]=0;//增量已经传递过了
    }


    void update(int _l, int _r, ll k, int f, int l, int r)//[_l, _r]为修改区间
    {
        if(l>_r || r<_l)return;
        if(l>=_l && r<=_r){
            tree[f]+=k*(r-l+1);
            tag[f]+=k;
            return;//整个区间直接修改，不需要向下传递
        }
        push_down(f, l, r);

        int mid=(l+r)/2;
        int lc=f*2, rc=f*2+1;
        update(_l, _r, k, lc, l, mid);
        update(_l, _r, k, rc, mid+1, r);

        push_up(f);
    }

    ll quiry(int _l, int _r, int f, int l, int r)//查询[_l, _r]
    {
        ll ans=0;
        if(l>=_l && r<=_r)return tree[f];
        if(l>_r || r<_l)return 0;
        int mid=(l+r)/2;
        int lc=f*2, rc=f*2+1;
        push_down(f, l, r);//将当前结点的信息tag[]传给孩子
        ans+=quiry(_l, _r, lc, l, mid);
        ans+=quiry(_l, _r, rc, mid+1, r);
        return ans;
    }
public:
    void add(int _l, int _r, ll k) {
        update(_l, _r, k, 1, 1, N);
    }
    ll sum(int _l, int _r){
        return quiry(_l, _r, 1, 1, N);
    }
};

int main()
{
    //freopen("segtree2.in", "r", stdin);
    //freopen("segtree.ans", "w", stdout);
    
    return 0;
}