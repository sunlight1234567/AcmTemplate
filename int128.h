/*
本文件为128位整数的模板
128位整数的范围大概为 e32
c++标准库对于__int128并没有配套的函数操作，仅支持+、-、*、/ 四则运算
对于输入、输出、开方等操作需要自己实现
*/
#include <bits/stdc++.h>
using namespace std;
__int128 mod=1e30;
__int128 powmod(__int128 a,__int128 b) {__int128 res=1;a%=mod; assert(b>=0); for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}
__int128 gcd(__int128 a,__int128 b) { return b?gcd(b,a%b):a;}

__int128 str_i128(const string& s) {
    __int128 res=0;
    for(auto si:s){
        res = 10*res + si-'0';
    }
    return res;
}

string i128_str(__int128 x) {
    string res="";
    bool sign =(x<0? 1 : 0);
    if(sign)x *= -1;
    while(x){
        res.append(1, '0'+x%10);//?
        x/=10;
    }
    if(sign)res.append(1, '1');
    reverse(res.begin(),res.end());
    return res;
}
__int128 read128() {
    string tem;
    cin>>tem;
    return str_i128(tem);
}
void put(__int128 x) {
    stack<int> sta;
    char sign = x>=0? '+':'-';
    while(x){
        sta.push(x%10);
        x/=10;
    }
    if(sta.size()==0){cout<<'0';return;}
    if(sign=='-')cout<<sign;
    while(sta.size()){
        cout<<sta.top();
        sta.pop();
    }
}
/*
void put(__int128 x) {
    cout<<i128_str(x);
}
*/
__int128 Sqrt(const __int128& x) {
    //x>=0
    __int128 L=0, R=1e19;
    while(L + 1 < R){
        __int128 M = (L + R)/2;
        if(M * M <=x){
            L = M;
        }
        else R = M;
    }
    return L;
}