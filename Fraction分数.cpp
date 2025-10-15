#include<bits/stdc++.h>
int gcd(int a,int b){if(b==0) return a;return gcd(b,a%b);}
int lcm(int a,int b){return a/gcd(a,b)*b;}

class Fraction {
public:
    int a,b;
    int sign(int x) {return (x>0?1:-1);}
    Fraction():a(0),b(1){}
    Fraction(int x):a(x),b(1){}
    Fraction(int x,int y)
    {
        int m = gcd(abs(x),abs(y));
        a = x/m*sign(y);
        if(a==0)b=1;else b = abs(y/m);
    }
    int get_denominator() {return b;}
    int get_numerator() {return a;}
    Fraction operator+(const Fraction &f)
    {
        int m = gcd(b,f.b);
        return Fraction(f.b/m*a+b/m*f.a,b/m*f.b);
    }
    Fraction operator-(const Fraction &f)
    {
        int m = gcd(b,f.b);
        return Fraction(f.b/m*a-b/m*f.a,b/m*f.b);
    }
    bool operator<(const Fraction &f)
    {
        int m = gcd(b,f.b);
        return f.b/m*a-b/m*f.a < 0;
    }
    bool operator<=(const Fraction &f)
    {
        int m = gcd(b,f.b);
        return f.b/m*a-b/m*f.a <= 0;
    }
    bool operator>(const Fraction &f)
    {
        int m = gcd(b,f.b);
        return f.b/m*a-b/m*f.a > 0;
    }
    bool operator>=(const Fraction &f)
    {
        int m = gcd(b,f.b);
        return f.b/m*a-b/m*f.a >= 0;
    }
    Fraction operator*(const Fraction &f)
    {
        int m1 = gcd(abs(a),f.b);
        int m2 = gcd(b,abs(f.a));
        return Fraction((a/m1)*(f.a/m2),(b/m2)*(f.b/m1));
    }
    Fraction operator/(const Fraction &f)
        {return (*this)*Fraction(f.b,f.a);}    
    friend ostream &operator << (ostream &out,const Fraction &f)
    {
        if(f.a==0) cout << 0;
        else if(f.b==1) cout << f.a;
        else cout << f.a << '/' << f.b;
        return out; 
    }
};