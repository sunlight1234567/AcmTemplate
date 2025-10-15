using mint = modnum<mod>;
struct Combi{
    vector<mint> _fac, _ifac;
    int n;
    
    Combi() {
        n = 1;
        _fac.assign(n + 1, 1);
        _ifac.assign(n + 1, 1);
    }
    
    void check_size(int m){
        int need = n;
        while (need < m) need *= 2;
        m = need;
        if (m <= n) return;
        
        _fac.resize(m + 1);
        _ifac.resize(m + 1);
        for (int i = n + 1; i <= m; i++) _fac[i] = i * _fac[i - 1];
        
        _ifac[m] = _fac[m].inv();
        for (int i = m - 1; i > n; i--) _ifac[i] = _ifac[i + 1] * (i + 1);
        n = m;
    }
    
    mint fac(int m){
        check_size(m);
        return _fac[m];
    }
    
    mint ifac(int m){
        check_size(m);
        return _ifac[m];
    }
    
    mint ncr(int n, int r){//comb n选r
        if (n < r || r < 0) return 0;
        
        return fac(n) * ifac(n - r) * ifac(r);
    }
    
    mint npr(int n, int r){
        if (n < r || r < 0) return 0;
        
        return fac(n) * ifac(n - r);
    }
} comb;