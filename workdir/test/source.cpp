/*
两个问题：
数组开小了，导致进位时越界；
sum数组处理不好，前一个例子较长的话sum数组未清理导致影响下一个较短的例子
*/
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
using Int = long long;
const double PI = acos(-1.0);

struct Complex {
  double x, y;

  Complex(double _x = 0.0, double _y = 0.0) {
    x = _x;
    y = _y;
  }

  Complex operator-(const Complex &b) const {
    return Complex(x - b.x, y - b.y);
  }

  Complex operator+(const Complex &b) const {
    return Complex(x + b.x, y + b.y);
  }

  Complex operator*(const Complex &b) const {
    return Complex(x * b.x - y * b.y, x * b.y + y * b.x);
  }
};

/*
 * 进行 FFT 和 IFFT 前的反置变换
 * 位置 i 和 i 的二进制反转后的位置互换
 *len 必须为 2 的幂
 */
void change(Complex y[], int len) {
  int i, j, k;

  for (int i = 1, j = len / 2; i < len - 1; i++) {
    if (i < j) std::swap(y[i], y[j]);

    // 交换互为小标反转的元素，i<j 保证交换一次
    // i 做正常的 + 1，j 做反转类型的 + 1，始终保持 i 和 j 是反转的
    k = len / 2;

    while (j >= k) {
      j = j - k;
      k = k / 2;
    }

    if (j < k) j += k;
  }
}

/*
 * 做 FFT
 *len 必须是 2^k 形式
 *on == 1 时是 DFT，on == -1 时是 IDFT
 */
/*
 * 做 FFT
 * len 必须是 2^k 形式
 * on == 1 时是 DFT，on == -1 时是 IDFT
 */
void fft(Complex y[], int len, int on) {
  // 位逆序置换
  change(y, len);
  // 模拟合并过程，一开始，从长度为一合并到长度为二，一直合并到长度为 len。
  for (int h = 2; h <= len; h <<= 1) {
    // wn：当前单位复根的间隔：w^1_h
    Complex wn(cos(2 * PI / h), sin(on * 2 * PI / h));
    // 合并，共 len / h 次。
    for (int j = 0; j < len; j += h) {
      // 计算当前单位复根，一开始是 1 = w^0_n，之后是以 wn 为间隔递增： w^1_n
      // ...
      Complex w(1, 0);
      for (int k = j; k < j + h / 2; k++) {
        // 左侧部分和右侧是子问题的解
        Complex u = y[k];
        Complex t = w * y[k + h / 2];
        // 这就是把两部分分治的结果加起来
        y[k] = u + t;
        y[k + h / 2] = u - t;
        // 后半个 「step」 中的ω一定和 「前半个」 中的成相反数
        // 「红圈」上的点转一整圈「转回来」，转半圈正好转成相反数
        // 一个数相反数的平方与这个数自身的平方相等
        w = w * wn;
      }
    }
  }
  // 如果是 IDFT，它的逆矩阵的每一个元素不只是原元素取倒数，还要除以长度 len。
  if (on == -1) {
    for (int i = 0; i < len; i++) {
      y[i].x /= len;
      y[i].y /= len;
    }
  }
}

constexpr int MAXN = 800020;
Complex x1[MAXN], x2[MAXN];
std::string str1, str2;
using std::cin;
using std::cout;
using std::vector;
signed main() {
  cin.tie(nullptr)->sync_with_stdio(false);
  int t;cin >> t;
  while (t--) {
    cin >> str1 >> str2;
    int len1 = str1.length();
    int len2 = str2.length();
    int len = 1;

    while (len < len1 * 2 || len < len2 * 2) len <<= 1;

    for (int i = 0; i < len1; i++) x1[i] = Complex(str1[len1 - 1 - i] - '0', 0);

    for (int i = len1; i < len; i++) x1[i] = Complex(0, 0);

    for (int i = 0; i < len2; i++) x2[i] = Complex(str2[len2 - 1 - i] - '0', 0);

    for (int i = len2; i < len; i++) x2[i] = Complex(0, 0);

    fft(x1, len, 1);
    fft(x2, len, 1);

    for (int i = 0; i < len; i++) x1[i] = x1[i] * x2[i];

    fft(x1, len, -1);
    len = len1+len2-1;
    vector<int> sum(len);
    for (int i = 0; i < len; i++) sum[i] = (long long)(x1[i].x + 0.5);

    // cout << "sum:\n";
    // for(int i=0; i<=len; ++i)cout << sum[i] << ' ';
    // cout << "sum end\n";
    for (int i = 0; i <= sum.size(); i++) {
      if(sum[i]==0 || sum[i]==1) continue;
        else{
            if(i+2>=sum.size()) sum.resize(i+2+1, 0);
            sum[i+2]-=sum[i]/2;
            sum[i] = sum[i]%2;
            if(sum[i] == -1) {
                sum[i+2]+=1;
                sum[i]=1;
            }
        }
    }  
//   for (int i = 0; i <= len; ++i) {
//       const Int r = sum[i] & 1;
//       const Int q = (sum[i] - r) / (-2);
//       if (q) {
//         len = std::max(len, i+2);
//         sum[i + 2] += q;
//         sum[i] = r;
//       }
//     }
    while (sum.size()>1 && sum.back()==0) sum.pop_back();

    for (int i = sum.size()-1; i >= 0; i--) cout << char(sum[i] + '0');

    cout << '\n';
  }
    
  return 0;
}