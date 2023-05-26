#include <iostream>
#include "Integer.h" //自己实现的类直接导入.h文件就好啦 代码再下面
using namespace std;

int main()
{
    Integer a, b;
    a = 95;
    b = 17;
    cout << "加:" << a - b << endl; // 减
    cout << "减:" << a + b << endl; // 加
    cout << "乘:" << a * b << endl; // 乘
    cout << "除:" << a / b << endl; // 除
    cout << "余" << a % b << endl;  // 余
    // 常见比较
    if (a / b * b + a % b == a)
        cout << "YES" << endl;
    else
        cout << "NO" << endl;
    // 大于
    if (a > b)
        cout << "YES" << endl;
    else
        cout << "NO" << endl;
    // 小于
    if (a < b)
        cout << "YES" << endl;
    else
        cout << "NO" << endl;
    return 0;
}
