#include <iostream>
#include "Integer.h" //自己实现的类直接导入.h文件就好啦 代码再下面
using namespace std;

Integer A[100010];

int fac(int n)
{
	// 给每个对象赋值
	for (int i = 1; i <= n; i++)
	{
		A[i] = i;
	}
	// 用分治的方法来计算(可以多线程)
	// 每次乘法为O(nlogn) 然后分治为O(logn) 总的复杂度为O(nlogn*logn) 十万以内几秒吧
	for (int i = 1; i <= n; i *= 2)
	{
		for (int j = 1; j <= n; j += 2 * i)
		{
			if (j + i <= n)
				A[j] = A[j] * A[j + i];
		}
	}
	// 结果
	cout << A[1] << endl;
	// cout << A[1].size() << endl;
	return 0;
}

int main()
{

	int n;
	// cin >> n;
	n = 100; // 懒得一次次输入啦

	fac(n);

	return 0;
}
