#include <iostream>
using namespace std;

typedef unsigned long long bitmask_t;
typedef unsigned long halfmask_t;

int main()
{
    #define maxDim (8*sizeof(bitmask_t))
    cout << "Print Var maxDim:" << maxDim << endl;
    bitmask_t coord[maxDim], coordPrev[maxDim];
    // nDims:  维度
    // nBits:  没维的字节数量
    // nPrints: 需要打印数量
    // orderCheck: TODO need to repair
    unsigned nDims, nBits, nPrints, orderCheck, i;
    bitmask_t r, r1;

    for (;;)
    {
    }

    return 0;
}
