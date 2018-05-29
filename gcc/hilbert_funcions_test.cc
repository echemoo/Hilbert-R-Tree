#include <iostream>
using namespace std;

int main()
{
  int rotation = 1;
  int nDims = 2;
  int bits = 0;
  int nd1Ones = 1;
 
  cout << rotation << "\t" << nDims << "\t" << bits << "\t" << nd1Ones << "\r\n";


  do {
    bits &= -bits & nd1Ones;
    while (bits)
      bits >> 1, ++rotation;
    if ( ++rotation >= nDims)
      rotation -= nDims;
  } while(0);

  cout << rotation << "\t" << nDims << "\t" << bits << "\t" << nd1Ones << "\r\n";


  return 0;
}
