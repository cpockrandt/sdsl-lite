#include <sdsl/vectors.hpp>
#include <sdsl/rank_support_int.hpp>
#include <sdsl/rank_support_int_scan.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    int_vector<2> iv = {0,1,2,3,0,1,2,3};
    // bit_vector b = {0,1,1,1,1,0,0,0,0,0};
    // rank_support_v<1,1> rb(&b);
    // cout << rb(5) << endl;
    rank_support_int_scan<2,0> rb(&iv);
    for (unsigned i = 0; i < iv.size() + 1; ++i)
        cout << "rank(" << i << "):" << rb(i) << '\n';
}
