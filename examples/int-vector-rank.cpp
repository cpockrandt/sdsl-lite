#include <sdsl/vectors.hpp>
#include <sdsl/rank_support_int.hpp>
#include <sdsl/rank_support_int_scan.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main()
{
    // bit_vector bv = {1,1,1,1};
    // rank_support_v<> rb(&bv);
    // for (unsigned i = 0; i < bv.size() + 1; ++i)
    //     std::cout << rb(i) << std::endl;
    int_vector<2> iv = {0,1,2,3,0,1,2,3};
    rank_support_int_scan<2, 4> rb(&iv);
    rank_support_int_trait<2, 4>::init();
    cout << "IV: " << iv << endl << endl;

    for (unsigned i = 0; i < iv.size() + 1; ++i)
    {
        cout << i << "\t";
        for (unsigned v = 0; v < 4; ++v)
            cout << rb(i, v) << "\t";
        cout << endl;
    }

    for (unsigned i = 0; i < iv.size() + 1; ++i)
    {
        cout << i << "\t";
        for (unsigned v = 0; v < 4; ++v)
            cout << rb(i, v) << "\t";
        cout << endl;
    }

    for (unsigned i = 0; i < iv.size() + 1; ++i)
    {
        cout << i << "\t";
        for (unsigned v = 0; v < 4; ++v)
            cout << rb(i, v) << "\t";
        cout << endl;
    }
    // for (unsigned v = 0; v < 4; ++v)
    // {
    //     unsigned cnt = 0;
    //     for (unsigned i = 0; i < iv.size() + 1; ++i)
    //     {
    //         if (i > 0 && iv[i - 1] <= v)
    //             ++cnt;
    //         // uint64_t tmp234 = rb(i, v);
    //         // cout << (1 == rb(i, v)) << "";
    //         if (rb(i, v) != cnt)
    //         {
    //             std::cout << "ERROR!!!!!!!!!!!!!!!\n";
    //             std::cout << "i = " << i << ", v = " << v << ", exp = " << cnt << ", actual = " << rb(i, v) << "\n";
    //         }
    //     }
    // }
}
