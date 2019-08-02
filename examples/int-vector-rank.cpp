#include <sdsl/vectors.hpp>
#include <sdsl/rank_support_int.hpp>
#include <sdsl/rank_support_int_scan.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

template <typename rank_t>
void test_int_vector(uint8_t const width, uint32_t const size)
{
    uint32_t const alph_size = 1 << width;

    int_vector<0> iv;
    iv.width(width);
    //iv = {0,1,2,3,0,1,2,3}; // does not work because it sets the width to 64
    for (uint32_t i = 0; i < size; ++i)
        iv.push_back(i % 4);

    rank_t rb(&iv);
    for (unsigned v = 0; v < alph_size; ++v)
    {
        unsigned cnt = 0;
        for (unsigned i = 0; i < iv.size() + 1; ++i)
        {
            if (i > 0 && iv[i - 1] == v)
                ++cnt;

            // auto const rank = rb(i, v);
            // if (rank != cnt)
            // {
            //     std::cout << iv << endl;
            //     std::cout << "WIDTH: " << (unsigned)width << endl;
            //     std::cout << "RANK ERROR!!!!!!!!!!!!!!!\n";
            //     std::cout << "i = " << i << ", v = " << v << ", exp = " << cnt << ", actual = " << rank << "\n";
            //     exit(1);
            // }

            auto const prefix_rank = rb.prefix_rank(i, v);
            auto const exp_prefix_rank = cnt + (v == 0 ? 0 : rb.prefix_rank(i, v - 1));
            if (prefix_rank != exp_prefix_rank)
            {
                std::cout << iv << endl;
                std::cout << "WIDTH: " << (unsigned)width << endl;
                std::cout << "PREFIX RANK ERROR!!!!!!!!!!!!!!!\n";
                std::cout << "i = " << i << ", v = " << v << ", exp = " << exp_prefix_rank << ", actual = " << prefix_rank << "\n";
                exit(2);
            }
        }
    }
}

int main()
{
    // TODO: test with random data, short data, edge cases (multiple of 32 chars, =/- 1 for DNA)

    // std::cout << "Rank support scan:\n";
    // test_int_vector<rank_support_int_scan>(2, 10000);
    // test_int_vector<rank_support_int_scan>(4, 10000);
    //test_int_vector<rank_support_int_scan>(5, 10000); // does not work as expected since values are overlapping 64bit-words (we don't want that ...)

    std::cout << "Rank support v:\n";
    test_int_vector<rank_support_int_v<1,3>>(16, 4*9);
    // test_int_vector<rank_support_int_v<1,3>>(2, 100);
    // test_int_vector<rank_support_int_v<1,3>>(4, 100);
    //test_int_vector<rank_support_int_v>(5, 10000); // does not work as expected since values are overlapping 64bit-words (we don't want that ...)
}
