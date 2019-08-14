#include <sdsl/vectors.hpp>
#include <sdsl/rank_support_int.hpp>
#include <sdsl/rank_support_int_scan.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

template <typename rank_t>
void test_int_vector(uint8_t const width, uint8_t const max_val, uint32_t const size)
{
    int_vector<0> iv;
    iv.width(width);
    //iv = {0,1,2,3,0,1,2,3}; // does not work because it sets the width to 64
    for (uint32_t i = 0; i < size; ++i)
        iv.push_back(i % 4);

    rank_t rb(&iv, max_val);
    for (unsigned v = 0; v < max_val; ++v)
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

    csa_wt<wt_blcd<>, 8, 8, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet> blcd;
    csa_wt<wt_pc_epr<>, 8, 8, sa_order_sa_sampling<>, isa_sampling<>, byte_alphabet> epr;

    auto const seed{1565824982/*time(NULL)*/};
    std::cout << "Seed: " << seed << '\n';
    srand(seed);

    while (true)
    {
        int_vector<8> text;//{1, 2, 3, 1, 2, 3, 1, 2, 3};
        int_vector<8> pattern;//{1, 2, 3};

        uint64_t text_len{std::rand() % 10000};
        uint64_t pattern_len{std::rand() % 15};
        text.resize(text_len);
        pattern.resize(pattern_len);

        for (uint32_t i = 0; i < text_len; ++i)
        {
            text[i] = (std::rand() % 3) + 1; // no 0s allowed. produces 1, 2 or 3.
            std::cout << (unsigned)text[i];
        }
        std::cout << std::endl;
        for (uint32_t i = 0; i < pattern_len; ++i)
        {
            pattern[i] = (std::rand() % 3) + 1; // no 0s allowed. produces 1, 2 or 3.
            std::cout << (unsigned)pattern[i];
        }
        std::cout << std::endl;

        construct_im(blcd, text, 0);
        construct_im(epr , text, 0);

        uint64_t lb{0}, rb{blcd.size()-1};      //, rev_lb{0}, rev_rb{epr.size()-1};
        // bidirectional_search(blcd, 'A', lb, rb, rev_lb, rev_rb);
        std::cout << lb << "," << rb << '\t';
        std::cout.flush();
        for (int32_t i = pattern.size() - 1; i >= 0; --i)
        {
            if (!backward_search(blcd, lb, rb, pattern[i], lb, rb))
                break;
            std::cout << lb << "," << rb << '\t';
            std::cout.flush();
        }
        std::cout << std::endl;

        lb = 0;
        rb = epr.size()-1;
        std::cout << lb << "," << rb << '\t';
        std::cout.flush();
        for (int32_t i = pattern.size() - 1; i >= 0; --i)
        {
            if (!backward_search(epr, lb, rb, pattern[i], lb, rb))
                break;
            std::cout << lb << "," << rb << '\t';
            std::cout.flush();
        }
        std::cout << std::endl << std::endl;

        // return 0;
    }

    return 0;

    std::cout << "Rank support v:\n";
    for (auto width : {2, 4, 8, 16})
    {
        for (uint64_t len = 0; len < 1000; ++len)
        {
            test_int_vector<rank_support_int_v<3,10>>(3, 6, len);

            test_int_vector<rank_support_int_v<1,2>>(width, 3, len);
            test_int_vector<rank_support_int_v<1,3>>(width, 3, len);
            test_int_vector<rank_support_int_v<1,4>>(width, 3, len);
            test_int_vector<rank_support_int_v<1,5>>(width, 3, len);

            test_int_vector<rank_support_int_v<2,2>>(width, 3, len);
            test_int_vector<rank_support_int_v<2,3>>(width, 3, len);
            test_int_vector<rank_support_int_v<2,4>>(width, 3, len);
            test_int_vector<rank_support_int_v<2,5>>(width, 3, len);

            test_int_vector<rank_support_int_v<3,2>>(width, 3, len);
            test_int_vector<rank_support_int_v<3,3>>(width, 3, len);
            test_int_vector<rank_support_int_v<3,4>>(width, 3, len);
            test_int_vector<rank_support_int_v<3,5>>(width, 3, len);
        }
    }

    // test_int_vector<rank_support_int_v<1,3>>(2, 100);
    // test_int_vector<rank_support_int_v<1,3>>(4, 100);
    //test_int_vector<rank_support_int_v>(5, 10000); // does not work as expected since values are overlapping 64bit-words (we don't want that ...)
}
