// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_v.hpp
    \brief rank_support_v.hpp contains rank_support_v.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_INT_V
#define INCLUDED_SDSL_RANK_SUPPORT_INT_V

#include "rank_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl {


//! A rank structure proposed by Sebastiano Vigna
/*! \par Space complexity
 *  \f$ 0.25n\f$ for a bit vector of length n bits.
 *
 * The superblock size is 512. Each superblock is subdivided into 512/64 = 8
 * blocks. So absolute counts for the superblock add 64/512 bits on top of each
 * supported bit. Since the first of the 8 relative count values is 0, we can
 * fit the remaining 7 (each of width log(512)=9) in a 64bit word. The relative
 * counts add another 64/512 bits on top of each supported bit.
 * In total this results in 128/512=25% overhead.
 *
 * \tparam t_b       Bit pattern `0`,`1`,`10`,`01` which should be ranked.
 * \tparam t_pat_len Length of the bit pattern.
 *
 * \par Reference
 *    Sebastiano Vigna:
 *    Broadword Implementation of Rank/Select Queries.
 *    WEA 2008: 154-168
 *
 * @ingroup rank_support_group
 */
class rank_support_int_v : public rank_support_int {
public:
	typedef int_vector<> int_vector_type;
	typedef typename rank_support_int::size_type size_type;
	typedef typename rank_support_int::value_type value_type;
	// enum { bit_pat = t_b };
	// enum { bit_pat_len = t_pat_len };
	// static_assert(t_b == 2, "Currently only for t_b = 2 implemented");
private:
	// using rank_support_int::m_v;
	// basic block for interleaved storage of superblockrank and blockrank
	int_vector<64> m_basic_block;
	// constexpr static uint8_t t_v = 1ULL << t_b;

public:

	// TODO: set_vector needed by util::init_support!!!

	void printWord(uint64_t x) const
	{
		std::bitset<64> b(x);
    	for (signed i = 63; i >= 0; --i) {
			std::cout << b[i];
			if (i > 0 && i % 8 == 0)
				std::cout << ".";
		}
	}

	explicit rank_support_int_v(const int_vector<>* v = nullptr, unsigned max_val = 3) : rank_support_int(v, max_val)
	{
		// m_v = v;
		if (v == nullptr) {
			return;
		} else if (v->empty()) {
			m_basic_block = int_vector<64>(2 * (t_v - 1), 0); // resize structure for basic_blocks
			return;
		}
		// 8 * 64 bit = 512 bit abdeckung je superblock
		// 8 (eigentlich 7) blöcke à 8 bits in 1 64 wort. superblock 64 bit. interleaved
		// format:
		//             64 bit       64 bit         64 bit       64 bit         64 bit       64 bit
		// ..... |---------------------------|---------------------------|---------------------------| .....
		// ..... | superblock a | 8 blöcke a | superblock b | 8 blöcke b | superblock c | 8 blöcke c | .....
		// ..... |---------------------------|---------------------------|---------------------------| .....

		// (t_v - 1) * (64 + 64) + 8 * 64 = 3 * 128 + 512

		size_type basic_block_size = ((v->capacity() >> 9) + 1) * 2 * (t_v - 1);
		m_basic_block.resize(basic_block_size); // resize structure for basic_blocks
		const uint64_t* data = m_v->data();

		size_type i = 1, j = 0;
		uint64_t b_cnt[t_v - 1]= {0};
		uint64_t b_cnt_word[t_v - 1] = {0};

		for (value_type v = 0; v < t_v - 1; ++v) {
			m_basic_block[2*v] = m_basic_block[2*v + 1] = 0;
			b_cnt[v] = full_word_prefix_rank(data, 0, v);
		}

		for (; i < (m_v->capacity() >> 6) + 1; ++i) {
			for (value_type v = 0; v < t_v - 1; ++v) { // TODO: maybe switch loop over v and if statement over i % 8 (order)
				if (likely(i & 0x7)) { // if i%8!=0
					b_cnt_word[v] |= b_cnt[v] << ((i & 0x7) << 3); // 48 40 32 24 16 8 0
				} else {
					m_basic_block[j + 1] = b_cnt_word[v];
					m_basic_block[j + 2 * (t_v - 1)] = m_basic_block[j] + b_cnt[v];
					b_cnt[v] = b_cnt_word[v] = 0;
					j += 2;
				}
				b_cnt[v] += full_word_prefix_rank(data, i, v);
			}
		}

		for (value_type v = 0; v < t_v - 1; ++v, j += 2) {
			m_basic_block[j + 1] = b_cnt_word[v];
		}

		// for (unsigned x = 0; x < m_basic_block.size(); ++x)
		// {
		// 	if (x % 6 == 0)
		// 		std::cout << "---------------------------------------------------------------------\n";
		// 	if (x % 2 == 0)
		// 		std::cout << "sb ";
		// 	else
		// 		std::cout << "b  ";
		//
		// 	if (x % 6 == 0 || x % 6 == 1)
		// 		std::cout << "A  ";
		// 	else if (x % 6 == 2 || x % 6 == 3)
		// 		std::cout << "C  ";
		// 	else
		// 		std::cout << "G  ";
		//
		// 	printWord(m_basic_block[x]);
		// 	std::cout /*<< std::bitset<64>(m_basic_block[x])*/ << std::endl;
		// }
	}

	rank_support_int_v(const rank_support_int_v&) = default;
	rank_support_int_v(rank_support_int_v&&)	  = default;
	rank_support_int_v& operator=(const rank_support_int_v&) = default;
	rank_support_int_v& operator=(rank_support_int_v&&) = default;

	size_type rank(size_type idx, const value_type v) const
	{
		assert(m_v != nullptr);
		assert(idx <= m_v->size());

		if (unlikely(v == 0))
			return prefix_rank(idx, v);

		uint64_t word_pos = (2 * (t_v - 1) * (((idx * t_b) >> 9))) + 2 * v;
		const uint64_t* p = m_basic_block.data() + word_pos; // 2*(idx*t_b/512) + 2*v

		size_type result = 0;
		if (unlikely(v == t_v - 1)) // TODO: test effect of likely/unlikely
			result = idx;
		else
		 	result = *p + ((*(p + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF);
		result -= *(p - 2) + ((*(p - 2 + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF);

		if (likely(idx & 0x1F)) { // TODO: ein word_rank!
			if (likely(v != t_v - 1)) // if (idx % 32 != 0) nur für DNA-alphabet
				result += word_rank(m_v->data(), idx, v);
			else
				result -= word_prefix_rank(m_v->data(), idx, v - 1);
		}
		return result;
	}

	inline size_type operator()(size_type idx, const value_type v) const { return rank(idx, v); }

	size_type prefix_rank(size_type idx, const value_type v) const
	{
		assert(m_v != nullptr);
		assert(idx <= m_v->size());

		if (unlikely(v == t_v - 1)) // TODO: test effect of likely/unlikely
			return idx;

		uint64_t word_pos = (2 * (t_v - 1) * (((idx * t_b) >> 9))) + 2 * v;
		const uint64_t* p = m_basic_block.data() + word_pos; // 2*(idx*t_b/512) + 2*v

		// if (idx == 512 && v == 0) {
		// 	std::cout << "rank(" << idx << ", " << (unsigned) v
		// 			  << ") - sb: " << *p
		// 			  //<< " shift by " << ((8 * (((idx * t_b) & 0x1FF) / 64))) << " ... "
		// 			  << ", b : " << ((*(p + 1) >> (8 * (((idx * t_b) & 0x1FF) / 64))) & 0b11111111);
		// 	if (idx & (0x1F))
		// 		std::cout << ", pc: " << trait_type::word_rank(m_v->data(), idx, v);
		// 	std::cout << std::endl;
		// 	std::cout << "word_pos: " << word_pos << "\n";
		// }

		// TODO: test effect of likely/unlikely
		if (likely(idx & (0x1F))) // if (idx % 32 != 0) nur für DNA-alphabet
			return *p + ((*(p + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF) +
				   word_prefix_rank(m_v->data(), idx, v);
		else
			return *p + ((*(p + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF);
	}

	size_type size() const { return m_v->size(); }

	size_type
	serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		size_type			 written_bytes = 0;
		structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		written_bytes += m_basic_block.serialize(out, child, "cumulative_counts");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	void load(std::istream& in, const int_vector<>* v = nullptr)
	{
		m_v = v;
		m_basic_block.load(in);
		init(v, 0);
	}

	void set_vector(const int_vector<>* v = nullptr) { m_v = v; }
};

} // end namespace sds

#endif // end file
