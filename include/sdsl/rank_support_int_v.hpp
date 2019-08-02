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
// NOTE: a word is 64 bit
template <uint8_t words_per_block = 1, uint8_t blocks_per_superblock = 4/*, uint8_t alphabet_size = 0*/>
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
	int_vector<0> m_block;
	int_vector<64> m_superblock; // TODO: set width at runtime? will probably affect speed
	// constexpr static uint8_t t_v = 1ULL << t_b;

	int_vector<64> m_basic_block;

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

	explicit rank_support_int_v(const int_vector<>* v = nullptr, unsigned max_val = 0) : rank_support_int(v, max_val)
	{
		uint8_t const t_v_decr = t_v - 1; // max_val not needed because of t_v?

		if (v == nullptr) {
			return;
		} else if (v->empty()) {
			m_block.resize(t_v_decr, 0); // TODO: check when rank/prefix-rank is done!
			m_superblock.resize(t_v_decr, 0); // TODO: check when rank/prefix-rank is done!
			// m_basic_block = int_vector<64>(2 * (t_v - 1), 0); // resize structure for basic_blocks
			return;
		}

		// TODO: use __builtin_clz instead of log
		// uint32_t const bits_per_value{std::ceil(log2(max_val))};
		// uint32_t const bits_per_value{std::ceil(log2(max_val))};
		constexpr uint64_t words_per_superblock = words_per_block * blocks_per_superblock;

		uint64_t const values_per_word{64 / v->width()};
		uint64_t const values_per_block{words_per_block * values_per_word};
		uint64_t const values_per_superblock{blocks_per_superblock * values_per_block};
		uint64_t const new_width{std::ceil(log2(values_per_superblock))};
		m_block.width(new_width);
		// TODO: could also set block width of superblocks. check running time impact!
		std::cout << "Bits per value: " << (unsigned)v->width() << std::endl;
		std::cout << "Values per word: " << values_per_word << std::endl;
		std::cout << "Values per superblock: " << values_per_superblock << std::endl;
		std::cout << "Block width: " << (unsigned)new_width << std::endl;

		uint64_t const word_count = ((m_v->size() - 1) / values_per_word) + 1; // equivalent to ceil(m_v->size() / values_per_word)
		uint64_t const block_count = ((word_count - 1) / words_per_block) + 1; // equivalent to ceil(word_count / words_per_block)

		size_type const block_size = (((block_count - 1) / blocks_per_superblock) + 1) * (blocks_per_superblock - 1) * t_v_decr; // ceil(v->size() / values_per_block) * t_v_decr
		size_type const superblock_size = (((v->size() - 1) / values_per_superblock) + 1) * t_v_decr; // ceil(v->size() / values_per_superblock) * t_v_decr
		std::cout << "t_v_decr: " << (unsigned)t_v_decr << std::endl;
		std::cout << "Elements: " << v->size() << std::endl;
		std::cout << "block count: " << block_count << std::endl;
		std::cout << "word count: " << word_count << std::endl;
		std::cout << "block_size: " << block_size << std::endl;
		std::cout << "superblock_size: " << superblock_size << std::endl;
		m_block.resize(block_size); // TODO: adjust
		m_superblock.resize(superblock_size); // TODO: adjust

		uint64_t const * data = m_v->data();
		std::vector<uint64_t> buf_blocks(t_v_decr, 0); // TODO: get rid of these objects
		std::vector<uint64_t> buf_superblocks(t_v_decr, 0);

		for (uint64_t v = 0; v < t_v_decr; ++v)
			m_superblock[v] = 0;

		for (uint64_t word_id = 0, block_id = 0, superblock_id = t_v_decr; word_id < word_count; word_id += t_v_decr)
		{
			// TODO: for loop missing over multiple words per block
			for (uint64_t v = 0; v < t_v_decr; ++v)
				buf_blocks[v] += full_word_prefix_rank(data, word_id, v);

			if (word_id % words_per_block) // divisor is constexpr, i.e., modulo operation is expected to be cheap
			{
				if (word_id % words_per_superblock != 0) // divisor is constexpr, i.e., modulo operation is expected to be cheap
				{ // don't store block information for the last block in the superblock!
					for (uint64_t v = 0; v < t_v_decr; ++v)
						m_block[block_id + v] = buf_blocks[v];
					block_id += t_v_decr;
				}
				else
				{
					for (uint64_t v = 0; v < t_v_decr; ++v)
					{
						buf_superblocks[v] += buf_blocks[v];
						m_superblock[superblock_id + v] = buf_superblocks[v];
						buf_blocks[v] = 0; // reset blocks
					}
					superblock_id += t_v_decr;
				}
			}
		}

        std::cout << "\nBlocks:\n";
        for (uint64_t i = 0; i < m_block.size(); ++i)
            std::cout << (unsigned)m_block[i] << ' ';
        std::cout << "\nSuperBlocks:\n";
        for (uint64_t i = 0; i < m_superblock.size(); ++i)
            std::cout << (unsigned)m_superblock[i] << ' ';
        std::cout << "\n\n";

		// size_type basic_block_size = ((v->capacity() >> 9) + 1) * 2 * (t_v - 1);
		// m_basic_block.resize(basic_block_size); // resize structure for basic_blocks
		// const uint64_t* data = m_v->data();
		//
		// size_type i = 1, j = 0;
		// uint64_t b_cnt[t_v - 1]= {0};
		// uint64_t b_cnt_word[t_v - 1] = {0};
		//
		// for (value_type v = 0; v < t_v - 1; ++v) {
		// 	m_basic_block[2*v] = m_basic_block[2*v + 1] = 0;
		// 	b_cnt[v] = full_word_prefix_rank(data, 0, v);
		// }
		//
		// for (; i < (m_v->capacity() >> 6) + 1; ++i) {
		// 	for (value_type v = 0; v < t_v - 1; ++v) { // TODO: maybe switch loop over v and if statement over i % 8 (order)
		// 		if (likely(i & 0x7)) { // if i%8!=0
		// 			b_cnt_word[v] |= b_cnt[v] << ((i & 0x7) << 3); // 48 40 32 24 16 8 0
		// 		} else {
		// 			m_basic_block[j + 1] = b_cnt_word[v];
		// 			m_basic_block[j + 2 * (t_v - 1)] = m_basic_block[j] + b_cnt[v];
		// 			b_cnt[v] = b_cnt_word[v] = 0;
		// 			j += 2;
		// 		}
		// 		b_cnt[v] += full_word_prefix_rank(data, i, v);
		// 	}
		// }
		//
		// for (value_type v = 0; v < t_v - 1; ++v, j += 2) {
		// 	m_basic_block[j + 1] = b_cnt_word[v];
		// }




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

		return 0;

		// uint64_t word_pos = (2 * (t_v - 1) * (((idx * t_b) >> 9))) + 2 * v;
		// const uint64_t* p = m_basic_block.data() + word_pos; // 2*(idx*t_b/512) + 2*v
		//
		// size_type result = 0;
		// if (unlikely(v == t_v - 1)) // TODO: test effect of likely/unlikely
		// 	result = idx;
		// else
		//  	result = *p + ((*(p + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF);
		// result -= *(p - 2) + ((*(p - 2 + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF);
		//
		// if (likely(idx & 0x1F)) { // TODO: ein word_rank!
		// 	if (likely(v != t_v - 1)) // if (idx % 32 != 0) nur für DNA-alphabet
		// 		result += word_rank(m_v->data(), idx, v);
		// 	else
		// 		result -= word_prefix_rank(m_v->data(), idx, v - 1);
		// }
		// return result;
	}

	inline size_type operator()(size_type idx, const value_type v) const { return rank(idx, v); }

	size_type prefix_rank(size_type idx, const value_type v) const
	{
		assert(m_v != nullptr);
		assert(idx <= m_v->size());

		if (unlikely(v == t_v - 1)) // TODO: test effect of likely/unlikely
			return idx;

		uint8_t const t_v_decr = t_v - 1;

		uint32_t const bits_per_value{m_v->width()};
		uint32_t const values_per_word{64 / bits_per_value};
		uint32_t const values_per_block{words_per_block * values_per_word};
		// uint64_t const values_per_superblock{blocks_per_superblock * values_per_block};

		// constexpr uint16_t words_per_superblock = words_per_block * blocks_per_superblock;

		// uint64_t word_pos = (2 * (t_v - 1) * (((idx * t_b) >> 9))) + 2 * v;
		// uint64_t* const p = m_basic_block.data() + word_pos; // 2*(idx*t_b/512) + 2*v

		// size_type const word_id = idx / values_per_word; // TODO: expensive!
		size_type const block_id = idx / values_per_block;
		size_type const superblock_id = block_id / blocks_per_superblock; // NOTE: former idx / values_per_superblock

		// superblock_id * (blocks_per_superblock - 1) + (word_id % words_per_block); // TODO: expensive!

		size_type res = m_superblock[t_v_decr * superblock_id + v];

		uint16_t const block_id_in_superblock = block_id % blocks_per_superblock;
		if (block_id_in_superblock != 0)
			res += m_block[t_v_decr * (blocks_per_superblock - 1) * block_id_in_superblock + v];

		// TODO: loop for multiple words in a block
		if (idx % values_per_block != 0)
			res += word_prefix_rank(m_v->data(), idx, v);

		return res;

		// uint64_t const superblock = m_superblock[idx / values_per_superblock + v];
		// uint64_t const block = m_block[idx / values_per_block + v];
		//
		// if (likely(idx & (0x1F))) // TODO: if (idx % 32 != 0) nur für DNA-alphabet
		// 	return superblock + block + word_prefix_rank(m_v->data(), idx, v);
		// else
		// 	return superblock + block;

		// uint64_t word_pos = (2 * (t_v - 1) * (((idx * t_b) >> 9))) + 2 * v;
		// const uint64_t* p = m_basic_block.data() + word_pos; // 2*(idx*t_b/512) + 2*v
		//
		// // if (idx == 512 && v == 0) {
		// // 	std::cout << "rank(" << idx << ", " << (unsigned) v
		// // 			  << ") - sb: " << *p
		// // 			  //<< " shift by " << ((8 * (((idx * t_b) & 0x1FF) / 64))) << " ... "
		// // 			  << ", b : " << ((*(p + 1) >> (8 * (((idx * t_b) & 0x1FF) / 64))) & 0b11111111);
		// // 	if (idx & (0x1F))
		// // 		std::cout << ", pc: " << trait_type::word_rank(m_v->data(), idx, v);
		// // 	std::cout << std::endl;
		// // 	std::cout << "word_pos: " << word_pos << "\n";
		// // }
		//
		// // TODO: test effect of likely/unlikely
		// if (likely(idx & (0x1F))) // if (idx % 32 != 0) nur für DNA-alphabet
		// 	return *p + ((*(p + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF) +
		// 		   word_prefix_rank(m_v->data(), idx, v);
		// else
		// 	return *p + ((*(p + 1) >> ((((idx * t_b) & 0x1FF) >> 6) << 3)) & 0xFF);
	}

	size_type size() const { return m_v->size(); }

	size_type
	serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		size_type			 written_bytes = 0;
		structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		written_bytes += m_block.serialize(out, child, "prefix_block_counts");
		written_bytes += m_superblock.serialize(out, child, "prefix_superblock_counts");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	void load(std::istream& in, const int_vector<>* v = nullptr)
	{
		m_v = v;
		m_block.load(in);
		m_superblock.load(in);
		init(v, 0);
	}

	void set_vector(const int_vector<>* v = nullptr) { m_v = v; }
};

} // end namespace sds

#endif // end file
