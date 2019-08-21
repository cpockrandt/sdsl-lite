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
template <uint8_t alphabet_size, uint8_t words_per_block = 1, uint8_t blocks_per_superblock = 4/*, uint8_t alphabet_size = 0*/>
class rank_support_int_v : public rank_support_int<alphabet_size> {
public:
	typedef int_vector<> int_vector_type;
	typedef typename rank_support_int<alphabet_size>::size_type size_type;
	typedef typename rank_support_int<alphabet_size>::value_type value_type;

private:
	int_vector<0> m_block;
	int_vector<64> m_superblock; // TODO: reduce width (at runtime)? benchmark space consumption and running time

	static constexpr uint64_t values_per_word{static_cast<uint64_t>(64) / rank_support_int<alphabet_size>::t_b};

public:
	explicit rank_support_int_v(const int_vector<>* v = nullptr) : rank_support_int<alphabet_size>(v)
	{
	    static_assert(blocks_per_superblock > 1, "There must be at least two blocks per superblock!");
		constexpr uint8_t t_v_decr{this->t_v - 1};

		if (v == nullptr)
		{
			return;
		}
		else if (v->empty())
		{
			m_block.resize(t_v_decr, 0);
			m_superblock.resize(t_v_decr, 0);
			return;
		}

		constexpr uint64_t words_per_superblock = words_per_block * blocks_per_superblock;

		uint64_t const values_per_block{words_per_block * values_per_word};
		uint64_t const values_per_superblock{blocks_per_superblock * values_per_block};
		uint64_t const new_width{ceil_log2(values_per_superblock)};
		m_block.width(new_width);
		// TODO: could also set block width of superblocks. check running time impact!
        // std::cout << "words_per_block: " << (unsigned)words_per_block << '\n'
        //           << "blocks_per_superblock: " << (unsigned)blocks_per_superblock << std::endl << '\n'
		//           << "Bits per value: " << (unsigned)v->width() << '\n'
		//           << "Values per word: " << values_per_word << '\n'
		//           << "Values per superblock: " << values_per_superblock << '\n'
		//           << "Block width: " << (unsigned)new_width << std::endl;

		// NOTE: number of elements is artificially increased because rank can be called on [size()]
		uint64_t const word_count = ((this->m_v->size() - 1 + 1) / values_per_word) + 1; // equivalent to ceil(m_v->size() / values_per_word)
		uint64_t const block_count = ((word_count - 1) / words_per_block) + 1; // equivalent to ceil(word_count / words_per_block)

		// for each superblock we only need `blocks_per_superblock-1` instead of `blocks_per_superblock` blocks.
		// for the last superblock we can subtract the last unused blocks.
        size_type const blocks_needed = (((block_count - 1) / blocks_per_superblock) + 1) * (blocks_per_superblock - 1) - ((blocks_per_superblock - (block_count % blocks_per_superblock)) % blocks_per_superblock);
		size_type const block_size = blocks_needed * t_v_decr;

		size_type const superblock_size = (((word_count - 1) / words_per_superblock) + 1) * t_v_decr; // ceil(word_count / words_per_superblock) * t_v_decr
		// std::cout << "t_v_decr: " << (unsigned)t_v_decr << '\n'
		//           << "Elements: " << v->size() << '\n'
		//           << "block count: " << block_count << '\n'
		//           << "word count: " << word_count << '\n'
        //           << "blocks_needed: " << blocks_needed << '\n'
        //           << "block_size: " << block_size << '\n'
		//           << "superblock_size: " << superblock_size << std::endl;
		m_block.resize(block_size); // TODO: adjust
		m_superblock.resize(superblock_size); // TODO: adjust

		uint64_t const * data = this->m_v->data();
		std::vector<uint64_t> buf_blocks(t_v_decr, 0); // TODO: get rid of these objects
		std::vector<uint64_t> buf_superblocks(t_v_decr, 0);

		for (uint64_t v = 0; v < t_v_decr; ++v)
			m_superblock[v] = 0;

		for (uint64_t word_id = 0, block_id = 0, superblock_id = t_v_decr; word_id < word_count; ++word_id)
		{
			for (uint64_t v = 0; v < t_v_decr; ++v)
				buf_blocks[v] += this->full_word_prefix_rank(data, word_id, v);

			// counted the values in the last word of the current block
			if (word_id % words_per_block == (words_per_block - 1)) // divisor is constexpr, i.e., modulo operation is expected to be cheap
			{
                if (word_id % words_per_superblock != (words_per_superblock - 1)) // divisor is constexpr, i.e., modulo operation is expected to be cheap
                {
                    if (block_id < m_block.size()) // TODO: bloed
                    {
                        for (uint64_t v = 0; v < t_v_decr; ++v)
                            m_block[block_id + v] = buf_blocks[v];
                        block_id += t_v_decr;
                    }
                }
                else
                { // don't store block information for the last block in the superblock!
                    if (superblock_id < m_superblock.size()) // TODO: bloed
                    {
                        for (uint64_t v = 0; v < t_v_decr; ++v)
                        {
                            buf_superblocks[v] += buf_blocks[v];
                            m_superblock[superblock_id + v] = buf_superblocks[v];
                            buf_blocks[v] = 0; // reset blocks
                        }
                    }
                    superblock_id += t_v_decr;
                }
			}
		}

        // std::cout << "\nBlocks:\n";
        // for (uint64_t i = 0; i < m_block.size(); i += t_v_decr)
		// {
		// 	for (uint64_t v = 0; v < t_v_decr; ++v)
	    //         std::cout << (unsigned)m_block[i + v] << ' ';
        //     std::cout << "| ";
		// }
        // std::cout << "\nSuperBlocks:\n";
        // for (uint64_t i = 0; i < m_superblock.size(); i += t_v_decr)
		// {
		// 	for (uint64_t v = 0; v < t_v_decr; ++v)
        //     	std::cout << (unsigned)m_superblock[i + v] << ' ';
        //     std::cout << "| ";
		// }
        // std::cout << "\n\n";
	}

	rank_support_int_v(const rank_support_int_v&) = default;
	rank_support_int_v(rank_support_int_v&&)	  = default;
	rank_support_int_v& operator=(const rank_support_int_v&) = default;
	rank_support_int_v& operator=(rank_support_int_v&&) = default;

	size_type rank(size_type idx, const value_type v) const
	{
		assert(this->m_v != nullptr);
		assert(idx <= this->m_v->size());

		// TODO: optimize?
		if (unlikely(v == 0))
			return prefix_rank(idx, v);

		return prefix_rank(idx, v) - prefix_rank(idx, v - 1);
	}

	inline size_type operator()(size_type idx, const value_type v) const { return rank(idx, v); }

	size_type prefix_rank(size_type idx, const value_type v) const
	{
		assert(this->m_v != nullptr);
		assert(idx <= this->m_v->size());
        assert(v <= this->t_v);

		if (unlikely(v == this->t_v - 1)) // TODO: test effect of likely/unlikely
			return idx;

		constexpr uint8_t t_v_decr{this->t_v - 1};

		uint32_t const values_per_block{words_per_block * values_per_word};

		// size_type const word_id = idx / values_per_word; // TODO: expensive!
		size_type const block_id = idx / values_per_block;
		size_type const superblock_id = block_id / blocks_per_superblock; // NOTE: former idx / values_per_superblock

		// superblock_id * (blocks_per_superblock - 1) + (word_id % words_per_block); // TODO: expensive!

        size_type res = m_superblock[t_v_decr * superblock_id + v];

        uint16_t const block_id_in_superblock = block_id % blocks_per_superblock;
        if (block_id_in_superblock > 0)
            res += m_block[t_v_decr * superblock_id * (blocks_per_superblock - 1) + (block_id_in_superblock - 1)  * t_v_decr + v];

        if (words_per_block > 1)
        {
            size_type const word_id = idx / values_per_word;
            uint64_t w = word_id - (word_id % words_per_block);
            while (w < word_id)
            {
                res += this->full_word_prefix_rank(this->m_v->data(), w, v);
                ++w;
            }
        }

		if (idx % values_per_block != 0)
			res += this->word_prefix_rank(this->m_v->data(), idx, v);

		return res;
	}

	size_type size() const { return this->m_v->size(); }

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
		this->m_v = v;
		m_block.load(in);
		m_superblock.load(in);
		this->init(v);
	}

	void set_vector(const int_vector<>* v = nullptr) { this->m_v = v; }
};

} // end namespace sdsl

#endif // end file
