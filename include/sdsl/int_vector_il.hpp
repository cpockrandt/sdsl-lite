// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*!\file bit_vector_il.hpp
   \brief bit_vector_il.hpp contains the sdsl::bit_vector_il class, and
          classes which support rank and select for bit_vector_il.
   \author Matthias Petri, Simon Gog
*/
#ifndef SDSL_INT_VECTOR_IL
#define SDSL_INT_VECTOR_IL

#include "int_vector.hpp"
#include "util.hpp"
#include "iterators.hpp"
#include "rank_support_int.hpp"
#include <bitset>

//! Namespace for the succinct data structure library
namespace sdsl {

// forward declaration needed for friend declaration
class rank_support_int_il;							// in bit_vector_il

//! A bit vector which interleaves the original bit_vector with rank information.
/*!
 * This class is a uncompressed bit vector representation. It copies the original
 * bit_vector and interleaves the data every t_bs bits with a cumulative
 * sum of set bits before the current position. Each cumulative sum is stored
 * in a 64 bit word.
 *
 * \tparam t_bs Block size in bits. t_bs has to be a power of 2 and t_bs >= 64.
 */

class int_vector_il {

public:
	typedef int_vector<>::size_type						size_type;
	typedef size_type									value_type;
	// typedef bit_vector::difference_type					difference_type;
	// typedef random_access_const_iterator<bit_vector_il> iterator;
	// typedef iterator									const_iterator;
	// typedef bv_tag										index_category;

	friend class rank_support_int_il;

	typedef rank_support_int_il   rank_1_type;

private:
	size_type	  m_size		 = 0; //!< Size of the original bitvector
	int_vector<64> m_data;		   //!< Data container

protected:
	uint8_t t_b;
	uint8_t t_v;
	uint64_t even_mask;
	uint64_t carry_select_mask;
	std::vector<uint64_t> masks;

	template <typename uintX_t>
	uintX_t bm_rec(const uintX_t w, const uint8_t length, const uint8_t max_length) const
	{
		return (length >= max_length) ? w : bm_rec(w | (w << length), length << 1, max_length);
	}

	uint64_t set_positions(uint64_t w, const value_type v) const
	{
		uint64_t w_even = even_mask & w;
		uint64_t w_odd = even_mask & (w >> t_b);
		uint64_t res = ((masks[v] - w_even) & ~(masks[v - 1] - w_even)) & carry_select_mask;
		res |= (((masks[v] - w_odd) & ~(masks[v-1] - w_odd)) & carry_select_mask) << 1;
		return res;
	}

	uint64_t set_positions_prefix(uint64_t w, const value_type v) const
	{
		uint64_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
	}

	// assumptions?
	uint32_t word_rank(const uint64_t* data, const size_type idx, const value_type v) const
	{
		size_type bit_pos = idx * t_b;
		uint64_t w = *(data + (bit_pos >> 6));
		return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
	}

	// assumptions?
	uint32_t full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v) const
	{
		uint64_t w = *(data + word_pos);
		return bits::cnt(set_positions(w, v));
	}

	// assumptions?
	uint32_t word_prefix_rank(const uint64_t* data, const size_type idx, const value_type v) const
	{
		size_type bit_pos = idx * t_b;
		uint64_t w = *(data + (bit_pos >> 6));
		return bits::cnt(set_positions_prefix(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
	}

	// assumptions?
	uint32_t full_word_prefix_rank(const uint64_t* data, const size_type word_pos, const value_type v) const
	{
		uint64_t w = *(data + word_pos);
		return bits::cnt(set_positions_prefix(w, v));
	}

	void init(/*const int_vector<>* v, unsigned max_val*/)
	{
		// if (v != nullptr) {
			// m_v = v;
			t_b = 2;//m_v->width();
			// max_val = 3;
			// if (max_val == 0)
			// {
			// 	for (unsigned i = 0; i < m_v->size(); ++i)
			// 		if ((*m_v)[i] > max_val)
			// 			max_val = (*m_v)[i];
			// }

			t_v = 1ULL << t_b;//max_val + 1;

			even_mask = bm_rec<uint64_t>(bits::lo_set[t_b], t_b * 2, 64);
			carry_select_mask = bm_rec<uint64_t>(1ULL << t_b, t_b * 2, 64);

			masks.resize(t_v);
			for (value_type v = 0; v < t_v; ++v)
			{
				masks[v] = v;
				for (uint8_t i = t_b * 2; i < 64; i <<= 1)
					masks[v] |= masks[v] << i;
			}

			uint64_t tmp_carry = masks[1];
			for (value_type v = 0; v < t_v; ++v)
				masks[v] |= tmp_carry << t_b;

			masks.shrink_to_fit();
		// }
	}

public:
	int_vector_il() {}
	int_vector_il(const int_vector_il&) = default;
	int_vector_il(int_vector_il&&)		= default;
	int_vector_il& operator=(const int_vector_il&) = default;
	int_vector_il& operator=(int_vector_il&&) = default;

	int_vector_il(const int_vector<>& bv)
	{
		m_size = bv.size();
		m_data.bit_resize(2*((bv.bit_size() + 256) / 256) * 256);

		size_type blocks = (bv.size() + 31) / 32; // assume width of 2
		const uint64_t* bvp = bv.data();

		init();

		uint64_t supercounts[3] = {0};
		size_type i = 0, j = 0;

		while (i < blocks) {
			uint64_t counts[3][3] = {0};
			while (i < blocks) {
				m_data[j] = bv.get_int(i * 64, 64);

				if (j % 4 != 3) {
					for (unsigned v = 0; v < 3; ++v) {
						if ((j & 3) > 0)
							counts[v][j & 3] += counts[v][(j & 3) - 1];
						counts[v][j & 3] += bits::cnt(set_positions_prefix(m_data[j], v));
					}
				}

				// if (j % 4 != 3) {
				// 	for (unsigned v = 0; v < 3; ++v)
				// 		std::cout << "counts[" << v << "][" << (j & 3) << "] = " << counts[v][j & 3] << "\n";
				// }

				++i; ++j;
				if (!(j & 3)) // j % 4 == 0
				 	break;
			}

			while (j & 3) // j % 4 != 0
				++j;

			m_data[j    ] = supercounts[0];
			m_data[j + 1] = supercounts[1];
			m_data[j + 2] = supercounts[2];
			m_data[j + 3] = 0ULL;

			m_data[j + 3] |= counts[2][2] << (64 - 1*7);
			m_data[j + 3] |= counts[2][1] << (64 - 2*7);
			m_data[j + 3] |= counts[2][0] << (64 - 3*7);
			m_data[j + 3] |= counts[1][2] << (64 - 4*7);
			m_data[j + 3] |= counts[1][1] << (64 - 5*7);
			m_data[j + 3] |= counts[1][0] << (64 - 6*7);
			m_data[j + 3] |= counts[0][2] << (64 - 7*7);
			m_data[j + 3] |= counts[0][1] << (64 - 8*7);
			m_data[j + 3] |= counts[0][0] << (64 - 9*7);

			supercounts[0] += counts[0][2] + bits::cnt(set_positions_prefix(m_data[j - 1], 0));
			supercounts[1] += counts[1][2] + bits::cnt(set_positions_prefix(m_data[j - 1], 1));
			supercounts[2] += counts[2][2] + bits::cnt(set_positions_prefix(m_data[j - 1], 2));
			j += 4;
		}

		if (blocks % 4 == 0 && m_size % 128 == 0) {
			j += 4;
			m_data[j    ] = supercounts[0];
			m_data[j + 1] = supercounts[1];
			m_data[j + 2] = supercounts[2];
			m_data[j + 3] = 0ULL;
		}

		//for (unsigned i = 0; i < m_data.size(); ++i)
		//	std::cout << "--- " << std::bitset<64>(m_data[i]) << "     " << m_data[i] << "\n";
	}

	//! Accessing the i-th element of the original bit_vector
	/*! \param i An index i with \f$ 0 \leq i < size()  \f$.
         *  \return The i-th bit of the original bit_vector
         *  \par Time complexity
         *     \f$ \Order{1} \f$
         */
	value_type operator[](size_type i) const
	{
		assert(i < m_size);
		// assume width of 2
		size_t pos = ((i >> 7) << 3) + ((i & 127) >> 5);
		size_t id = (i & 31) << 1;
		return (m_data[pos] >> id) & 3ULL;
		// return m_data.get_int(..., 2);
	}

	//! Get the integer value of the binary string of length len starting at position idx.
	/*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *   \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
	// uint64_t get_int(size_type idx, uint8_t len = 64) const
	// {
	// 	assert(idx + len - 1 < m_size);
	// 	size_type bs	  = idx >> m_block_shift;
	// 	size_type b_block = bs + (idx >> 6) + 1;
	// 	bs				  = (idx + len - 1) >> m_block_shift;
	// 	size_type e_block = bs + ((idx + len - 1) >> 6) + 1;
	// 	if (b_block == e_block) { // spans on block
	// 		return (m_data[b_block] >> (idx & 63)) & bits::lo_set[len];
	// 	} else { // spans two blocks
	// 		uint8_t b_len = 64 - (idx & 63);
	// 		return (m_data[b_block] >> (idx & 63)) |
	// 			   (m_data[e_block] & bits::lo_set[len - b_len]) << b_len;
	// 	}
	// }

	//! Returns the size of the original bit vector.
	size_type size() const { return m_size; }

	//! Serializes the data structure into the given ostream
	size_type
	serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		size_type			 written_bytes = 0;
		written_bytes += write_member(m_size, out, child, "size");
		written_bytes += m_data.serialize(out, child, "data");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	//! Loads the data structure from the given istream.
	void load(std::istream& in)
	{
		read_member(m_size, in);
		m_data.load(in);
		init();
	}

	// iterator begin() const { return iterator(this, 0); }
	//
	// iterator end() const { return iterator(this, size()); }


	bool operator==(const int_vector_il& v) const
	{
		return m_size == v.m_size && m_data == v.m_data;
	}

	bool operator!=(const int_vector_il& v) const { return !(*this == v); }
};

class rank_support_int_il {

public:
	typedef int_vector<>::size_type size_type;
	typedef size_type value_type;
	typedef int_vector_il   int_vector_type;

	friend class int_vector_il;

private:
	const int_vector_type* m_v;

public:
	rank_support_int_il(const int_vector_type* v = nullptr)
	{
		set_vector(v);
		// m_block_shift	= bits::hi(t_bs);
		// m_block_mask	 = t_bs - 1;
		// m_block_size_U64 = bits::hi(t_bs >> 6);
	}

	size_type rank(size_type idx, const value_type v) const
	{
		assert(m_v != nullptr);
		assert(idx <= m_v->size());

		if (unlikely(v == 0))
			return prefix_rank(idx, v);

		size_type block_512 = (idx >> 7) << 3; // (i / 128) * 8, beginning of one 512 bit block (256 bit bv and 256 bit rank)
		size_type block_pos = (idx & 127) >> 5;
		size_type p = m_v->m_data[block_512 + block_pos]; // (i % 128) / 32

		size_type result;
		if (unlikely(v == m_v->t_v - 1)) // TODO: test effect of likely/unlikely
			result = idx;
		else { // (likely(v != m_v->t_v - 1)) {
		 	result = m_v->m_data[block_512 + 4 + v];
			if (block_pos > 0) // sb[v] + b[v] - b[v-1]
				result += ((m_v->m_data[block_512 + 7] >> (7 * block_pos - 6 + 21 * v)) & 0x7F); // (1 + 7 * ((block_pos - 1) + 3 * v)) .......... get last 7 bits (127 in hex)
		}
		result -= m_v->m_data[block_512 + 3 + v]; // superblockcount of v-1
		if (block_pos > 0)
			result -= (m_v->m_data[block_512 + 7] >> (7 * block_pos + 21 * v - 27)) & 0x7F;

		uint8_t pos_in_word = (idx * 2) & 0x3F;
		if (likely(pos_in_word)) { // TODO: ein word_rank!
			if (likely(v != m_v->t_v - 1)) // if (idx % 32 != 0) nur für DNA-alphabet
				// result += word_rank(m_v->data(), idx, v);
				result += bits::cnt(m_v->set_positions(p, v) & bits::lo_set[pos_in_word + 1]);
			else
				// result -= word_prefix_rank(m_v->data(), idx, v - 1);
				result -= bits::cnt(m_v->set_positions_prefix(p, v - 1) & bits::lo_set[pos_in_word + 1]);
		}
		return result;
	}

	//! Returns the position of the i-th occurrence in the bit vector.
	size_type prefix_rank(size_type idx, const value_type v) const
	{
		assert(m_v != nullptr);
		assert(idx <= m_v->size());

		if (unlikely(v == m_v->t_v - 1)) // TODO: test effect of likely/unlikely
			return idx;

		size_type block_512 = (idx >> 7) << 3; // (i / 128) * 8, beginning of one 512 bit block (256 bit bv and 256 bit rank)
		size_type block_pos = (idx & 127) >> 5;
		size_type p = m_v->m_data[block_512 + block_pos]; // (i % 128) / 32
		size_type superblockcount = m_v->m_data[block_512 + 4 + v];
		size_type blockcount = 0;
		if (likely(block_pos > 0))
			blockcount = (m_v->m_data[block_512 + 7] >> (7 * block_pos + 21 * v - 6)) & 0x7F; // (1 + 7 * ((block_pos - 1) + 3 * v)) .... get last 7 bits

		uint8_t pos_in_word = (idx * 2) & 0x3F;
		if (likely(pos_in_word)) // if (idx % 32 != 0) nur für DNA-alphabet
			return superblockcount + blockcount
				+ bits::cnt(m_v->set_positions_prefix(p, v) & bits::lo_set[pos_in_word + 1]);
		else
			return superblockcount + blockcount;
	}

	size_type operator()(size_type i, value_type v) const { return rank(i, v); }

	size_type size() const { return m_v->size(); }

	void set_vector(const int_vector_type* v = nullptr) { m_v = v; }

	rank_support_int_il& operator=(const rank_support_int_il& rs)
	{
		if (this != &rs) {
			set_vector(rs.m_v);
		}
		return *this;
	}

	void load(std::istream&, const int_vector_type* v = nullptr) { set_vector(v); }

	size_type
	serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		return serialize_empty_object(out, v, name, this);
	}
};

} // end namespace sdsl
#endif
