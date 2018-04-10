// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_int.hpp
    \brief rank_support_int.hpp contains classes that support a sdsl::int_vector with constant time rank information.
	\author Christopher Pockrandt
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_INT
#define INCLUDED_SDSL_RANK_SUPPORT_INT

/** \defgroup rank_support_group Rank Support (RS)
 * This group contains data structures which support an sdsl::bit_vector with the rank method.
 */

 #include "int_vector.hpp"
 #include "uint128_t.hpp"

#include <bitset>

//! Namespace for the succinct data structure library.
namespace sdsl {

//! The base class of classes supporting rank_queries for a sdsl::bit_vector in constant time.
/*!
*/
template <uint8_t t_b, uint8_t t_v>
class rank_support_int {
protected:
	const int_vector<t_b>* m_v; //!< Pointer to the rank supported bit_vector

public:
	typedef typename int_vector<t_b>::size_type size_type;
	typedef typename int_vector<t_b>::value_type value_type;

	//! Constructor
	/*! \param v The supported bit_vector.
         */
	rank_support_int(const int_vector<t_b>* v = nullptr);
	//! Copy constructor
	rank_support_int(const rank_support_int&) = default;
	rank_support_int(rank_support_int&&)	  = default;
	rank_support_int& operator=(const rank_support_int&) = default;
	rank_support_int& operator=(rank_support_int&&) = default;
	//! Destructor
	virtual ~rank_support_int() {}

	//! Answers rank queries for the supported bit_vector.
	/*!	\param i Argument for the length of the prefix v[0..i-1].
        	\returns Number of 1-bits in the prefix [0..i-1] of the supported bit_vector.
        	\note Method init has to be called before the first call of rank.
        	\sa init
         */
	virtual size_type rank(size_type i, const value_type v) const = 0;
	//! Alias for rank(i)
	virtual size_type operator()(size_type idx, const value_type v) const = 0;
	//! Serializes rank_support.
	/*! \param out Out-Stream to serialize the data to.
        */
	virtual size_type
	serialize(std::ostream& out, structure_tree_node* v, std::string name) const = 0;
	//! Loads the rank_support.
	/*! \param in In-Stream to load the rank_support data from.
            \param v The supported bit_vector.
         */
	virtual void load(std::istream& in, const int_vector<t_b>* v = nullptr) = 0;
	//! Sets the supported bit_vector to the given pointer.
	/*! \param v The new bit_vector to support.
         *  \note Method init has to be called before the next call of rank.
         *  \sa init, rank
         */
	virtual void set_vector(const int_vector<t_b>* v = nullptr) = 0;
};

template <uint8_t t_b, uint8_t t_v>
inline rank_support_int<t_b, t_v>::rank_support_int(const int_vector<t_b>* v) { m_v = v; }

//----------------------------------------------------------------------

template <typename uintX_t>
constexpr uintX_t bm_rec(const uintX_t w, const uint8_t length, const uint8_t max_length)
{
    return (length >= max_length) ? w
        : bm_rec(w | (w << length), length << 1, max_length);
}

// template <uint8_t t_bits, uint8_t t_value>
// constexpr uint64_t create_bitmask()
// {
//     return create_bitmask_recursive(0);
// }

template <uint8_t t_b, uint8_t t_v, typename T = void>
struct rank_support_int_trait;

template <uint8_t t_b, uint8_t t_v>
struct rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b == 0> > {
protected:
	static uint64_t masks[t_v];
	static uint64_t even_mask;
	static uint64_t carry_select_mask;

public:
	typedef typename rank_support_int<t_b, t_v>::size_type size_type;
	typedef typename rank_support_int<t_b, t_v>::value_type value_type;

	static void init()
	{
		even_mask = bits::lo_set[t_b];
		for (uint8_t i = t_b * 2; i < 64; i <<= 1)
			even_mask |= even_mask << i;

		// even_mask = 0;//bm_rec(bits::lo_set[t_b], t_b * 2 + 3, 64);

		// std::cout << std::bitset<64>(even_mask) << std::endl <<  << std::endl << std::endl;

		for (value_type v = 0; v < t_v; ++v)
		{
			masks[v] = v;
			for (uint8_t i = t_b * 2; i < 64; i <<= 1)
				masks[v] |= masks[v] << i;
		}

		carry_select_mask = masks[1] << t_b;

		uint64_t tmp_carry = masks[1];
		for (value_type v = 0; v < t_v; ++v)
			masks[v] |= tmp_carry << t_b;
	}

    static size_type set_positions(uint64_t w, const value_type v)
    {
		uint64_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
    }

	static uint32_t word_rank(const uint64_t* data, const size_type idx, const value_type v)
	{
		size_type bit_pos = idx * t_b;
		uint64_t w = *(data + (bit_pos >> 6));
		return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
	}

	static uint32_t full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	{
		uint64_t w = *(data + word_pos);
		return bits::cnt(set_positions(w, v));
	}
};

template <uint8_t t_b, uint8_t t_v>
struct rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b != 0> > {
protected:
	static constexpr uint8_t cyclic_shifts = t_b / util::gcd(t_b, 64 % t_b);
	static uint8_t offsets[cyclic_shifts];
	static uint128_t masks[t_v];
	static uint128_t even_mask;
	static uint128_t carry_select_mask;

public:
	typedef typename rank_support_int<t_b, t_v>::size_type size_type;
	typedef typename rank_support_int<t_b, t_v>::value_type value_type;

	static void init()
	{
		even_mask = bits::lo_set[t_b];
		for (uint8_t i = t_b * 2; i < 128; i <<= 1)
			even_mask |= even_mask << i;

		for (value_type v = 0; v < t_v; ++v)
		{
			masks[v] = v;
			for (uint8_t i = t_b * 2; i < 128; i <<= 1)
				masks[v] |= masks[v] << i;
		}

		for (uint8_t i = 1; i < cyclic_shifts; ++i)
			offsets[i] = t_b - ((64 - offsets[i - 1]) % t_b);

		carry_select_mask = masks[1] << t_b;

		uint64_t tmp_carry = masks[1];
		for (value_type v = 0; v < t_v; ++v)
			masks[v] |= tmp_carry << t_b;
	}

    static size_type set_positions(uint64_t w, const value_type v)
    {
		// TODO: 128bit-masken auf 64bit casten?
		uint64_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
    }

    static uint128_t set_positions(uint128_t w, const value_type v, const size_type word_pos)
    {
		uint8_t offset = 64 - (t_b - offsets[word_pos % cyclic_shifts]);
		w >>= offset;
		uint8_t offset_for_even_and_carry = t_b - offsets[1] + ((64 / t_b) & 0x01) * t_b;
		uint128_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
    }

	static void set(uint128_t & w, const uint64_t hi, const uint64_t lo)
	{
		w = hi;
		w <<= 64;
		w |= lo;
	}

	static uint32_t word_rank(const uint64_t* data, const size_type idx, const value_type v)
	{
		size_type bit_pos = idx * t_b;
		size_type word_pos = bit_pos >> 6;
		if (word_pos % cyclic_shifts == 0)
		{
			uint64_t w = *(data + word_pos);
			return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
		}
		else
		{
			size_type in_word_pos = bit_pos & 0x3F;

			uint128_t w;
			set(w, *(data + word_pos), *(data + word_pos - 1));

			uint8_t w_width = in_word_pos + (t_b - offsets[(word_pos) % cyclic_shifts]) + 1; // + 1 because carry bit will be one position further to the left

			uint128_t lo_set128;
			set(lo_set128, bits::lo_set[std::max((signed) 0, w_width - 64)], bits::lo_set[std::min((uint8_t)64, w_width)]);
			uint128_t res = set_positions(w, v, word_pos) & lo_set128;
			return bits::cnt(res >> 64) + bits::cnt(res);
		}
	}

	// NOTE(cpockrandt): we only count full blocks!
	static uint32_t full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	{
		if (word_pos % cyclic_shifts == 0)
		{
			uint64_t w = *(data + word_pos);
			return bits::cnt(set_positions(w, v));
		}
		else
		{
			uint128_t w;
			set(w, *(data + word_pos), *(data + word_pos - 1));

			// uint8_t w_width = 64 + (t_b - offsets[word_pos % cyclic_shifts]) - (t_b - offsets[(word_pos + 1) % cyclic_shifts]);
			// wrong for last cyclic shift: uint8_t w_width = 64 - offsets[word_pos % cyclic_shifts] + offsets[(word_pos + 1) % cyclic_shifts] + 1;
			uint8_t w_width = 64 + (t_b - offsets[word_pos % cyclic_shifts]) + 1; // + 1 because carry bit will be one position further to the left
			if ((word_pos + 1) % cyclic_shifts) // offsets[(word_pos + 1) % cyclic_shifts] > 0
			 	w_width -= (t_b - offsets[(word_pos + 1) % cyclic_shifts]);

			uint128_t lo_set128;
			set(lo_set128, bits::lo_set[std::max((signed) 0, w_width - 64)], bits::lo_set[std::min((uint8_t)64, w_width)]);
			uint128_t res = set_positions(w, v, word_pos) & lo_set128;
			return bits::cnt(res >> 64) + bits::cnt(res);
		}
	}
};

template <uint8_t t_b, uint8_t t_v>
uint64_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b == 0> >::masks[t_v];
template <uint8_t t_b, uint8_t t_v>
uint128_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b != 0> >::masks[t_v];

template <uint8_t t_b, uint8_t t_v>
uint8_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b != 0> >::offsets[rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b != 0> >::cyclic_shifts];

template <uint8_t t_b, uint8_t t_v>
uint64_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b == 0> >::even_mask;
template <uint8_t t_b, uint8_t t_v>
uint128_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b != 0> >::even_mask;

template <uint8_t t_b, uint8_t t_v>
uint64_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b == 0> >::carry_select_mask;
template <uint8_t t_b, uint8_t t_v>
uint128_t rank_support_int_trait<t_b, t_v, std::enable_if_t<64 % t_b != 0> >::carry_select_mask;

} // end namespace sdsl

// #include "rank_support_int_v.hpp"
#include "rank_support_int_scan.hpp"

#endif // end file
