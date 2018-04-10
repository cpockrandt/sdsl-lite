// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_int.hpp
    \brief rank_support_int.hpp contains classes that support a sdsl::int_vector with constant time rank information.
	\author Simon Gog
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

// constexpr uint64_t create_bitmask_recursive(unsigned const bitsUsed, unsigned const blocksLeft, unsigned const blocksProcessed, unsigned const blocksize, TWord const value1, TWord const value2, TWord const result)
// {
//     return (blocksLeft == 0)
//         ? result
//         : create_bitmask_recursive(bitsUsed,
//                                    blocksLeft - 1,
//                                    blocksProcessed + 1,
//                                    blocksize,
//                                    value1,
//                                    value2,
//                                    result | (((blocksProcessed % 2) ? value1 : value2)
//                                            << (bitsUsed - blocksProcessed * blocksize)));
// }
//
// template <uint8_t t_bits, uint8_t t_value>
// constexpr uint64_t create_bitmask()
// {
//     return create_bitmask_recursive(0);
// }

// https://stackoverflow.com/questions/30953248/why-doesnt-sfinae-enable-if-work-for-member-functions-of-a-class-template

template <uint8_t t_b, uint8_t t_v>
struct rank_support_int_trait {
protected:
	static constexpr uint8_t cyclic_shifts = t_b / util::gcd(t_b, 64 % t_b);
	static uint8_t offsets[cyclic_shifts + 1]; // TODO(cpockrandt) + 1
	static uint64_t masks[t_v]/*[cyclic_shifts + 1]*/; // TODO(cpockrandt) + 1
	static uint64_t even_mask; // TODO: 128bit even_mask
	static uint64_t carry_select_mask; // TODO: 128bit

public:
	typedef typename rank_support_int<t_b, t_v>::size_type size_type;
	typedef typename rank_support_int<t_b, t_v>::value_type value_type;

	static void init()
	{
		even_mask = bits::lo_set[t_b];
		for (uint8_t i = t_b * 2; i < 64; i <<= 1)
		{
			even_mask |= even_mask << i;
		}

		for (value_type v = 0; v < t_v; ++v)
		{
			masks[v] = v;
			for (uint8_t i = t_b * 2; i < 64; i <<= 1)
			{
				masks[v] |= masks[v] << i;
			}
			// uint8_t tmp_n = 0;
			// util::cyclic_shifts(masks[v], tmp_n, v, t_b, offsets);
			// for (unsigned i = 0; i < cyclic_shifts; ++i)
			// {
				masks[v]/*[i]*/ &= even_mask;
			// }
		}

		uint64_t tmp_ignore[cyclic_shifts + 1];
		uint8_t tmp_n;
		util::cyclic_shifts(tmp_ignore, tmp_n, 0, t_b, offsets);

		// std::cout << "offsets: ";
		// for (unsigned i = 0; i < cyclic_shifts; ++i)
		// 	std::cout << (unsigned) offsets[i] << ", ";
		// std::cout << std::endl;

		carry_select_mask = masks[1]/*[0]*/ << t_b;

		uint64_t tmp_carry = masks[1]/*[0]*/; // TODO: offi
		for (value_type v = 0; v < t_v; ++v)
		{
			for (unsigned i = 0; i < cyclic_shifts; ++i)
			{
				masks[v]/*[i]*/ |= tmp_carry << t_b;
				// std::cout << std::bitset<64>(masks[v][i]) << std::endl;
			}
			// std::cout << "..............................................................................\n";
		}
	}

    static size_type set_positions(uint64_t w, const value_type v)
    {
		// std::cout << "w      : " << std::bitset<64>(w) << std::endl;
		// std::cout << "mask   : " << std::bitset<64>(masks[v]) << std::endl;
		// std::cout << "w_even : " << std::bitset<64>(even_mask & w) << std::endl;
		// std::cout << "diff   : " << std::bitset<64>(masks[v]/*[0]*/ - (even_mask & w)) << std::endl;
		// std::cout << "c_mask : " << std::bitset<64>(carry_select_mask) << std::endl;
		// std::cout << "res1   : " << std::bitset<64>((masks[v]/*[0]*/ - (even_mask & w)) & carry_select_mask) << std::endl;

		uint64_t res = ((masks[v]/*[0]*/ - (even_mask & w)) & carry_select_mask);
		res |= ((masks[v]/*[0]*/ - (even_mask & (w >> t_b))) & carry_select_mask) << 1;

		// std::cout << "w>>t_b : " << std::bitset<64>(w >> t_b) << std::endl;
		// std::cout << "mask   : " << std::bitset<64>(masks[v]) << std::endl;
		// std::cout << "w_even : " << std::bitset<64>(even_mask & (w >> t_b)) << std::endl;
		// std::cout << "diff   : " << std::bitset<64>(masks[v]/*[0]*/ - (even_mask & (w >> t_b))) << std::endl;
		// std::cout << "c_mask : " << std::bitset<64>(carry_select_mask) << std::endl;
		// std::cout << "res2   : " << std::bitset<64>((masks[v]/*[0]*/ - (even_mask & (w >> t_b))) & carry_select_mask) << std::endl;
		// std::cout << "res2<<1: " << std::bitset<64>(((masks[v]/*[0]*/ - (even_mask & (w >> t_b))) & carry_select_mask) << 1) << std::endl;
		// std::cout << "res    : " << std::bitset<64>(res) << std::endl;

		return res;
    }

    static uint128_t set_positions(uint128_t w, const value_type v, const size_type word_pos, const bool debug = false)
    {
		uint8_t offset = 64 - (t_b - offsets[word_pos % cyclic_shifts]);
		w >>= offset;
		if (debug)
		{
			std::cout << "w (offset " << (unsigned) offset << "): "; print(w);
		}

		uint8_t offset_for_even_and_carry = t_b - offsets[1] + ((64 / t_b) & 0x01) * t_b;

		uint128_t masks128;
		set(masks128, /*masks[v][1]*/masks[v] >> offset_for_even_and_carry, masks[v]/*[0]*/);
		if (debug)
		{
			std::cout << "masks        : "; print(masks128);
		}

		uint128_t even_mask128;
		set(even_mask128, even_mask >> offset_for_even_and_carry, even_mask);
		if (debug)
		{
			std::cout << "even_mask    : "; print(even_mask128);
		}

		uint128_t carry_select_mask128;
		set(carry_select_mask128, carry_select_mask >> offset_for_even_and_carry, carry_select_mask);
		if (debug)
		{
			std::cout << "carry_select : "; print(carry_select_mask128);
		}

		uint128_t res1 = ((masks128 - (even_mask128 & w)) & carry_select_mask128);
		if (debug)
		{
			std::cout << "res1         : "; print(res1);
		}

		uint128_t res2 = ((masks128 - (even_mask128 & (w >> t_b))) & carry_select_mask128);
		if (debug)
		{
			std::cout << "res2         : "; print(res2);
		}

		uint128_t res = res1 | (res2 << 1);
		if (debug)
		{
			std::cout << "res          : "; print(res);
		}

		return res;
    }

	static int print(uint128_t x)
	{
		std::cout << std::bitset<64>(x >> 64) << " " << std::bitset<64>(x) << std::endl;
	}

	static void set(uint128_t & w, const uint64_t hi, const uint64_t lo)
	{
		w = hi;
		w <<= 64;
		w |= lo;
	}

	static std::enable_if_t<64 % t_b != 0, uint32_t> word_rank(const uint64_t* data, const size_type idx, const value_type v)
	{
		bool debug = false;
		if (debug)
			std::cout << ".................................................................\nword_rank\n";

		size_type bit_pos = idx * t_b;
		size_type word_pos = bit_pos >> 6;
		if (word_pos % cyclic_shifts == 0)
		{
			uint64_t w = *(data + word_pos);
			return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
		}
		else
		{
			// for (unsigned xxx = 0; xxx < 3; ++xxx)
			// 	std::cout << std::bitset<64>(*(data + xxx)) << std::endl;
			if (debug)
				std::cout << "Complex case.\n";
			//std::cout << "div          : " << std::bitset<64>(~0ULL) << std::bitset<64>(0ULL) << std::endl;
			//std::cout << "w manual     : " << std::bitset<64>(*(data + word_pos)) << std::bitset<64>(*(data + word_pos - 1)) << std::endl;

			size_type in_word_pos = bit_pos & 0x3F;

			uint128_t w;
			set(w, *(data + word_pos), *(data + word_pos - 1));
			if (debug)
			{
				std::cout << "w            : "; print(w);
			}

			// uint8_t offset_shift = 64 - (t_b - offsets[(word_pos) % cyclic_shifts]);
			// if (debug)
			// 	std::cout << "offset_shift : " << (unsigned) offset_shift << std::endl;
			// w >>= offset_shift;
			//
			// if (debug)
			// 	std::cout << "w nach shift : "; print(w);

			if (debug)
			{
				std::cout << "in_word_pos  : " << in_word_pos << std::endl;
			}

			// uint8_t w_width = 64 + (t_b - offsets[word_pos % cyclic_shifts]) - (t_b - offsets[(word_pos + 1) % cyclic_shifts]);
			uint8_t w_width = in_word_pos + (t_b - offsets[(word_pos) % cyclic_shifts]) + 1;
			if (debug)
			{
				std::cout << "w_width      : " << (unsigned) w_width << std::endl; // added plus one because:
			}
			// even:   0111000111000111000111000111000111000111000111000111000111000111
			// lo_set: 0111111111111111111111111111111111111111111111111111111111111111
			// but carry bit will be on the very first position!

			uint128_t lo_set128;
			set(lo_set128, bits::lo_set[std::max((signed) 0, w_width - 64)], bits::lo_set[std::min((uint8_t)64, w_width)]);
			if (debug)
			{
				std::cout << "lo_set       : "; print(lo_set128);
			}

			uint128_t res = set_positions(w, v, word_pos, debug) & lo_set128;
			if (debug)
			{
				std::cout << "lo_set & res : "; print(res);
			}

			if (debug)
			{
				std::cout << "popcount(res): " << (bits::cnt(res >> 64) + bits::cnt(res)) << std::endl;
			}

			return bits::cnt(res >> 64) + bits::cnt(res);
		}
	}

	// NOTE(cpockrandt): we only count full blocks!
	static std::enable_if_t<64 % t_b != 0, uint32_t> full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	{
		bool debug = false;
		if (debug)
			std::cout << ".................................................................\nfull_word_rank\n";

		if (word_pos % cyclic_shifts == 0)
		{
			// word_pos can also be 0
			uint64_t w = *(data + word_pos);
			if (debug)
				std::cout << "Simple case: " << bits::cnt(set_positions(w, v)) << ".\n";
			return bits::cnt(set_positions(w, v));
		}
		else
		{
			// for (unsigned xxx = 0; xxx < 3; ++xxx)
			// 	std::cout << std::bitset<64>(*(data + xxx)) << std::endl;
			if (debug)
				std::cout << "Complex case.\n";
			//std::cout << "div          : " << std::bitset<64>(~0ULL) << std::bitset<64>(0ULL) << std::endl;
			//std::cout << "w manual     : " << std::bitset<64>(*(data + word_pos)) << std::bitset<64>(*(data + word_pos - 1)) << std::endl;

			uint128_t w;
			set(w, *(data + word_pos), *(data + word_pos - 1));
			if (debug)
			{
				std::cout << "w            : "; print(w);
			}

			// uint8_t w_width = 64 + (t_b - offsets[word_pos % cyclic_shifts]) - (t_b - offsets[(word_pos + 1) % cyclic_shifts]);
			// wrong for last cyclic shift: uint8_t w_width = 64 - offsets[word_pos % cyclic_shifts] + offsets[(word_pos + 1) % cyclic_shifts] + 1;
			uint8_t w_width = 64 + (t_b - offsets[word_pos % cyclic_shifts]) + 1;
			if ((word_pos + 1) % cyclic_shifts) // offsets[(word_pos + 1) % cyclic_shifts] > 0
			 	w_width -= (t_b - offsets[(word_pos + 1) % cyclic_shifts]);
			if (debug)
				std::cout << "w_width      : " << (unsigned) w_width << std::endl; // added plus one because:
			// even:   0111000111000111000111000111000111000111000111000111000111000111
			// lo_set: 0111111111111111111111111111111111111111111111111111111111111111
			// but carry bit will be on the very first position!

			uint128_t lo_set128;
			set(lo_set128, bits::lo_set[std::max((signed) 0, w_width - 64)], bits::lo_set[std::min((uint8_t)64, w_width)]);
			if (debug)
			{
				std::cout << "lo_set       : "; print(lo_set128);
			}

			uint128_t res = set_positions(w, v, word_pos, debug) & lo_set128;
			if (debug)
			{
				std::cout << "lo_set & res : "; print(res);
			}

			if (debug)
			{
				std::cout << "popcount(res): " << (bits::cnt(res >> 64) + bits::cnt(res)) << std::endl;
			}

			return bits::cnt(res >> 64) + bits::cnt(res);
		}
	}

	// static std::enable_if_t<64 % t_b == 0, uint32_t> word_rank(const uint64_t* data, const size_type idx, const value_type v)
	// {
	// 	size_type bit_pos = idx * t_b;
	// 	uint64_t w = *(data + (bit_pos >> 6));
	// 	return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
	// }
	//
	// static std::enable_if_t<64 % t_b == 0, uint32_t> full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	// {
	// 	uint64_t w = *(data + word_pos);
	// 	return bits::cnt(set_positions(w, v));
	// }
};

template <uint8_t t_b, uint8_t t_v>
uint64_t rank_support_int_trait<t_b, t_v>::masks[t_v]/*[rank_support_int_trait<t_b, t_v>::cyclic_shifts + 1]*/; // TODO(cpockrandt) + 1

template <uint8_t t_b, uint8_t t_v>
uint8_t rank_support_int_trait<t_b, t_v>::offsets[rank_support_int_trait<t_b, t_v>::cyclic_shifts + 1]; // TODO(cpockrandt) + 1

template <uint8_t t_b, uint8_t t_v>
uint64_t rank_support_int_trait<t_b, t_v>::even_mask;

template <uint8_t t_b, uint8_t t_v>
uint64_t rank_support_int_trait<t_b, t_v>::carry_select_mask;

} // end namespace sdsl

// #include "rank_support_int_v.hpp"
#include "rank_support_int_scan.hpp"

#endif // end file
