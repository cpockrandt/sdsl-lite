// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_int.hpp
    \brief rank_support_int.hpp contains classes that support a sdsl::int_vector with constant time rank information.
           Rank is defined as the number of occurrences lexicographically smaller or equal than a given element up to a
           given position.
	\author Christopher Pockrandt
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_INT
#define INCLUDED_SDSL_RANK_SUPPORT_INT

/** \defgroup rank_support_group Rank Support (RS)
 * This group contains data structures which support an sdsl::int_vector with the rank method.
 */

 #include "int_vector.hpp"
 #include "uint128_t.hpp"

#include <bitset>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

//! Namespace for the succinct data structure library.
namespace sdsl {

//! The base class of classes supporting rank_queries for a sdsl::int_vector in constant time.
/*!
*/
template <uint8_t t_b>
class rank_support_int {
protected:
	const int_vector<t_b>* m_v; //!< Pointer to the rank supported bit_vector
    static_assert(t_b > 1, "rank_support_int and deriving classes are not supported on sdsl::bit_vector. Please use an sdsl::int_vector<t_width> with 1 < t_width < 8.");
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
	//! Answers rank queries for the supported bit_vector.
	/*!	\param i Argument for the length of the prefix v[0..i-1].
        	\returns Number of 1-bits in the prefix [0..i-1] of the supported bit_vector.
        	\note Method init has to be called before the first call of rank.
        	\sa init
         */
	virtual size_type prefix_rank(size_type i, const value_type v) const = 0;
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
	// //! Sets the supported bit_vector to the given pointer.
	// /*! \param v The new bit_vector to support.
    //      *  \note Method init has to be called before the next call of rank.
    //      *  \sa init, rank
    //      */
	virtual void set_vector(const int_vector<t_b>* v = nullptr) = 0;
};

template <uint8_t t_b>
inline rank_support_int<t_b>::rank_support_int(const int_vector<t_b>* v) { m_v = v; }

//----------------------------------------------------------------------

template <typename uintX_t>
constexpr uintX_t bm_rec(const uintX_t w, const uint8_t length, const uint8_t max_length)
{
    return (length >= max_length) ? w : bm_rec(w | (w << length), length << 1, max_length);
}

template <uint8_t t_b, typename T = void>
struct rank_support_int_trait;

template <uint8_t t_b>
struct rank_support_int_trait<t_b, std::enable_if_t<64 % t_b == 0> > {
protected:
    constexpr static uint8_t t_v = 1ULL << t_b;
	static uint64_t masks[t_v];
	constexpr static uint64_t even_mask = bm_rec<uint64_t>(bits::lo_set[t_b], t_b * 2, 64);
	constexpr static uint64_t carry_select_mask = bm_rec<uint64_t>(1ULL << t_b, t_b * 2, 64);

public:
	typedef typename rank_support_int<t_b>::size_type size_type;
	typedef typename rank_support_int<t_b>::value_type value_type;

    // https://stackoverflow.com/questions/2226291/is-it-possible-to-create-and-initialize-an-array-of-values-using-template-metapr
	// TODO(cpockrandt): make bitmasks all constexpr and remove this method
	static void init()
	{
		for (value_type v = 0; v < t_v; ++v)
		{
			masks[v] = v;
			for (uint8_t i = t_b * 2; i < 64; i <<= 1)
				masks[v] |= masks[v] << i;
		}

		uint64_t tmp_carry = masks[1];
		for (value_type v = 0; v < t_v; ++v)
			masks[v] |= tmp_carry << t_b;
	}

    static uint64_t set_positions(uint64_t w, const value_type v)
    {
		// uint64_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		// res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
        uint64_t w_even = even_mask & w;
        uint64_t w_odd = even_mask & (w >> t_b);
		uint64_t res = ((masks[v] - w_even) & ~(masks[v - 1] - w_even)) & carry_select_mask;
		res |= (((masks[v] - w_odd) & ~(masks[v-1] - w_odd)) & carry_select_mask) << 1;
		return res;
    }

    static uint64_t set_positions_prefix(uint64_t w, const value_type v)
    {
		uint64_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
    }

        // assumptions?
	static uint32_t word_rank(const uint64_t* data, const size_type idx, const value_type v)
	{
		size_type bit_pos = idx * t_b;
		uint64_t w = *(data + (bit_pos >> 6));
		return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
	}

    // assumptions?
	static uint32_t full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	{
		uint64_t w = *(data + word_pos);
		return bits::cnt(set_positions(w, v));
	}

    // assumptions?
	static uint32_t word_prefix_rank(const uint64_t* data, const size_type idx, const value_type v)
	{
		size_type bit_pos = idx * t_b;
		uint64_t w = *(data + (bit_pos >> 6));
		return bits::cnt(set_positions_prefix(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
	}

    // assumptions?
	static uint32_t full_word_prefix_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	{
		uint64_t w = *(data + word_pos);
		return bits::cnt(set_positions_prefix(w, v));
	}
};

template <uint8_t t_b>
struct rank_support_int_trait<t_b, std::enable_if_t<64 % t_b != 0> > {
protected:
    constexpr static uint8_t t_v = 1ULL << t_b;
	static constexpr uint8_t cyclic_shifts = t_b / util::gcd(t_b, 64 % t_b);
	static uint8_t offsets[cyclic_shifts];
	static uint128_t masks[t_v];
	constexpr static uint128_t even_mask = bm_rec<uint128_t>(bits::lo_set[t_b], t_b * 2, 128);
	constexpr static uint128_t carry_select_mask = bm_rec<uint128_t>(1ULL << t_b, t_b * 2, 128);

public:
	typedef typename rank_support_int<t_b>::size_type size_type;
	typedef typename rank_support_int<t_b>::value_type value_type;

	static void init()
	{
		for (value_type v = 0; v < t_v; ++v)
		{
			masks[v] = v;
			for (uint8_t i = t_b * 2; i < 128; i <<= 1)
				masks[v] |= masks[v] << i;
		}

		uint128_t tmp_carry = masks[1];
		for (value_type v = 0; v < t_v; ++v)
			masks[v] |= tmp_carry << t_b;

		for (uint8_t i = 1; i < cyclic_shifts; ++i)
			offsets[i] = t_b - ((64 - offsets[i - 1]) % t_b);
	}

    static uint64_t set_positions(uint64_t w, const value_type v)
    {
		// TODO: cast 128bit-bitmasks down to 64bit for efficiency reasons?
		uint64_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
    }

    static uint128_t set_positions(uint128_t w, const value_type v, const size_type word_pos)
    {
		w >>= 64 - (t_b - offsets[word_pos % cyclic_shifts]); // offset

		uint128_t res = (masks[v] - (even_mask & w)) & carry_select_mask;
		res |= ((masks[v] - (even_mask & (w >> t_b))) & carry_select_mask) << 1;
		return res;
    }

	static void set(uint128_t & w, const uint64_t hi, const uint64_t lo)
	{
		w = hi;
		w = (w <<= 64) | lo;
	}

	static uint32_t word_rank(const uint64_t* data, const size_type idx, const value_type v)
	{
		size_type bit_pos = idx * t_b;
		size_type word_pos = bit_pos >> 6;
		uint8_t word_shift = word_pos % cyclic_shifts;
		if (word_shift == 0)
		{
			uint64_t w = *(data + word_pos);
			return bits::cnt(set_positions(w, v) & bits::lo_set[(bit_pos & 0x3F) + 1]);
		}
		else
		{
			uint128_t w;
			set(w, *(data + word_pos), *(data + word_pos - 1));

			uint8_t in_word_pos = bit_pos & 0x3F;
			uint8_t w_width = in_word_pos + (t_b - offsets[word_shift]) + 1; // + 1 because carry bit will be one position further to the left

			uint128_t lo_set128;
			set(lo_set128, bits::lo_set[std::max((signed) 0, w_width - 64)], bits::lo_set[std::min((uint8_t)64, w_width)]);
			uint128_t res = set_positions(w, v, word_pos) & lo_set128;
			return bits::cnt((res >> 64) | res);
		}
	}

	// NOTE(cpockrandt): we only count full blocks!
	static uint32_t full_word_rank(const uint64_t* data, const size_type word_pos, const value_type v)
	{
		uint8_t word_shift = word_pos % cyclic_shifts;
		if (word_shift == 0)
		{
			uint64_t w = *(data + word_pos);
			return bits::cnt(set_positions(w, v));
		}
		else
		{
			uint128_t w;
			set(w, *(data + word_pos), *(data + word_pos - 1));

			uint8_t w_width = 64 + (t_b - offsets[word_shift]) + 1; // + 1 because carry bit will be one position further to the left
			if (word_shift != cyclic_shifts - 1) // (word_pos + 1) has offset
			 	w_width -= (t_b - offsets[word_shift + 1]);

			uint128_t lo_set128; // TODO(cpockrandt): construct global array
			set(lo_set128, bits::lo_set[std::max((signed) 0, w_width - 64)], bits::lo_set[std::min((uint8_t)64, w_width)]);
			uint128_t res = set_positions(w, v, word_pos) & lo_set128;
			return bits::cnt((res >> 64) | res);
		}
	}
};

template <uint8_t t_b>
uint64_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b == 0> >::masks[t_v];
template <uint8_t t_b>
uint128_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b != 0> >::masks[t_v];

template <uint8_t t_b>
uint8_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b != 0> >::offsets[rank_support_int_trait<t_b, std::enable_if_t<64 % t_b != 0> >::cyclic_shifts];

template <uint8_t t_b>
constexpr uint64_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b == 0> >::even_mask;
template <uint8_t t_b>
constexpr uint128_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b != 0> >::even_mask;

template <uint8_t t_b>
constexpr uint64_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b == 0> >::carry_select_mask;
template <uint8_t t_b>
constexpr uint128_t rank_support_int_trait<t_b, std::enable_if_t<64 % t_b != 0> >::carry_select_mask;

} // end namespace sdsl

#include "rank_support_int_v.hpp"
#include "rank_support_int_scan.hpp"

#endif // end file
