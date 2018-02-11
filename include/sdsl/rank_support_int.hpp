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
	virtual size_type rank(size_type i) const = 0;
	//! Alias for rank(i)
	virtual size_type operator()(size_type idx) const = 0;
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

constexpr uint64_t create_bitmask_recursive(unsigned const bitsUsed, unsigned const blocksLeft, unsigned const blocksProcessed, unsigned const blocksize, TWord const value1, TWord const value2, TWord const result)
{
    return (blocksLeft == 0)
        ? result
        : create_bitmask_recursive(bitsUsed,
                                   blocksLeft - 1,
                                   blocksProcessed + 1,
                                   blocksize,
                                   value1,
                                   value2,
                                   result | (((blocksProcessed % 2) ? value1 : value2)
                                           << (bitsUsed - blocksProcessed * blocksize)));
}

template <uint8_t t_bits, uint8_t t_value>
constexpr uint64_t create_bitmask()
{
    return create_bitmask_recursive(0);
}

template <uint8_t t_bits, uint8_t t_value>
struct rank_support_int_trait {
protected:
    static constexpr uint64_t bitmask = create_bitmask<t_bits, t_value>();
    // TODO: eine bitmask für t_v. build in constructor
public:
	typedef typename rank_support_int<t_bits, t_value>::size_type size_type;

    // TWord const erg1 = (TRD::_CHAR_BITMASKS[ordValue(c)] - (word & TRD::_SELECT_BITMASK)) & TRD::_TRUNC_BITMASKS[(posInWord + 1)/2];
    // TWord const erg2 = (TRD::_CHAR_BITMASKS[ordValue(c)] - ((word >> TRD::_BITS_PER_VALUE) & TRD::_SELECT_BITMASK)) & TRD::_TRUNC_BITMASKS[(posInWord/2) + 1];
    // return popCount(erg1 | (erg2 << 1));

	static uint32_t word_rank(const uint64_t* data, size_type idx)
	{
        std::cout << "<" << (unsigned) t_bits << ", " << (unsigned) t_value << ">: " << std::bitset<64>(bitmask) << '\n';
		return bits::cnt((~*(data + (idx >> 6))) & bits::lo_set[idx & 0x3F]);
	}

	static uint32_t full_word_rank(const uint64_t* data, size_type idx)
	{
		return bits::cnt((~*(data + (idx >> 6))));
	}
};

// template <uint8_t t_b, uint8_t t_v>
// constexpr uint64_t rank_support_int_trait<t_b, t_v>::bitmask;

} // end namespace sdsl

// #include "rank_support_int_v.hpp"
#include "rank_support_int_scan.hpp"

#endif // end file
