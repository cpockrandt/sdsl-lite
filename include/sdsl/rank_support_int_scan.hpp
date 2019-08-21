// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file rank_support_int_scan.hpp
    \brief rank_support_int_scan.hpp contains rank_support_int_scan that support a sdsl::bit_vector with linear time rank information.
    \author Simon Gog
*/
#ifndef INCLUDED_SDSL_RANK_SUPPORT_INT_SCAN
#define INCLUDED_SDSL_RANK_SUPPORT_INT_SCAN

#include "rank_support_int.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl {

//! A class supporting rank queries in linear time.
/*! \par Space complexity
 *       Constant.
 *  \par Time complexity
 *       Linear in the size of the supported vector.
 *
 *  \tparam t_b       Bit pattern which should be supported. Either `0`,`1`,`10`,`01`.
 *  \tparam t_pat_len Length of the bit pattern.
 * @ingroup rank_support_group
 */

template <uint8_t alphabet_size>
class rank_support_int_scan : public rank_support_int<alphabet_size> {
public:
	typedef int_vector<> int_vector_type;
    typedef typename rank_support_int<alphabet_size>::size_type size_type;
    typedef typename rank_support_int<alphabet_size>::value_type value_type;
	// enum { bit_pat = t_b };
	// enum { bit_pat_len = t_pat_len };
private:
	// using rank_support_int::m_v;
	// using rank_support_int::t_b;
	// using rank_support_int::t_v;

public:
	explicit rank_support_int_scan(const int_vector<>* v = nullptr) : rank_support_int<alphabet_size>(v){
    	// rank_support_int_trait<>::init();
	};
	rank_support_int_scan(const rank_support_int_scan& rs) = default;
	rank_support_int_scan(rank_support_int_scan&& rs)	  = default;
	rank_support_int_scan& operator=(const rank_support_int_scan& rs) = default;
	rank_support_int_scan& operator=(rank_support_int_scan&& rs) = default;
	size_type rank(size_type idx, const value_type v) const;
	size_type operator()(size_type idx, const value_type v) const { return rank(idx, v); };
	size_type prefix_rank(size_type idx, const value_type v) const;
	size_type size() const { return this->m_v->size(); };
	size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		return serialize_empty_object(out, v, name, this);
	}
	void load(std::istream&, const int_vector<>* v = nullptr) { this->m_v = v; this->init(v); }
	void set_vector(const int_vector<>* v = nullptr) { this->m_v = v; }

	void print() const {
		std::cout << "t_b = " << (unsigned)this->t_b << "\n";
		std::cout << "t_v = " << (unsigned)this->t_v << "\n";
		std::cout << std::bitset<64>(this->even_mask) << "\n";
		std::cout << std::bitset<64>(this->carry_select_mask) << "\n";
		for (unsigned i = 0; i < this->masks.size(); ++i)
			std::cout << std::bitset<64>(this->masks[i]) << "\n";
	}
};

template <uint8_t alphabet_size>
inline typename rank_support_int_scan<alphabet_size>::size_type
rank_support_int_scan<alphabet_size>::rank(size_type idx, const value_type v) const
{
	assert(v < this->t_v);
	assert(this->m_v != nullptr);
	assert(idx <= this->m_v->size());

	if (unlikely(v == 0))
		return prefix_rank(idx, v);

	const uint64_t* p = this->m_v->data();
	size_type i = 0;
	size_type result = 0;
	size_type word_pos = (idx * this->t_b) >> 6;
	while (i < word_pos) {
		result += this->full_word_rank(p, i, v);
				// - rank_support_int_trait<t_b>::full_word_prefix_rank(p, i, v - 1);
		++i;
	}
	// result += rank_support_int_trait<t_b>::word_prefix_rank(p, idx, v);
			// - rank_support_int_trait<t_b>::word_prefix_rank(p, idx, v - 1);
	return result + word_rank(p, idx, v);
}

template <uint8_t alphabet_size>
inline typename rank_support_int_scan<alphabet_size>::size_type
rank_support_int_scan<alphabet_size>::prefix_rank(size_type idx, const value_type v) const
{
	assert(v < this->t_v);
	assert(this->m_v != nullptr);
	assert(idx <= this->m_v->size());

	// std::cout << "t_b = " << (unsigned) t_b << ", t_v = " << (unsigned) t_v << "\n";
	if (unlikely(v == this->t_v - 1))
		return idx;

	const uint64_t* p	    = this->m_v->data();
	size_type		word_pos = (idx * this->t_b) >> 6;
	size_type       i       = 0;
	size_type		result  = 0;
	// std::cout << "full_word_prefix_rank(" << p << ", " << i << ", " << (unsigned)v << "): ";
	while (i < word_pos) {
		result += this->full_word_prefix_rank(p, i, v);
		// std::cout << full_word_prefix_rank(p, i, v) << " + ";
		++i;
	}
	// std::cout << word_prefix_rank(p, idx, v) << "\n";
	return result + this->word_prefix_rank(p, idx, v);
}

} // end namespace sds

#endif // end file
