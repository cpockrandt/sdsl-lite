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
template <uint8_t t_b>
class rank_support_int_scan : public rank_support_int<t_b> {
public:
	typedef int_vector<t_b> int_vector_type;
    typedef typename rank_support_int<t_b>::size_type size_type;
    typedef typename rank_support_int<t_b>::value_type value_type;
	// enum { bit_pat = t_b };
	// enum { bit_pat_len = t_pat_len };
private:
	using rank_support_int<t_b>::m_v;
    constexpr static uint8_t t_v = 1ULL << t_b;

public:
	explicit rank_support_int_scan(const int_vector<t_b>* v = nullptr) : rank_support_int<t_b>(v){};
	rank_support_int_scan(const rank_support_int_scan& rs) = default;
	rank_support_int_scan(rank_support_int_scan&& rs)	  = default;
	rank_support_int_scan& operator=(const rank_support_int_scan& rs) = default;
	rank_support_int_scan& operator=(rank_support_int_scan&& rs) = default;
	size_type rank(size_type idx, const value_type v) const;
	size_type operator()(size_type idx, const value_type v) const { return rank(idx, v); };
	size_type					   size() const { return m_v->size(); };
	size_type
	serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		return serialize_empty_object(out, v, name, this);
	}
	void load(std::istream&, const int_vector<t_b>* v = nullptr) { m_v = v; }
};

template <uint8_t t_b>
inline typename rank_support_int_scan<t_b>::size_type
rank_support_int_scan<t_b>::rank(size_type idx, const value_type v) const
{
	assert(m_v != nullptr);
	assert(idx <= m_v->size());
	const uint64_t* p	    = m_v->data();
	size_type		word_pos = (idx * t_b) >> 6;
	size_type       i       = 0;
	size_type		result  = 0;
	while (i < word_pos) {
		result += rank_support_int_trait<t_b>::full_word_rank(p, i, v);
		++i;
	}
	return result + rank_support_int_trait<t_b>::word_rank(p, idx, v);
}

} // end namespace sds

#endif // end file
