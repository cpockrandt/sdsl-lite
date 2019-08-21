// Copyright (c) 2016, the SDSL Project Authors.  All rights reserved.
// Please see the AUTHORS file for details.  Use of this source code is governed
// by a BSD license that can be found in the LICENSE file.
/*! \file wt_pc.hpp
    \brief wt_pc.hpp contains a class for the wavelet tree of byte sequences.
           The wavelet tree shape is parametrized by a prefix code.
    \author Simon Gog, Timo Beller
*/
#ifndef INCLUDED_SDSL_WT_EPR
#define INCLUDED_SDSL_WT_EPR

#include "bit_vectors.hpp"
#include "rank_support.hpp"
#include "rank_support_int.hpp"
#include "select_support.hpp"
#include "wt_helper.hpp"
#include <vector>
#include <utility>
#include <tuple>

//! Namespace for the succinct data structure library.
namespace sdsl {

//! A prefix code-shaped wavelet.
/*!
 * \tparam t_shape       Shape of the tree ().
 * \tparam t_bitvector   Underlying bitvector structure.
 * \tparam t_rank        Rank support for pattern `1` on the bitvector.
 * \tparam t_select      Select support for pattern `1` on the bitvector.
 * \tparam t_select_zero Select support for pattern `0` on the bitvector.
 * \tparam t_tree_strat  Tree strategy determines alphabet and the tree
 *                       class used to navigate the WT.
 *
 *  @ingroup wt
 */
template <uint8_t alphabet_size,
		  class bit_vector_type  = int_vector<>,
		  class rank_1_type      = rank_support_int_v<alphabet_size>,
		  class t_tree_strat     = byte_tree<>
		  >
class wt_pc_epr {
public:
	// typedef typename t_tree_strat::template type<wt_pc> tree_strat_type;
	typedef int_vector<>::size_type						size_type;
	typedef int_vector<>::value_type					value_type;
	// typedef typename t_bitvector::difference_type		difference_type;
	// typedef random_access_const_iterator<wt_pc>			const_iterator;
	// typedef const_iterator								iterator;
	// typedef int_vector<>								bit_vector_type;
	// typedef rank_support_int_v						rank_1_type;
	typedef wt_tag										index_category;
	typedef byte_alphabet_tag /*typename tree_strat_type::alphabet_category*/ alphabet_category;
	enum { lex_ordered = true };

private:

	size_type		m_size  = 0;  // original text size
	size_type		m_sigma = 0;  // alphabet size
	bit_vector_type m_bv;		  // bit vector to store the wavelet tree
	rank_1_type		m_bv_rank;	// rank support for the wavelet tree bit vector

	void construct_init_rank_select()
	{
		// util::init_support(m_bv_rank, &m_bv);
        rank_1_type temp(&m_bv);			 // generate a temporary support object
        m_bv_rank = std::move(temp); // swap its content with the target object
        m_bv_rank.set_vector(&m_bv);	 // set the support object's  pointer to x
	}

	// recursive internal version of the method interval_symbols
	// void _interval_symbols(size_type				i,
	// 					   size_type				j,
	// 					   size_type&				k,
	// 					   std::vector<value_type>& cs,
	// 					   std::vector<size_type>&  rank_c_i,
	// 					   std::vector<size_type>&  rank_c_j,
	// 					   node_type				v) const
	// {
	// 	// invariant: j>i
	// 	size_type i_new = (m_bv_rank(m_tree.bv_pos(v) + i) - m_tree.bv_pos_rank(v));
	// 	size_type j_new = (m_bv_rank(m_tree.bv_pos(v) + j) - m_tree.bv_pos_rank(v));
	// 	// goto left child
	// 	i -= i_new;
	// 	j -= j_new;
	// 	if (i != j) {
	// 		node_type v_new = m_tree.child(v, 0);
	// 		if (!m_tree.is_leaf(v_new)) {
	// 			_interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, v_new);
	// 		} else {
	// 			rank_c_i[k] = i;
	// 			rank_c_j[k] = j;
	// 			cs[k++]		= m_tree.bv_pos_rank(v_new);
	// 		}
	// 	}
	// 	// goto right child
	// 	if (i_new != j_new) {
	// 		node_type v_new = m_tree.child(v, 1);
	// 		if (!m_tree.is_leaf(v_new)) {
	// 			_interval_symbols(i_new, j_new, k, cs, rank_c_i, rank_c_j, v_new);
	// 		} else {
	// 			rank_c_i[k] = i_new;
	// 			rank_c_j[k] = j_new;
	// 			cs[k++]		= m_tree.bv_pos_rank(v_new);
	// 		}
	// 	}
	// }

public:
	const size_type&	   sigma = m_sigma;
	const bit_vector_type& bv	= m_bv;

	// Default constructor
	wt_pc_epr(){};

	//! Construct the wavelet tree from a sequence defined by two interators
	/*!
         * \param begin Iterator to the start of the input.
         * \param end   Iterator one past the end of the input.
         * \par Time complexity
         *      \f$ \Order{n\log|\Sigma|}\f$, where \f$n=size\f$
         */
	template <typename t_it>
	wt_pc_epr(t_it begin, t_it end) : m_size(std::distance(begin, end))
	{
		if (0 == m_size) return;
		// O(n + |\Sigma|\log|\Sigma|) algorithm for calculating node sizes
		// TODO: C should also depend on the tree_strategy. C is just a mapping
		// from a symbol to its frequency. So a map<uint64_t,uint64_t> could be
		// used for integer alphabets...
		std::vector<size_type> C;
		// 1. Count occurrences of characters
		calculate_character_occurences(begin, end, C);
		// 2. Calculate effective alphabet size
		calculate_effective_alphabet_size(C, m_sigma);
		// 4. Generate wavelet tree bit sequence m_bv

		int_vector<> m_bv_tmp;
		m_bv_tmp.width(std::ceil(std::log2(m_sigma)));
		m_bv_tmp.resize(m_size);

		// TODO: use byte_alphabet to reduce alphabet to [0..sigma-1]? but dont use by default, since inefficient!
		// TODO: was wenn ein buchstaben gesucht wird, der nicht vorkommt?
		size_type idx = 0;
		for (t_it it = begin; it != end; ++it) {
			m_bv_tmp[idx++] = *it;
		}
		bit_vector_type m_bv2(m_bv_tmp);
		m_bv = std::move(m_bv2); // swap its content with the target object
		// 5. Initialize rank and select data structures for m_bv
		construct_init_rank_select();
	}

	template <typename t_it>
	wt_pc_epr(t_it begin, t_it end, std::string) : wt_pc_epr(begin, end)
	{
	}

	//! Copy constructor
	wt_pc_epr(const wt_pc_epr& wt)
		: m_size(wt.m_size)
		, m_sigma(wt.m_sigma)
		, m_bv(wt.m_bv)
		, m_bv_rank(wt.m_bv_rank)
	{
		m_bv_rank.set_vector(&m_bv);
	}

	wt_pc_epr(wt_pc_epr&& wt)
		: m_size(wt.m_size)
		, m_sigma(wt.m_sigma)
		, m_bv(std::move(wt.m_bv))
		, m_bv_rank(std::move(wt.m_bv_rank))
	{
		m_bv_rank.set_vector(&m_bv);
	}

	//! Assignment operator
	wt_pc_epr& operator=(const wt_pc_epr& wt)
	{
		if (this != &wt) {
			wt_pc_epr tmp(wt);			// re-use copy-constructor
			*this = std::move(tmp); // re-use move-assignment
		}
		return *this;
	}

	//! Move assignment operator
	wt_pc_epr& operator=(wt_pc_epr&& wt)
	{
		if (this != &wt) {
			m_size	= wt.m_size;
			m_sigma   = wt.m_sigma;
			m_bv	  = std::move(wt.m_bv);
			m_bv_rank = std::move(wt.m_bv_rank);
			m_bv_rank.set_vector(&m_bv);
		}
		return *this;
	}

	//! Returns the size of the original vector.
	size_type size() const { return m_size; }

	//! Returns whether the wavelet tree contains no data.
	bool empty() const { return m_size == 0; }

	//! Recovers the i-th symbol of the original vector.
	/*!
         * \param i Index in the original vector.
         * \return The i-th symbol of the original vector.
         * \par Time complexity
         *      \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *      zero order entropy of the sequence
         *
         * \par Precondition
         *      \f$ i < size() \f$
         */
	value_type operator[](size_type i) const
	{
		assert(i < size());
		return m_bv[i];
	};

	//! Calculates how many symbols c are in the prefix [0..i-1].
	/*!
         * \param i Exclusive right bound of the range.
         * \param c Symbol c.
         * \return Number of occurrences of symbol c in the prefix [0..i-1].
         * \par Time complexity
         *      \f$ \Order{H_0} \f$ on average, where \f$ H_0 \f$ is the
         *      zero order entropy of the sequence
         *
         * \par Precondition
         *      \f$ i \leq size() \f$
         */

	size_type rank(size_type i, value_type c) const
	{
		assert(i <= size());
		return m_bv_rank.rank(i, c);
	};

	//! Calculates how many times symbol wt[i] occurs in the prefix [0..i-1].
	/*!
         * \param i The index of the symbol.
         * \return  Pair (rank(wt[i],i),wt[i])
         * \par Time complexity
         *      \f$ \Order{H_0} \f$
         *
         * \par Precondition
         *      \f$ i < size() \f$
         */
	std::pair<size_type, value_type> inverse_select(size_type i) const
	{
		assert(i < size());
		return std::make_pair(m_bv_rank.rank(i, m_bv[i]), m_bv[i]);
	}


	//! For each symbol c in wt[i..j-1] get rank(i,c) and rank(j,c).
	/*!
         * \param i        The start index (inclusive) of the interval.
         * \param j        The end index (exclusive) of the interval.
         * \param k        Reference for number of different symbols in [i..j-1].
         * \param cs       Reference to a vector that will contain in
         *                 cs[0..k-1] all symbols that occur in [i..j-1] in
         *                 arbitrary order (if lex_ordered = false) and ascending
         *                 order (if lex_ordered = true).
         * \param rank_c_i Reference to a vector which equals
         *                 rank_c_i[p] = rank(i,cs[p]), for \f$ 0 \leq p < k \f$.
         * \param rank_c_j Reference to a vector which equals
         *                 rank_c_j[p] = rank(j,cs[p]), for \f$ 0 \leq p < k \f$.
         * \par Time complexity
         *      \f$ \Order{\min{\sigma, k \log \sigma}} \f$
         *
         * \par Precondition
         *      \f$ i \leq j \leq size() \f$
         *      \f$ cs.size() \geq \sigma \f$
         *      \f$ rank_{c_i}.size() \geq \sigma \f$
         *      \f$ rank_{c_j}.size() \geq \sigma \f$
         */
	// void interval_symbols(size_type				   i,
	// 					  size_type				   j,
	// 					  size_type&			   k,
	// 					  std::vector<value_type>& cs,
	// 					  std::vector<size_type>&  rank_c_i,
	// 					  std::vector<size_type>&  rank_c_j) const
	// {
	// 	assert(i <= j and j <= size());
	// 	if (i == j) {
	// 		k = 0;
	// 	} else if (1 == m_sigma) {
	// 		k			= 1;
	// 		cs[0]		= m_tree.bv_pos_rank(m_tree.root());
	// 		rank_c_i[0] = std::min(i, m_size);
	// 		rank_c_j[0] = std::min(j, m_size);
	// 	} else if ((j - i) == 1) {
	// 		k			= 1;
	// 		auto rc		= inverse_select(i);
	// 		rank_c_i[0] = rc.first;
	// 		cs[0]		= rc.second;
	// 		rank_c_j[0] = rank_c_i[0] + 1;
	// 	} else if ((j - i) == 2) {
	// 		auto rc		= inverse_select(i);
	// 		rank_c_i[0] = rc.first;
	// 		cs[0]		= rc.second;
	// 		rc			= inverse_select(i + 1);
	// 		rank_c_i[1] = rc.first;
	// 		cs[1]		= rc.second;
	//
	// 		if (cs[0] == cs[1]) {
	// 			k			= 1;
	// 			rank_c_j[0] = rank_c_i[0] + 2;
	// 		} else {
	// 			k = 2;
	// 			if (lex_ordered and cs[0] > cs[1]) {
	// 				std::swap(cs[0], cs[1]);
	// 				std::swap(rank_c_i[0], rank_c_i[1]);
	// 			}
	// 			rank_c_j[0] = rank_c_i[0] + 1;
	// 			rank_c_j[1] = rank_c_i[1] + 1;
	// 		}
	// 	} else {
	// 		k = 0;
	// 		_interval_symbols(i, j, k, cs, rank_c_i, rank_c_j, 0);
	// 	}
	// }


	//! How many symbols are lexicographic smaller/greater than c in [i..j-1].
	/*!
         * \param i       Start index (inclusive) of the interval.
         * \param j       End index (exclusive) of the interval.
         * \param c       Symbol c.
         * \return A triple containing:
         *         * rank(i,c)
         *         * #symbols smaller than c in [i..j-1]
         *         * #symbols greater than c in [i..j-1]
         *
         * \par Precondition
         *       \f$ i \leq j \leq size() \f$
         * \note
         * This method is only available if lex_ordered = true
         */
	template <class t_ret_type = std::tuple<size_type, size_type, size_type>>
	t_ret_type lex_count(size_type i, size_type j, value_type c) const
	{
		assert(i <= j and j <= size());
		// size_type smaller = 0;
		// size_type greater = (j - i) - (m_bv_rank.prefix_rank(j, c) - m_bv_rank.prefix_rank(i, c));
		// if (c > 0)
		// 	smaller = m_bv_rank.prefix_rank(j, c-1) - m_bv_rank.prefix_rank(i, c-1);
		// size_type rank = m_bv_rank.rank(i, c);

		size_type smaller = 0;
		size_type prefix_i_c = m_bv_rank.prefix_rank(i, c);
		size_type prefix_i_c_1 = 0;
		size_type greater = j - i - m_bv_rank.prefix_rank(j, c) + prefix_i_c;
		if (c > 0)
		{
			prefix_i_c_1 = m_bv_rank.prefix_rank(i, c-1);
			smaller = m_bv_rank.prefix_rank(j, c-1) - prefix_i_c_1;
		}
		size_type rank = prefix_i_c - prefix_i_c_1;

		return t_ret_type{rank, smaller, greater};
	}

	//! How many symbols are lexicographic smaller than c in [0..i-1].
	/*!
         * \param i Exclusive right bound of the range.
         * \param c Symbol c.
         * \return A tuple containing:
         *         * rank(i,c)
         *         * #symbols smaller than c in [0..i-1]
         * \par Precondition
         *       \f$ i \leq size() \f$
         * \note
         * This method is only available if lex_ordered = true
         */
	template <class t_ret_type = std::tuple<size_type, size_type>>
	t_ret_type lex_smaller_count(size_type i, value_type c) const
	{
		assert(i <= size());
		size_type prefix_c_1 = m_bv_rank.prefix_rank(i, c - 1);
		return t_ret_type{m_bv_rank.prefix_rank(i, c) - prefix_c_1, prefix_c_1};
	}

	//! Returns a const_iterator to the first element.
	// const_iterator begin() const { return const_iterator(this, 0); }
	//
	// //! Returns a const_iterator to the element after the last element.
	// const_iterator end() const { return const_iterator(this, size()); }

	//! Serializes the data structure into the given ostream
	size_type
	serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
	{
		structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
		size_type			 written_bytes = 0;
		written_bytes += write_member(m_size, out, child, "size");
		written_bytes += write_member(m_sigma, out, child, "sigma");
		written_bytes += m_bv.serialize(out, child, "bv");
		written_bytes += m_bv_rank.serialize(out, child, "bv_rank");
		structure_tree::add_size(child, written_bytes);
		return written_bytes;
	}

	//! Loads the data structure from the given istream.
	void load(std::istream& in)
	{
		read_member(m_size, in);
		read_member(m_sigma, in);
		m_bv.load(in);
		m_bv_rank.load(in, &m_bv);
	}

	//! Returns for a symbol c the next larger or equal symbol in the WT.
	/*! \param c the symbol
         *  \return A pair. The first element of the pair consititues if
         *          a valid answer was found (true) or no valid answer (false)
         *          could be found. The second element contains the found symbol.
         */
	// std::pair<bool, value_type> symbol_gte(value_type c) const { return m_tree.symbol_gte(c); }

	//! Returns for a symbol c the previous smaller or equal symbol in the WT.
	/*! \param c the symbol
         *  \return A pair. The first element of the pair consititues if
         *          a valid answer was found (true) or no valid answer (false)
         *          could be found. The second element contains the found symbol.
         */
	// std::pair<bool, value_type> symbol_lte(value_type c) const { return m_tree.symbol_lte(c); }

};
}

#endif
