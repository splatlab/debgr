/*
 * =====================================================================================
 *
 *       Filename:  hashutil.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2016 04:49:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#include<set>
#include<map>
#include<utility>
#include<unordered_set>
#include<iostream>
#include<cassert>
#include <fstream>
#include <time.h>
#include <sys/time.h>

#include "debruijn_graph.h"
#include "hashutil.h"

static void print_time_elapsed(std::string desc,
															 struct timeval* start,
															 struct timeval* end);

#define BITMASK(nbits) ((nbits) == 64								\
												? 0xffffffffffffffff				\
												: (1ULL << (nbits)) - 1ULL)

using namespace std;

namespace dna {

	/////////////// bases /////////////////
	base operator-(base b) {
		return (base)((~((uint64_t)b)) & 0x3ULL);
	}
	const base bases[4] = {C, A, T, G};
	const map<char, base> base_from_char = { {'A', A}, {'C', C},
		{'G', G}, {'T', T},
		{'N', A}};
	const map<base, char> base_to_char = { {A, 'A'}, {C, 'C'},
		{G, 'G'}, {T, 'T'}};

	///////////// kmers /////////////////////
	kmer::kmer(void) : len(0), val(0) {}
	kmer::kmer(base b) : len(1), val((uint64_t)b) {}
	kmer::kmer(int l, uint64_t v) : len(l), val(v & BITMASK(2*l)) {
		assert(l <= 32);
	}

	static uint64_t string_to_kmer_val(string s) {
		uint64_t val = 0;
		for (auto c : s)
			val = (val << 2) | ((uint64_t)(base_from_char.at(c)));
		return val;
	}

	kmer::kmer(string s) : len(s.size()), val(string_to_kmer_val(s)) {
		assert(s.size() <= 32);
	}

	// Convert to string
	kmer::operator string() const {
		string s;
		for (auto i = 1; i < len+1; i++)
			s = s + base_to_char.at((base)((val >> (2*(len - i))) & BITMASK(2)));
		return s;
	}

	bool operator<(kmer a, kmer b) {
		return a.len != b.len ? a.len < b.len : a.val < b.val;
	}
	bool operator==(kmer a, kmer b) {
		return a.len == b.len && a.val == b.val;
	}
	bool operator!=(kmer a, kmer b) {
		return !operator==(a, b);
	}

	// Return the reverse complement of k
	kmer operator-(kmer k) {
		uint64_t val = k.val;
		val =
			(val >> 32) |
			(val << 32);
		val =
			((val >> 16) & 0x0000ffff0000ffff) |
			((val << 16) & 0xffff0000ffff0000);
		val =
			((val >> 8) & 0x00ff00ff00ff00ff) |
			((val << 8) & 0xff00ff00ff00ff00);
		val =
			((val >> 4) & 0x0f0f0f0f0f0f0f0f) |
			((val << 4) & 0xf0f0f0f0f0f0f0f0);
		val =
			((val >> 2) & 0x3333333333333333) |
			((val << 2) & 0xcccccccccccccccc);
		val = ~val;
		val >>= 64-2*k.len;
		return kmer(k.len, val);
	}

	// backwards from standard definition to match kmer.h definition
	kmer canonicalize(kmer k) {
		return -k < k ? k : -k;
	}

	// Return the kmer of length |a| that results from shifting b into a
	// from the right
	kmer operator<<(kmer a, kmer b) {
		uint64_t val = ((a.val << (2*b.len)) | b.val) & BITMASK(2*a.len);
		return kmer(a.len, val);
	}
	// Return the kmer of length |b| that results from shifting b into a
	// from the left
	kmer operator>>(kmer a, kmer b) {
		uint64_t val
			= ((b.val >> (2*a.len)) | (a.val << (2*(b.len - a.len))))
			& BITMASK(2*b.len);
		return kmer(b.len, val);
	}
	// Append two kmers
	kmer operator+(kmer a, kmer b) {
		int len = a.len + b.len;
		assert(len <= 32);
		uint64_t val = (a.val << (2*b.len)) | b.val;
		return kmer(len, val);
	}

	kmer prefix(kmer k, int len) { return kmer(len, k.val >> (2*(k.len - len))); }
	kmer suffix(kmer k, int len) { return kmer(len, k.val & BITMASK(2*len)) ; }

	bool period_divides(kmer k, uint64_t periodicity)	{
		static const uint64_t multipliers[33] =
			{
				0,
				0x5555555555555555, // 1
				0x1111111111111111, // 2
				0x1041041041041041, // 3
				0x0101010101010101, // 4
				0x1004010040100401, // 5
				0x1001001001001001, // 6
				0x0100040010004001, // 7
				0x0001000100010001, // 8
				0x0040001000040001, // 9
				0x1000010000100001, // 10
				0x0000100000400001, // 11
				0x0001000001000001, // 12
				0x0010000004000001, // 13
				0x0100000010000001, // 14
				0x1000000040000001, // 15
				0x0000000100000001, // 16
				0x0000000400000001, // 17
				0x0000001000000001, // 18
				0x0000004000000001, // 19
				0x0000010000000001, // 20
				0x0000040000000001, // 21
				0x0000100000000001, // 22
				0x0000400000000001, // 23
				0x0001000000000001, // 24
				0x0004000000000001, // 25
				0x0010000000000001, // 26
				0x0040000000000001, // 27
				0x0100000000000001, // 28
				0x0400000000000001, // 29
				0x1000000000000001, // 30
				0x4000000000000001, // 31
				0x0000000000000001, // 32
			};
		uint64_t piece = k.val & BITMASK(2*periodicity);
		piece = piece * multipliers[periodicity];
		piece = piece & BITMASK(2*k.len);
		return piece == k.val;
	}
	
	uint64_t period(kmer k) {
		for (int i = 1; i <= k.len; i++) {
			if (period_divides(k, i))
				return i;
		}
		abort();
	}
	
	canonical_kmer::canonical_kmer(void) : kmer() {}
	canonical_kmer::canonical_kmer(base b) : kmer(canonicalize(kmer(b))) {}
	canonical_kmer::canonical_kmer(int l, uint64_t v)
		: kmer(canonicalize(kmer(l, v))) {}
	canonical_kmer::canonical_kmer(std::string s) : kmer(canonicalize(kmer(s))) {}
	canonical_kmer::canonical_kmer(kmer k) : kmer(canonicalize(k)) {}
}

namespace std {
	size_t hash<dna::kmer>::operator()(const dna::kmer & x) const {
		return hash<int>()(x.val);
	}
	size_t
	hash<dna::canonical_kmer>::operator()(const dna::canonical_kmer & x) const {
		return hash<int>()(x.val);
	}
}

/////////////// Interface for CQFs of kmers ///////////////
namespace dna {
	template<class kmer_type>
	kmer_counter_cqf<kmer_type>::kmer_counter_cqf(uint64_t l, QF cqf)
		: counts(cqf), len(l) {}
	
	template<class kmer_type>
	static uint64_t kmer_hash(kmer_type k, uint32_t seed, __uint128_t range) {
		return
			kmercounting::HashUtil::MurmurHash64A(&k.val,
																						sizeof(k.val),
																						seed)
			% range;
	}

	template<class kmer_type>
	uint64_t kmer_counter_cqf<kmer_type>::get_hash(kmer_type k) const {
		return kmer_hash(k, counts.metadata->seed, counts.metadata->range);
	}

	template<class kmer_type>
	void kmer_counter_cqf<kmer_type>::increment(kmer_type k, uint64_t amount) {
		qf_insert(&counts, get_hash(k), 0, amount, 0, 0);
	}

	template<class kmer_type>
	size_t kmer_counter_cqf<kmer_type>::count(kmer_type k) const {
		return qf_count_key_value(&counts, get_hash(k), 0);
	}

	template<class kmer_type>
	kmer_counter_lossless_cqf<kmer_type>::kmer_counter_lossless_cqf(uint64_t l,
																																	QF cqf)
		: kmer_counter_cqf<kmer_type>(l, cqf) {}

	template<class kmer_type>
	uint64_t kmer_invertible_hash(kmer_type k, uint32_t seed) {
		return kmercounting::HashUtil::hash_64(k.val, BITMASK(2*k.len));
	}

	template<class kmer_type>
	void kmer_counter_lossless_cqf<kmer_type>::increment(kmer_type k,
																											 uint64_t amount) {
		uint64_t h = kmer_invertible_hash<kmer_type>(k, this->counts.metadata->seed);
		qf_insert(&this->counts, h, 0, amount, 0, 0);
	}

	template<class kmer_type>
	size_t kmer_counter_lossless_cqf<kmer_type>::count(kmer_type k) const {
		uint64_t h = kmer_invertible_hash<kmer_type>(k, this->counts.metadata->seed);
		return qf_count_key_value(&this->counts, h, 0);
	}

	template<class kmer_type>
	uint64_t kmer_counter_lossless_cqf<kmer_type>::size(void) const {
		return this->counts.metadata->ndistinct_elts;
	}
	
	template<class kmer_type>
	kmer_counter_lossless_cqf<kmer_type>::iterator::iterator(uint64_t l, QFi it)
		: len(l), iter(it) {};

	template<class kmer_type>
	kmer_type
	kmer_counter_lossless_cqf<kmer_type>::iterator::operator*(void) const {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&iter, &key, &value, &count);
		return kmer_type(len, kmercounting::HashUtil::hash_64i(key, BITMASK(2*len)));
	}

	template<class kmer_type>
		void kmer_counter_lossless_cqf<kmer_type>::iterator::operator++(void) {
			qfi_next(&iter);
		}

	template<class kmer_type>
	bool
	operator!=(const typename kmer_counter_lossless_cqf<kmer_type>::iterator &a,
						 const typename kmer_counter_lossless_cqf<kmer_type>::iterator &b) {
		return !qfi_end(&a.iter) || !qfi_end(&b.iter);
	}

	template<class kmer_type>
	typename kmer_counter_lossless_cqf<kmer_type>::iterator
	kmer_counter_lossless_cqf<kmer_type>::begin(void) const {
		QFi qfi;
		qf_iterator(&this->counts, &qfi, 0);
		return iterator(this->len, qfi);
	}

	template<class kmer_type>
	typename kmer_counter_lossless_cqf<kmer_type>::iterator
	kmer_counter_lossless_cqf<kmer_type>::end(void) const {
		QFi qfi;
		qf_iterator(&this->counts, &qfi, 0xffffffffffffffff);
		return iterator(this->len, qfi);
	}

	/////////////// CQF-based de Bruijn graphs ///////////////

	debruijn_graph::debruijn_graph(int l,
																 kmer_counter_cqf<canonical_kmer> occ,
																 kmer_counter_lossless_cqf<canonical_kmer> begs,
																 kmer_counter_lossless_cqf<canonical_kmer> enz,
																 map<edge, uint64_t> cor
																)
		: len(l), occurences(occ), begins(begs), ends(enz), corrections(cor) {}

	uint64_t debruijn_graph::abundance(edge e) const {
		uint64_t o = get_occurences(e);
		if (o)
			return o - get_corrections(e);
		return 0;
	}

	void debruijn_graph::insert_read(string s) {
		// TODO: handle locking
		if (s.size() < len)
			return;
		kmer k(s.substr(0, len));
		node n = prefix(k, len-1);    
		record_end(n, k);
		record_occurence(k);
		for (auto it = s.begin() + len; it != s.end(); ++it) {
			node nextn = suffix(k, len-1);
			if (can_have_duplex_edges(nextn))
				record_end(nextn, k);
			k = k << base_from_char.at(*it);
			record_occurence(k);
			if (can_have_duplex_edges(nextn))
				record_end(nextn, k);
		}
		n = suffix(k, len-1);
		record_end(n, k);
	}

	const uint64_t debruijn_graph::infinite_depth = ((uint64_t)-1);

	template<template<class, class> class Q> 
	typename debruijn_graph::node_iterator<Q>::work_item
	debruijn_graph::node_iterator<Q>::front(const queue<work_item> &w) const {
		return w.front();
	}

	template<template<class, class> class Q> 
	typename debruijn_graph::node_iterator<Q>::work_item
	debruijn_graph::node_iterator<Q>::front(const stack<work_item> &w) const {
		return w.top();
	}

	template<template<class, class> class Q> 
	debruijn_graph::node_iterator<Q>::node_iterator(const debruijn_graph &g,
																									const std::unordered_set<node> starts,
																									uint64_t depth)
		: dbg(g) {
		for (const auto &n : starts) {
			work_item w;
			w.curr = n;
			w.depth = depth;
			work.push(w);
			visited.insert(n);
		}
	}

	template<template<class, class> class Q>
	debruijn_graph::node debruijn_graph::node_iterator<Q>::operator*(void) const {
		return front(work).curr;
	}

	template<template<class, class> class Q>
	void debruijn_graph::node_iterator<Q>::operator++(void) {
		work_item w;
		w = front(work);
		work.pop();
		
		uint64_t nsiblings = 0;
		if (w.depth)
			for (const auto &neighbor : dbg.neighbors(w.curr))
				if (neighbor != w.curr && neighbor != w.prev) {
					if (visited.count(neighbor))
						nsiblings++;
					else {
						work_item neww = { neighbor,
															 w.curr,
															 w.depth == infinite_depth
															 ? infinite_depth
															 : w.depth - 1 };
						visited.insert(neighbor);
						work.push(neww);
					}
				}
		if (w.depth == infinite_depth && nsiblings == 0)
			visited.erase(w.curr);
	}

	template<template<class, class> class Q> bool
	debruijn_graph::node_iterator<Q>::operator!=(const node_iterator<Q> &other) const {
		return !work.empty() || !other.work.empty();
	}

	template<template<class, class> class Q>
	debruijn_graph::nodes_container<Q>::nodes_container(const debruijn_graph &g,
																											const std::unordered_set<node> strts,
																											uint64_t dpth)
		: dbg(g), starts(strts), depth(dpth) {}

	template<template<class, class> class Q> 
	debruijn_graph::node_iterator<Q>
	debruijn_graph::nodes_container<Q>::begin(void) const {
		return node_iterator<Q>(dbg, starts, depth);
	}

	template<template<class, class> class Q>
	debruijn_graph::node_iterator<Q>
	debruijn_graph::nodes_container<Q>::end(void) const {
		return node_iterator<Q>(dbg, unordered_set<node>(), depth);
	}

	debruijn_graph::nodes_container<stack>
	debruijn_graph::reachable_nodes(void) const {
		return nodes_container<stack>(*this, get_all_ends(),
																	infinite_depth);
	}
	
	debruijn_graph::nodes_container<stack>
	debruijn_graph::connected_component(node n) const {
		return nodes_container<stack>(*this,
																	unordered_set<node>({n}),
																	infinite_depth);
	}

	debruijn_graph::nodes_container<queue>
	debruijn_graph::nearby_nodes(node n, uint64_t depth) const {
		return nodes_container<queue>(*this, unordered_set<node>({n}), depth);
	}

	template<template<class, class> class Q> 
	debruijn_graph::edge_iterator<Q>::edge_iterator(const debruijn_graph &g,
																									const unordered_set<node> strts,
																									uint64_t dpth)
		: dbg(g),
		anode(g, strts, dpth),
		nodeend(dbg, unordered_set<node>(), infinite_depth),
		anode_edges(),
		aedge(anode_edges.end())
	{
		while (anode != nodeend) {
			anode_edges = dbg.edges(*anode);
			aedge = anode_edges.begin();
			while (aedge != anode_edges.end() &&
						 !is_loop(*aedge) &&
						 is_right_node(*aedge, *anode)) {
				++aedge;
			}
			if (aedge != anode_edges.end())
				return;
			++anode;
		}
	}

	template<template<class, class> class Q>
	debruijn_graph::edge debruijn_graph::edge_iterator<Q>::operator*(void) const {
		return *aedge;
	}

	template<template<class, class> class Q>
		void debruijn_graph::edge_iterator<Q>::operator++(void) {
			++aedge;
			while (aedge != anode_edges.end() &&
						 !is_loop(*aedge) &&
						 is_right_node(*aedge, *anode))
				++aedge;
			while(aedge == anode_edges.end()) {
				++anode;
				if (anode != nodeend) {
					anode_edges = dbg.edges(*anode);
					aedge = anode_edges.begin();
					while (aedge != anode_edges.end() &&
								 !is_loop(*aedge) &&
								 is_right_node(*aedge, *anode))
						++aedge;
				} else {
					anode_edges = set<node>();
					aedge = anode_edges.end();
					break;
				}
			}
		}

	template<template<class, class> class Q>
	bool
	debruijn_graph::edge_iterator<Q>::operator!=(const debruijn_graph::edge_iterator<Q> &other) const {
			return aedge != anode_edges.end() || other.aedge != other.anode_edges.end();
		}

	template<template<class, class> class Q>
	debruijn_graph::edges_container<Q>::edges_container(const debruijn_graph &g,
																											const std::unordered_set<node> strts,
																											uint64_t dpth)
		:dbg(g), starts(strts), depth(dpth) {
		}

	template<template<class, class> class Q>
	debruijn_graph::edge_iterator<Q>
	debruijn_graph::edges_container<Q>::begin(void) const {
			return edge_iterator<Q>(dbg, starts, depth);
		}

	template<template<class, class> class Q>
	debruijn_graph::edge_iterator<Q>
	debruijn_graph::edges_container<Q>::end(void) const {
			return edge_iterator<Q>(dbg, unordered_set<node>(), depth);
		}

	debruijn_graph::edges_container<stack>
	debruijn_graph::reachable_edges(void) const {
		return edges_container<stack>(*this, get_all_ends(), infinite_depth);
	}

	debruijn_graph::edges_container<stack>
	debruijn_graph::connected_component_edges(node n) const {
		return edges_container<stack>(*this,
																	unordered_set<node>({n}),
																	infinite_depth);
	}

	debruijn_graph::edges_container<queue>
	debruijn_graph::nearby_edges(node n, uint64_t depth) const {
		return edges_container<queue>(*this, unordered_set<node>({n}), depth);
	}

	set<debruijn_graph::edge> debruijn_graph::left_edges(node n) const {
		set<edge> result;
		for (const auto b : bases) {
			edge e = b + n;
			if (abundance(e) && !is_duplex_edge(n, e))
				result.insert(e);
		}
		return result;
	}
	set<debruijn_graph::edge> debruijn_graph::right_edges(node n) const {
		set<edge> result;
		for (const auto b : bases) {
			edge e = n + b;
			if (abundance(e) && !is_duplex_edge(n, e))
				result.insert(e);
		}
		return result;
	}
	set<debruijn_graph::edge> debruijn_graph::duplex_edges(node n) const {
		set<edge> result;
		for (const auto b : bases) {
			edge e1 = b + n;
			if (abundance(e1) && is_duplex_edge(n, e1))
				result.insert(e1);
			edge e2 = n + b;
			if (abundance(e2) && is_duplex_edge(n, e2))
				result.insert(e2);
		}
		return result;
	}
	set<debruijn_graph::edge> debruijn_graph::edges(node n) const {
		set<edge> result;
		for (const auto b : bases) {
			if (abundance(n + b))
				result.insert(n + b);
			if (abundance(b + n))
				result.insert(b + n);
		}
		return result;
	}

	set<debruijn_graph::node> debruijn_graph::left_neighbors(node n) const {
		set<node> result;
		for (const auto b : bases)
			if (abundance(b + n) && !is_duplex_edge(n, b + n))
				result.insert(b >> n);
		return result;
	}
	set<debruijn_graph::node> debruijn_graph::right_neighbors(node n) const {
		set<node> result;
		for (const auto b : bases)
			if (abundance(n + b) && !is_duplex_edge(n, n + b))
				result.insert(n << b);
		return result;
	}
	set<debruijn_graph::node> debruijn_graph::duplex_neighbors(node n) const {
		set<node> result;
		for (const auto b : bases) {
			if (abundance(b + n) && is_duplex_edge(n, b + n))
				result.insert(b >> n);
			if (abundance(n + b) && is_duplex_edge(n, n + b))
				result.insert(n << b);
		}
		return result;
	}
	set<debruijn_graph::node> debruijn_graph::neighbors(node n) const {
		set<node> result;
		for (const auto b : bases) {
			if (abundance(b + n))
				result.insert(b >> n);
			if (abundance(n + b))
				result.insert(n << b);
		}
		return result;
	}

	uint64_t debruijn_graph::left_degree(node n) const {
		return left_edges(n).size();
		return left_edges(n).size();
	}
	uint64_t debruijn_graph::right_degree(node n) const {
		return right_edges(n).size();
	}
	uint64_t debruijn_graph::duplex_degree(node n) const {
		return duplex_edges(n).size();
	}
	uint64_t debruijn_graph::degree(node n) const {
		return edges(n).size();
	}

	uint64_t debruijn_graph::left_abundance(node n) const {
		uint64_t sum = 0;
		for (const auto &e : left_edges(n))
			sum += abundance(e);
		return sum;
	}
	uint64_t debruijn_graph::right_abundance(node n) const {
		uint64_t sum = 0;
		for (const auto &e : right_edges(n))
			sum += abundance(e);
		return sum;
	}
	uint64_t debruijn_graph::duplex_abundance(node n) const {
		uint64_t sum = 0;
		for (const auto &e : duplex_edges(n))
			sum += abundance(e);
		return sum;
	}

	bool debruijn_graph::is_leaf(node n) const {
		return degree(n) < 2;
	}
	bool debruijn_graph::is_leaf_edge(edge e) const {
		return is_leaf(left_node(e)) || is_leaf(right_node(e));
	}

	template<template<class, class> class Q>
		void debruijn_graph::print_graph(string name,
																		 nodes_container<Q> nodes,
																		 edges_container<Q> edges,
																		 ostream &os) const {
			os << "digraph " << name << " {" << endl;
			for (const auto & n : nodes) {
				os << "  "
					<< (string)n
					<< "[ label=\""
					<< (string)n
					<< "\\n("
					<< get_left_ends(n) << ", "
					<< get_right_ends(n) 
					<< ")\", style=filled, fillcolor="
					<< ((get_left_ends(n) + get_right_ends(n)) ? "grey" : "white")
					<< ", shape="
					<< (can_have_duplex_edges(n) ? "diamond" : "ellipse")
					<<" ];" << endl;
			}
			for (const auto & e : edges) {
				node ln = left_node(e);
				node rn = right_node(e);
				os << "  "
					<< (string)ln << " -> "
					<< (string)rn << "[ label=\"("
					<< (string)e << ", "
					<< get_occurences(e)
					<< ")\", tailport="
					<< (is_left_edge(ln, e) ? "w" : is_right_edge(ln, e) ? "e" : "s")
					<< ", headport="
					<< (is_left_edge(rn, e) ? "w" : is_right_edge(rn, e) ? "e" : "s")
					<<" ];" << endl;
			}
			os << "}" << endl;
		}

	void debruijn_graph::print_reachable_graph(std::string name,
																						 std::ostream &os) const {
		print_graph(name, reachable_nodes(), reachable_edges(), os);
	}

	void debruijn_graph::print_nearby_graph(std::string name,
																					std::ostream &os,
																					node n,
																					uint64_t dist) const {
		print_graph(name, nearby_nodes(n, dist), nearby_edges(n, dist), os);
	}

	void debruijn_graph::record_occurence(edge k) {
		occurences.increment(k, 1);
	}

	uint64_t debruijn_graph::get_occurences(edge e) const {
		return occurences.count(e);
	}

	void debruijn_graph::record_left_end(node n) {
		begins.increment(n, 1);
	}

	uint64_t debruijn_graph::get_left_ends(node n) const {
		return begins.count(n);
	}

	//set<debruijn_graph::node> debruijn_graph::get_all_left_ends(void) const {
	//set<node> result;
	//for (const auto & it = begins.begin(); it != begins.end(); ++it)
	//result.insert(*it);
	////for (const auto & p : left_ends)
	////result.insert(p.first);
	//return result;
	//}

	//set<debruijn_graph::node> debruijn_graph::get_all_right_ends(void) const {
	//set<node> result;
	//for (const auto & it = ends.begin(); it != ends.end(); ++it)
	//result.insert(*it);
	////for (const auto & p : right_ends)
	////result.insert(p.first);
	//return result;
	//}

	unordered_set<debruijn_graph::node> debruijn_graph::get_all_ends(void) const {
		unordered_set<node> result;
		for (auto it = begins.begin(); it != begins.end(); ++it)
			result.insert(*it);
		for (auto it = ends.begin(); it != ends.end(); ++it)
			result.insert(*it);
		//for (const auto & p : left_ends)
		//result.insert(p.first);
		//for (const auto & p : right_ends)
		//result.insert(p.first);
		return result;
	}

	void debruijn_graph::record_right_end(node n) {
		ends.increment(n, 1);
		//if (right_ends.count(n))
		//right_ends[n]++;
		//else
		//right_ends[n] = 1;
	}

	uint64_t debruijn_graph::get_right_ends(node n) const {
		return ends.count(n);
		//return right_ends.count(n) ? right_ends.at(n) : 0;
	}

	void debruijn_graph::record_end(node n, edge e) {
		if (is_left_edge(n, e)) {
			record_left_end(n);
			//cout << "Begin node: " << string(n) << endl;
			//cout << "Begin edge: " <<  string(e) << endl;
		} else if (is_right_edge(n, e)) {
			record_right_end(n);
			//cout << "End node: " << string(n) << endl;
			//cout << "End edge: " <<  string(e) << endl;
		}
	}

	int64_t debruijn_graph::left_invariant_value(node n) const
	{
		uint64_t sum = 0;
		for (const auto &e : left_edges(n))
			sum += (is_loop(e) ? 2 : 1) * abundance(e);
		return sum - get_left_ends(n);
	}
	int64_t debruijn_graph::right_invariant_value(node n) const
	{
		uint64_t sum = 0;
		for (const auto &e : right_edges(n))
			sum += (is_loop(e) ? 2 : 1) * abundance(e);
		return sum - get_right_ends(n);
	}
	int64_t debruijn_graph::invariant_discrepancy(node n) const {
		return right_invariant_value(n) - left_invariant_value(n);
	}

	int64_t debruijn_graph::invariant_implied_abundance_correction(node n,
																																 edge e) const {
		if (is_left_edge(n, e))
			return invariant_discrepancy(n);
		else if (is_right_edge(n, e))
			return -invariant_discrepancy(n);
		else // duplex edge
			return 0;
	}

	bool debruijn_graph::abundance_is_correct_whp(edge e) const {
		// FIXME: need to check for small cycles, and adjust 3
		for (const auto & n : nearby_nodes(left_node(e), 3))
			if (invariant_discrepancy(n) != 0)
				return false;
		return period(e) > 6;
	}

	debruijn_graph::unresolved_edge_set::unresolved_edge_set(void)
		: unordered_set<edge>() {}
	
	bool debruijn_graph::abundance_is_definitely_correct(edge e,
																											 unresolved_edges_map &mbi) const {
		if (abundance(e) == 0)
			return true;
		
		uint64_t hash = occurences.get_hash(e);
		return mbi.count(hash) == 0 || mbi[hash].count(e) == 0;
	}

	bool
	debruijn_graph::left_side_is_definitely_correct(node n,
																									unresolved_edges_map &mbi) const {
		for (const auto & e : left_edges(n)) {
			if (!abundance_is_definitely_correct(e, mbi))
				return false;
		}
		return true;
	}

	bool
	debruijn_graph::right_side_is_definitely_correct(node n,
																									 unresolved_edges_map &mbi) const {
		for (const auto & e : right_edges(n)) {
			if (!abundance_is_definitely_correct(e, mbi))
				return false;
		}
		return true;
	}

	debruijn_graph::edge
	debruijn_graph::find_single_error_edge(node n,
																				 unresolved_edges_map &mbi) const {
		vector<edge> error_edges;
		for (const auto & e : edges(n))
			if (!abundance_is_definitely_correct(e, mbi))
				error_edges.push_back(e);
		if (error_edges.size() == 1)
			return error_edges.front();
		return edge();
	}

	void debruijn_graph::add_correction(edge e, int64_t correction) {
		assert(correction < 0);
		corrections[e] -= correction;
	}

	uint64_t debruijn_graph::get_corrections(edge e) const {
		return corrections.count(e) ? corrections.at(e) : 0;
	}

	void debruijn_graph::mark_edge_as_correct(edge e,
																						unresolved_edges_map &mbi,
																						unordered_set<edge> &work_queue,
																						bool can_use_hashing_rules) {
		uint64_t hash = occurences.get_hash(e);

		// If the edge is currently in mbi
		if (mbi.count(hash) > 0 && mbi[hash].count(e) > 0) {

			// Remove it, update the amount of unresolved count, and add its
			// nodes to the work-queue
			mbi[hash].erase(e);
			for (const auto & n : nodes(e))
				work_queue.insert(n);

			if (can_use_hashing_rules) {
				mbi[hash].unresolved_count -= abundance(e);

				// If any edge currently wants more count than is available,
				// then we know it wants too much
				for (const auto & o : mbi[hash]) {
					if (abundance(o) > mbi[hash].unresolved_count) {
						add_correction(o, mbi[hash].unresolved_count - abundance(o));
						for (const auto & n : nodes(o))
							work_queue.insert(n);
					}
				}
			
				// Figure out how much count the remaining edges currently want
				uint64_t current_total_edge_count = 0;
				for (const auto & o : mbi[hash])
					current_total_edge_count += abundance(o);

				// If the amount of total wanted count is the same as the amount
				// of available count, then we know all the remaining edges are
				// correct.
				if (mbi[hash].unresolved_count == current_total_edge_count) {
					for (const auto & o : mbi[hash])
						for (const auto & n : nodes(o))
							work_queue.insert(n);
					mbi.erase(hash);
				}
			}
		}
	}

	void debruijn_graph::do_work(node n,
															 unresolved_edges_map &mbi,
															 unordered_set<node> &work_queue,
															 bool can_use_hashing_rules) {
		// If there's no invariant discrepancy, and we know that all the
		// edges on one side are correct, then we can conclude that all
		// the edges on the other side are also correct.
		if (!invariant_discrepancy(n)) {
			if (left_side_is_definitely_correct(n, mbi))
				for (const auto & re : right_edges(n))
					mark_edge_as_correct(re, mbi, work_queue, can_use_hashing_rules);
			else if (right_side_is_definitely_correct(n, mbi))
				for (const auto & le : left_edges(n))
					mark_edge_as_correct(le, mbi, work_queue, can_use_hashing_rules);

			// If we know that there are no "through paths" on the left
			// side, and we know there are no paths ending on the right
			// side, then there must not be _any_ paths on the right side.
		} else if (left_invariant_value(n) == 0 && get_right_ends(n) == 0) {
			for (const auto & re : right_edges(n)) {
				add_correction(re, -abundance(re));
				mark_edge_as_correct(re, mbi, work_queue, can_use_hashing_rules);
			}
			// And vice versa
		} else if (right_invariant_value(n) == 0 && get_left_ends(n) == 0) {
			for (const auto & le : left_edges(n)) {
				add_correction(le, -abundance(le));
				mark_edge_as_correct(le, mbi, work_queue, can_use_hashing_rules);
			}

			// If we know all but one edge is correct, then solve for the
			// abundance of the remaining edge, fix it, and mark it correct.
		} else {
			edge wrong_edge = find_single_error_edge(n, mbi);
			if (wrong_edge != edge()) {
				int64_t correction = invariant_implied_abundance_correction(n,
																																		wrong_edge);
				if (correction != 0) {
					if (is_loop(wrong_edge)) {
						assert(((-correction) % 2) == 0);
						correction /= 2;
					}
					add_correction(wrong_edge, correction);
					mark_edge_as_correct(wrong_edge, mbi, work_queue, can_use_hashing_rules);
				}
			}
		}
	}

	template<template<class, class> class Q>
	void debruijn_graph::update_corrections(nodes_container<Q> &nc,
																					bool is_entire_graph,
																					bool profile) {
		unordered_set<node> work_queue;
		unresolved_edges_map mbi;
		struct timeval start1, end1;

		gettimeofday(&start1, NULL);

		// Find all the edges that are near an invariant discrepancy or
		// might be in a small cycle.  Mark them as "might be incorrect"
		for (const auto & n : nc)
			if (invariant_discrepancy(n) || period(n) <= 6)
				for (const auto & e : nearby_edges(n, 3)) {
					uint64_t hash = occurences.get_hash(e);
					mbi[hash].insert(e);
				}

		// The hash-collision-based error correction rules are only valid
		// when we are working on the entire graph
		if (is_entire_graph) {
			// Ensure that all edges that hash-collide with an edge that
			// "might be incorrect" are also marked as "might be incorrect".
			for (const auto & n : nc)
				for (const auto & e : edges(n)) {
					uint64_t h = occurences.get_hash(e);
					if (mbi.count(h) > 0)
						mbi[h].insert(e);
				}
		
			// Now remove all the hash sets that have a single item -- the
			// edges in such singleton sets must be correct.
			unordered_set<uint64_t> to_be_deleted;
			for (const auto & p : mbi)
				if (p.second.size() <= 1)
					to_be_deleted.insert(p.first);
			for (const auto h : to_be_deleted)
				mbi.erase(h);

			// Finally, initialize the unresolved count info for each
			// unresolved edge set.
			for (auto & p : mbi)
				p.second.unresolved_count = get_occurences(*p.second.begin());
		}
		
		
		uint64_t mbi_hashes = 0;
		uint64_t mbi_edges = 0;
		// Prime the work queue with all nodes that touch an
		// edge that "might be incorrect"
		for (const auto & p : mbi) {
			for (const auto & e : p.second) {
				work_queue.insert(left_node(e));
				work_queue.insert(right_node(e));
				mbi_edges++;
			}
			mbi_hashes++;
		}

		if (profile) {
			cout << "Initial work queue size: " << work_queue.size() << endl;
			cout << "Num mbi hashes: " << mbi_hashes << endl;
			cout << "Num mbi edges: " << mbi_edges << endl;
		}

		// Process the work queue, looking for opportunities to correct
		// errors and declare edges as not "might be incorrect"
		uint64_t niterations = 0;
		while (!work_queue.empty()) {
			node n = *(work_queue.begin());
			work_queue.erase(n);
			do_work(n, mbi, work_queue, is_entire_graph);
			niterations++;
		}

		gettimeofday(&end1, NULL);
		if (profile)
			print_time_elapsed("", &start1, &end1);

		if (profile) {
			cout << "Number of iterations: " << niterations << endl;
			cout << "Number of corrections: " << corrections.size() << endl;
			cout << "Number of remaining suspicious edges: " << mbi.size() << endl;
			uint64_t num_invariant_discrepancies = 0;
			for (const auto & n : nc)
				if (invariant_discrepancy(n)) {
					num_invariant_discrepancies++;
					ofstream idfile("/tmp/id-" + (string)n + ".dot");
					print_nearby_graph((string)n, idfile, n, 6);
				}
			cout << "Number of remaining invariant discrepancies: "
					 << num_invariant_discrepancies << endl;
		}
	}

	void debruijn_graph::update_corrections(node n) {
		nodes_container<queue> nc = nearby_nodes(n, 3);
		update_corrections(nc, false, false);
	}
	void debruijn_graph::update_corrections(void) {
		nodes_container<stack> nc = reachable_nodes();
		update_corrections(nc, true, true);
	}

	debruijn_graph::node left_node(debruijn_graph::edge e) {
		return prefix(e, e.len-1);
	}
	debruijn_graph::node right_node(debruijn_graph::edge e) {
		return suffix(e, e.len-1);
	}
	set<debruijn_graph::node> nodes(debruijn_graph::edge e) {
		return set<debruijn_graph::node>({ left_node(e), right_node(e) });
	}

	bool is_left_edge(debruijn_graph::node n, debruijn_graph::edge e) {
		for (const auto b : bases)
			if (e == canonicalize(b + n))
				return !is_duplex_edge(n, e);
		return false;
	}
	bool is_right_edge(debruijn_graph::node n, debruijn_graph::edge e) {
		for (const auto b : bases)
			if (e == canonicalize(n + b))
				return !is_duplex_edge(n, e);
		return false;
	}
	bool is_duplex_edge(debruijn_graph::node n, debruijn_graph::edge e) {
		bool lefty = false, righty = false;
		for (const auto b : bases) {
			righty |= (e == canonicalize(n + b));
			lefty |= (e == canonicalize(b + n));
		}
		return righty && lefty;
	}  

	bool is_left_node(debruijn_graph::edge e, debruijn_graph::node n) {
		return n == canonicalize(prefix(e, e.len-1));
	}
	bool is_right_node(debruijn_graph::edge e, debruijn_graph::node n) {
		return n == canonicalize(suffix(e, e.len-1));
	}

	bool is_self_revcomp(debruijn_graph::node n) {
		return n == -n;
	}
	bool is_constant(debruijn_graph::node n) {
		uint64_t b = n.val & 3ULL;
		uint64_t val = b * 0x5555555555555555;
		debruijn_graph::node n2(n.len, val & BITMASK(2*n.len));
		return n == n2;
	}
	bool can_have_duplex_edges(debruijn_graph::node n) {
		return is_self_revcomp(n) || is_constant(n);
	}

	bool is_loop(debruijn_graph::edge e) {
		return nodes(e).size() == 1;
	}

};

//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#define K 28 // yuck

using namespace dna;

void kmer_tests(void)
{
	kmer aagt1("AAGT");
	kmer aagt2(A);
	aagt2 = aagt2 + A + G + T;
	assert(aagt1 == aagt2);
	assert(aagt1 + A != aagt1 + C);
	assert((aagt1 << A) == kmer("AGTA"));
	assert(-aagt1 == kmer("ACTT"));
	assert(aagt1 + (-aagt2) == kmer("AAGTACTT"));
	assert(suffix(aagt1, 3) == kmer("AGT"));
	assert(prefix(aagt1, 3) == kmer("AAG"));

	canonical_kmer can_aagt = aagt1;
	assert(can_aagt == kmer("ACTT") || can_aagt == kmer("AAGT"));
	canonical_kmer can_tttt("TTTT"), can_aaaa("AAAA");
	assert(can_tttt == can_aaaa);

	assert(period(kmer("AAAAAAAAAAAAA")) == 1);
	assert(period(kmer("ACACACACACACAC")) == 2);
	assert(period(kmer("ACTACTACTACTACTACT")) == 3);

	string s = "ACTGGTGCGCCCGCCTTTTTGGGGGGTCCGTCC";
	for (uint64_t periodicity = 1; periodicity <=32; periodicity++) {
		string base = s.substr(0, periodicity);
		string k;
		for (uint64_t i = 0; i < 32; i++)
			k += base;
		for (uint64_t len = periodicity; len <=32; len++) {
			kmer kk(k.substr(0,len));
			cout << (string)kk << " "
					 << periodicity << " " << len << " " << period(kk) << endl;
			assert(period(kk) == periodicity);
		}
	}
}

#if 0
void dbg_tests(void)
{
	debruijn_graph dbg(5, 1ULL << 60);

	assert(dbg.abundance(kmer("AAAGA")) == 0);
	dbg.insert_read("AAAGA");
	assert(dbg.abundance(kmer("AAAGA")) == 1);

	assert(dbg.is_loop(kmer("AAAAA")));
	assert(!dbg.is_loop(kmer("ACCGA")));

	dbg.insert_read("CGAATTGGA");

	assert(dbg.is_leaf(kmer("CGAA")));
	assert(dbg.is_leaf(kmer("TGGA")));
	assert(dbg.is_leaf(kmer("AAAG")));
	assert(dbg.is_leaf(kmer("AAGA")));
	assert(!dbg.is_leaf(kmer("AATT")));

	assert(dbg.is_leaf_edge(kmer("CGAAT")));
	assert(dbg.is_leaf_edge(kmer("TTGGA")));
	assert(dbg.is_leaf_edge(kmer("AAAGA")));
	assert(!dbg.is_leaf_edge(kmer("AATTG")));

	assert(dbg.can_have_duplex_edges(kmer("AAAA")));
	assert(dbg.can_have_duplex_edges(kmer("ACGT")));
	assert(dbg.can_have_duplex_edges(kmer("AATT")));
	assert(!dbg.can_have_duplex_edges(kmer("ACCA")));

	assert(dbg.is_constant(kmer("AAAA")));
	assert(!dbg.is_constant(kmer("AAAC")));

	assert(dbg.is_self_revcomp(kmer("ACGT")));
	assert(!dbg.is_self_revcomp(kmer("ACGA")));

	{
		debruijn_graph::edge e("GGGTG");
		debruijn_graph::node ln("GGGT");
		debruijn_graph::node rn("GGTG");
		assert(dbg.is_left_node(e, ln));
		assert(!dbg.is_left_node(e, rn));
		assert(dbg.is_right_node(e, rn));
		assert(!dbg.is_right_node(e, ln));
	}

	{
		debruijn_graph::node n("CATG");
		debruijn_graph::edge de1("CATGC");
		debruijn_graph::edge de2("ACATG");
		assert(dbg.is_duplex_edge(n, de1));
		assert(dbg.is_duplex_edge(n, de2));
		assert(!dbg.is_left_edge(n, de1));
		assert(!dbg.is_right_edge(n, de1));
	}

	{
		debruijn_graph::node n("CCCG");
		debruijn_graph::edge le(A + n);
		debruijn_graph::edge re(n + A);
		assert(!dbg.is_duplex_edge(n, le));
		assert(!dbg.is_duplex_edge(n, re));
		assert(dbg.is_left_edge(n, le));
		assert(!dbg.is_left_edge(n, re));
		assert(dbg.is_right_edge(n, re));
		assert(!dbg.is_right_edge(n, le));
	}

	{
		debruijn_graph::edge e("AAAGT"); 
		debruijn_graph::node n1("AAAG");
		debruijn_graph::node n2("AAGT");
		debruijn_graph::node n3("CCCC");
		assert(dbg.right_node(e) == n2);
		assert(dbg.left_node(e) == n1);
		assert(dbg.left_node(e) != n3);
		assert(dbg.right_node(e) != n3);
		assert(dbg.nodes(e) == set<debruijn_graph::node>({n1, n2}));
	}

	{
		debruijn_graph::node aleaf("TGGA");
		debruijn_graph::edge itsedge("TTGGA");
		debruijn_graph::node itsneighbor("TTGG");
		assert(dbg.left_edges(aleaf) == set<debruijn_graph::edge>({itsedge}));
		assert(dbg.right_edges(aleaf) == set<debruijn_graph::edge>());
		assert(dbg.edges(aleaf) == set<debruijn_graph::edge>({itsedge}));
		assert(dbg.left_degree(aleaf) == 1);
		assert(dbg.right_degree(aleaf) == 0);
		assert(dbg.duplex_degree(aleaf) == 0);
		assert(dbg.degree(aleaf) == 1);
		assert(dbg.left_abundance(aleaf) == 1);
		assert(dbg.right_abundance(aleaf) == 0);
		assert(dbg.duplex_abundance(aleaf) == 0);
		assert(dbg.left_neighbors(aleaf) == set<debruijn_graph::node>({itsneighbor}));
		assert(dbg.right_neighbors(aleaf) == set<debruijn_graph::node>());
		assert(dbg.duplex_neighbors(aleaf) == set<debruijn_graph::node>());
		assert(dbg.neighbors(aleaf) == set<debruijn_graph::node>({itsneighbor}));

	}

	{
		debruijn_graph::node aleaf2("TTCG");
		debruijn_graph::edge itsedge2("CGAAT");
		debruijn_graph::node itsneighbor2("GAAT");
		assert(dbg.left_edges(aleaf2) == set<debruijn_graph::edge>({itsedge2}));
		assert(dbg.right_edges(aleaf2) == set<debruijn_graph::edge>());
		assert(dbg.edges(aleaf2) == set<debruijn_graph::edge>({itsedge2}));
		assert(dbg.left_degree(aleaf2) == 1);
		assert(dbg.right_degree(aleaf2) == 0);
		assert(dbg.duplex_degree(aleaf2) == 0);
		assert(dbg.degree(aleaf2) == 1);
		assert(dbg.left_abundance(aleaf2) == 1);
		assert(dbg.right_abundance(aleaf2) == 0);
		assert(dbg.duplex_abundance(aleaf2) == 0);
		assert(dbg.left_neighbors(aleaf2) ==
					 set<debruijn_graph::node>({itsneighbor2}));
		assert(dbg.right_neighbors(aleaf2) == set<debruijn_graph::node>());
		assert(dbg.duplex_neighbors(aleaf2) == set<debruijn_graph::node>());
		assert(dbg.neighbors(aleaf2) == set<debruijn_graph::node>({itsneighbor2}));
	}

	{
		debruijn_graph::node anode("GAAT");
		debruijn_graph::edge ledge("CGAAT");
		debruijn_graph::edge redge("GAATT");
		debruijn_graph::node lneighbor("CGAA");
		debruijn_graph::node rneighbor("AATT");
		assert(dbg.left_edges(anode) == set<debruijn_graph::edge>({ledge}));
		assert(dbg.right_edges(anode) == set<debruijn_graph::edge>({redge}));
		assert(dbg.edges(anode) == set<debruijn_graph::edge>({ledge, redge}));
		assert(dbg.left_degree(anode) == 1);
		assert(dbg.right_degree(anode) == 1);
		assert(dbg.duplex_degree(anode) == 0);
		assert(dbg.degree(anode) == 2);
		assert(dbg.left_abundance(anode) == 1);
		assert(dbg.right_abundance(anode) == 1);
		assert(dbg.duplex_abundance(anode) == 0);
		assert(dbg.left_neighbors(anode) == set<debruijn_graph::node>({lneighbor}));
		assert(dbg.right_neighbors(anode) == set<debruijn_graph::node>({rneighbor}));
		assert(dbg.duplex_neighbors(anode) == set<debruijn_graph::node>());
		assert(dbg.neighbors(anode) ==
					 set<debruijn_graph::node>({lneighbor, rneighbor}));
	}

	{
		debruijn_graph::node dnode("AATT");
		debruijn_graph::edge dedge1("GAATT");
		debruijn_graph::edge dedge2("AATTG");
		debruijn_graph::node dnbr1("ATTG");
		debruijn_graph::node dnbr2("GAAT");
		assert(dbg.is_duplex_edge(dnode, dedge1));
		assert(dbg.is_duplex_edge(dnode, dedge2));
		assert(dbg.left_edges(dnode) == set<debruijn_graph::edge>());
		assert(dbg.right_edges(dnode) == set<debruijn_graph::edge>());
		assert(dbg.duplex_edges(dnode) ==
					 set<debruijn_graph::edge>({dedge1, dedge2}));
		assert(dbg.edges(dnode) == set<debruijn_graph::edge>({dedge1, dedge2}));
		assert(dbg.left_degree(dnode) == 0);
		assert(dbg.right_degree(dnode) == 0);
		assert(dbg.duplex_degree(dnode) == 2);
		assert(dbg.degree(dnode) == 2);
		assert(dbg.left_abundance(dnode) == 0);
		assert(dbg.right_abundance(dnode) == 0);
		assert(dbg.duplex_abundance(dnode) == 2);
		assert(dbg.left_neighbors(dnode) == set<debruijn_graph::node>());
		assert(dbg.right_neighbors(dnode) == set<debruijn_graph::node>());
		assert(dbg.duplex_neighbors(dnode) ==
					 set<debruijn_graph::node>({dnbr1, dnbr2}));
		assert(dbg.neighbors(dnode) == set<debruijn_graph::node>({dnbr1, dnbr2}));
	}
}

int test_main(int argc, char **argv)
{
	kmer_tests();
	dbg_tests();

	dna::debruijn_graph dbg(K, 1ULL << 5);

	dbg.insert_read("ACTGAATGC");
	// node revc    cano
	// -----------------
	// ACTG CAGT    ACTG
	// CTGA TCAG    TCAG
	// TGAA TTCA    TGAA
	// GAAT ATTC    GAAT
	// AATG CATT    AATG
	// ATGC GCAT    GCAT

	// edge   revc   canon
	// -------------------
	// ACTGA  TCAGT  TCAGT
	// CTGAA  TTCAG  TTCAG
	// TGAAT  ATTCA  TGAAT  TGAA -> GAAT
	// GAATG  CATTC  GAATG
	// AATGC  GCATT  GCATT

	ofstream before_file("/tmp/before.dot");
	dbg.print_graph("before",
									dbg.reachable_nodes(),
									dbg.reachable_edges(),
									before_file);

	for (const auto & e : dbg.nearby_edges(debruijn_graph::node("TGCA"), 3))
		cout << (string)e << endl;

	ofstream sub_file("/tmp/sub.dot");
	dbg.print_graph("sub",
									dbg.nearby_nodes(debruijn_graph::node("TGCA"), 3),
									dbg.nearby_edges(debruijn_graph::node("TGCA"), 3),
									sub_file);

	dbg.update_corrections();

	ofstream after_file("/tmp/after.dot");
	dbg.print_graph("after",
									dbg.reachable_nodes(),
									dbg.reachable_edges(),
									after_file);

	return 0;
}

int main(int argc, char **argv)
{
	kmer_tests();
}
#endif

/* Print elapsed time using the start and end timeval */
void print_time_elapsed(string desc, struct timeval* start, struct timeval* end)
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	cout << desc << "Total Time Elapsed: " << to_string(time_elapsed) << 
		"seconds" << endl;
}

int main ( int argc, char *argv[] )
{
	if (argc < 2) {
		cout << "Insufficient arguments!" << endl;
		abort();
	}
	struct timeval start1, end1;
	struct timezone tzp;

	// De-serialize QFs
	QF occurences, begins, ends, true_kmers;
	{
		string exact_ext(".exact");
		string first_ext(".first");
		string last_ext(".last");
		string ser_ext(".ser");
		string exact_cqf = string(argv[1]) + exact_ext;
		string main_cqf = string(argv[1]) + ser_ext;
		string first_cqf = string(argv[1]) + first_ext;
		string last_cqf = string(argv[1]) + last_ext;
		//Initialize the QF
		cout << "Reading CQFs off disk:" << endl;
		qf_deserialize(&true_kmers, exact_cqf.c_str());
		qf_deserialize(&occurences, main_cqf.c_str());
		qf_deserialize(&begins, first_cqf.c_str());
		qf_deserialize(&ends, last_cqf.c_str());
	}

	// Building first dBG
	dna::debruijn_graph
		dbg(K,
				dna::kmer_counter_cqf<dna::debruijn_graph::edge>(K, occurences),
				dna::kmer_counter_lossless_cqf<dna::debruijn_graph::node>(K-1, begins),
				dna::kmer_counter_lossless_cqf<dna::debruijn_graph::node>(K-1, ends),
				map<dna::debruijn_graph::node, uint64_t>());
	
	// ofstream before("/tmp/before.dot");
	// dbg.print_reachable_graph("reachable_graph", before);

	dna::debruijn_graph dbg_orig = dbg;

	cout << "Finding false positive kmers." << endl;
	unordered_set<debruijn_graph::edge> fp_kmers;
	unordered_set<debruijn_graph::edge> abundance_errors;
	gettimeofday(&start1, &tzp);
	for (const auto & e : dbg.reachable_edges()) {
		uint64_t hash = kmercounting::HashUtil::hash_64(e.val, BITMASK(2*K));
		if (qf_count_key_value(&true_kmers, hash, 0) == 0)
			fp_kmers.insert(e);
		if (qf_count_key_value(&true_kmers, hash, 0) != dbg.abundance(e))
			abundance_errors.insert(e);
	}
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);
	cout << "Number of reachable fp kmers: " << fp_kmers.size() << endl;
	cout << "Number of reachable abundance errors kmers: " << abundance_errors.size() << endl;

	cout << "Counting nodes." << endl;
	gettimeofday(&start1, &tzp);
	uint64_t nnodes = 0, nedges = 0;
	for (const auto & n : dbg.reachable_nodes())
		nnodes++;
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);
	cout << "Number of nodes: " << nnodes << endl;

	cout << "Counting edges." << endl;
	gettimeofday(&start1, &tzp);
	for (const auto & e : dbg.reachable_edges())
		nedges++;
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);
	cout << "Number of edges: " << nedges << endl;

	cout << "Correcting abundances:" << endl;
	dbg.update_corrections();

	uint64_t num_fp_kmers_not_found_by_correction = 0;
	for (const auto & e : fp_kmers) {
		if (dbg.abundance(e)) {
			num_fp_kmers_not_found_by_correction++;
			ofstream nbhd("/tmp/fp-" + ((string)e) + ".dot");
			dbg_orig.print_nearby_graph((string)e, nbhd, left_node(e), 6);
		}
	}
	cout << "Num of false-positive k-mers not found by correcting abundance: " << num_fp_kmers_not_found_by_correction << endl;

	uint64_t num_abundance_errors_not_found_by_correction = 0;
	uint64_t num_over_counts = 0;
	uint64_t num_under_counts = 0;
	for (const auto & e : abundance_errors) {
		uint64_t hash = kmercounting::HashUtil::hash_64(e.val, BITMASK(2*K));
		if (dbg.abundance(e) != qf_count_key_value(&true_kmers, hash, 0)) {
			num_abundance_errors_not_found_by_correction++;
			ofstream nbhd("/tmp/ae-" + ((string)e) + ".dot");
			dbg_orig.print_nearby_graph((string)e, nbhd, left_node(e), 6);
			if (dbg.abundance(e) > qf_count_key_value(&true_kmers, hash, 0))
				num_over_counts++;
			else
				num_under_counts++;
		}
	}
	cout << "Num of abundance errors not found by correcting abundance: " << num_abundance_errors_not_found_by_correction << endl;
	cout << "Num of over counts: " << num_over_counts << endl;
	cout << "Num of under counts: " << num_under_counts << endl;

	QFi cfi;
	uint64_t num_fn_kmers = 0;
	qf_iterator(&true_kmers, &cfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&cfi, &key, &value, &count);
		uint64_t item = kmercounting::HashUtil::hash_64i(key, BITMASK(2*K));

		dna::debruijn_graph::edge e = debruijn_graph::edge(K, item);
		if (dbg.abundance(e) == 0) {
			num_fn_kmers++;
			ofstream nbhd("/tmp/fn-" + ((string)e) + ".dot");
			dbg_orig.print_nearby_graph((string)e, nbhd, left_node(e), 6);
		}
	} while (!qfi_next(&cfi));
	cout << "Number of false negative kmers: " << num_fn_kmers << endl;

	cout << "Testing local correction performance" << endl;
	gettimeofday(&start1, &tzp);
	uint64_t random_query_abundance_errors = 0;
	for (uint64_t i = 0; i < 1000000; i++) {
		dna::debruijn_graph::edge e(K, random() & BITMASK(2*K));
		dbg.update_corrections(left_node(e));
		uint64_t hash = kmercounting::HashUtil::hash_64(e.val, BITMASK(2*K));
		assert(dbg.abundance(e) >= qf_count_key_value(&true_kmers, hash, 0));
		if (dbg.abundance(e) != qf_count_key_value(&true_kmers, hash, 0)) {
			random_query_abundance_errors++;		
		}
	}
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);
	cout << "Number of random query abundance errors: " << random_query_abundance_errors << endl;

	cout << "Testing navigation query performance." << endl;
	nedges = 0;
	gettimeofday(&start1, &tzp);
	for (const auto & e : dbg.reachable_edges())
		nedges++;
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);
	cout << "Number of edges: " << nedges << endl;

	cout << "Testing random query performance." << endl;
	gettimeofday(&start1, &tzp);
	for (uint64_t i = 0; i < 1000000; i++) {
		dna::debruijn_graph::edge e(K, random() & BITMASK(2*K));
		uint64_t hash = kmercounting::HashUtil::hash_64(e.val, BITMASK(2*K));
		assert(dbg.abundance(e) >= qf_count_key_value(&true_kmers, hash, 0));
		if (dbg.abundance(e) != qf_count_key_value(&true_kmers, hash, 0)) {
			random_query_abundance_errors++;		
		}
	}
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);
	cout << "Number of random query abundance errors: " << random_query_abundance_errors << endl;

	return 0;
}
