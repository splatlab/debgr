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

#include<map>
#include<cinttypes>
#include<string>
#include<set>
#include<utility>
#include<unordered_set>
#include<unordered_map>
#include<stack>
#include<queue>

#include "threadsafe-gqf/gqf.h"

namespace dna {

	/////////////// bases /////////////////
	enum base { C = 0, A = 1, T = 2, G = 3 };
	base operator-(base b); // return the complementary base
	extern const base bases[4];
	extern const std::map<char, base> base_from_char;
	extern const std::map<base, char> base_to_char;

	///////////// kmers /////////////////////
	class kmer {
		public:
			int len;
			uint64_t val;

			kmer(void);
			kmer(base b);
			kmer(int l, uint64_t v);
			kmer(std::string s);

			// Convert to string
			operator std::string() const;
	};

	bool operator<(kmer a, kmer b);
	bool operator==(kmer a, kmer b);
	bool operator!=(kmer a, kmer b);

	// Return the reverse complement of k
	kmer operator-(kmer k);

	kmer canonicalize(kmer k);

	// Return the kmer of length |a| that results from shifting b into a
	// from the right
	kmer operator<<(kmer a, kmer b);

	// Return the kmer of length |b| that results from shifting a into b
	// from the left
	kmer operator>>(kmer a, kmer b);

	// Append two kmers
	kmer operator+(kmer a, kmer b);

	kmer suffix(kmer k, int len);
	kmer prefix(kmer k, int len);

	// The purpose of this class is to enable us to declare containers
	// as holding canonical kmers, e.g. set<canonical_kmer>.  Then all
	// inserts/queries/etc will automatically canonicalize their
	// arguments.
	class canonical_kmer : public kmer {
		public:
			canonical_kmer(void);
			canonical_kmer(base b);
			canonical_kmer(int l, uint64_t v);
			canonical_kmer(std::string s);
			canonical_kmer(kmer k);
	};
}

namespace std {
	template <> struct hash<dna::kmer> {
		size_t operator()(const dna::kmer & x) const;
	};
	template <> struct hash<dna::canonical_kmer> {
		size_t operator()(const dna::canonical_kmer & x) const;
	};
}

namespace dna {
	/////////////// Interfaces for CQFs of kmers ///////////////

	template<class kmer_type> // kmer_type = kmer or canonical_kmer
		class kmer_counter_cqf {
			public:
				QF counts;
				const uint64_t len;

				kmer_counter_cqf(uint64_t l, QF cqf);

				uint64_t get_hash(kmer_type k) const;
				virtual void increment(kmer_type k, uint64_t amount);
				virtual size_t count(kmer_type k) const;
		};

	template<class kmer_type>
		class kmer_counter_lossless_cqf : kmer_counter_cqf<kmer_type> {
			public:
				kmer_counter_lossless_cqf(uint64_t l, QF cqf);

				void increment(kmer_type k, uint64_t amount);
				uint64_t count(kmer_type k) const;
				uint64_t size(void) const;
				
				class iterator {
					public:
						const int64_t len;
						QFi iter;
						iterator(uint64_t l, QFi it);

						kmer_type operator*(void) const;
						void operator++(void);
				};
				iterator begin(void) const;
				iterator end(void) const;
		};

	template<class kmer_type=canonical_kmer>
		bool operator!=(const typename kmer_counter_lossless_cqf<kmer_type>::iterator &a,
										const typename kmer_counter_lossless_cqf<kmer_type>::iterator &b);

	/////////////// CQF-based de Bruijn graphs ///////////////

	class debruijn_graph {
		public:
			typedef canonical_kmer edge; // k-mer
			typedef canonical_kmer node; // (k-1)-mer

			const unsigned int len;
			kmer_counter_cqf<edge> occurences;
			kmer_counter_lossless_cqf<node> begins;
			kmer_counter_lossless_cqf<node> ends;
			std::map<edge, uint64_t> corrections;
			
			//std::map<uint64_t, uint64_t> occurences;
			//std::map<node, uint64_t> left_ends;
			//std::map<node, uint64_t> right_ends;

			debruijn_graph(int len,
										 kmer_counter_cqf<edge> occ,
										 kmer_counter_lossless_cqf<node> begs,
										 kmer_counter_lossless_cqf<node> enz,
										 std::map<edge, uint64_t> cor);

			uint64_t abundance(edge e) const;

			// Simple API, but supports little parallelism
			void insert_read(std::string s);

			// This class is for performing highly parallel inserts
			/* class read_inserter { */
			/* private: */
			/*   debruijn_graph &dbg; */
			/*   read_inserter(debruijn_graph &g); */
			/* public: */
			/*   ~read_inserter(void); */
			/*   void insert_read(string s); */
			/* } */
			/* read_inserter get_read_inserter(void); */

			void update_corrections(void);
			void update_corrections(node);

			static const uint64_t infinite_depth;

			template<template<class, class> class Q> class node_iterator {
				private:
					class work_item {
						public:
							node curr;
							node prev;
							uint64_t depth;
					};

					const debruijn_graph &dbg;
					Q<work_item, std::deque<work_item> > work;
					std::unordered_set<node> visited;

					work_item front(const std::queue<work_item> &) const;
					work_item front(const std::stack<work_item> &) const;
				public:
					node_iterator(const debruijn_graph &dbg,
												const std::unordered_set<node> strts,
												uint64_t depth);
					node operator*(void) const;
					void operator++(void);
					bool operator!=(const node_iterator &other) const;
			};

			template<template<class, class> class Q> class nodes_container {
				private:
					const debruijn_graph &dbg;
					const std::unordered_set<node> starts;
					const uint64_t depth;
				public:
					nodes_container(const debruijn_graph &dbg,
													const std::unordered_set<node> strts,
													uint64_t depth);
					node_iterator<Q> begin(void) const; // may visit some edges twice
					node_iterator<Q> end(void) const;
			};

			nodes_container<std::stack> reachable_nodes(void) const; 
			nodes_container<std::stack> connected_component(node n) const;
			nodes_container<std::queue> nearby_nodes(node n, uint64_t depth) const;

			template<template<class, class> class Q> class edge_iterator {
				private:
					const debruijn_graph &dbg;
					node_iterator<Q> anode;
					node_iterator<Q> nodeend;
					std::set<edge> anode_edges;
					std::set<edge>::iterator aedge;
				public:
					edge_iterator(const debruijn_graph &dbg,
												const std::unordered_set<node> strts,
												uint64_t depth);
					edge operator*(void) const;
					void operator++(void);
					bool operator!=(const edge_iterator &other) const;
			};

			template<template<class, class> class Q> class edges_container {
				private:
					const debruijn_graph &dbg;
					const std::unordered_set<node> starts;
					const uint64_t depth;
				public:
					edges_container(const debruijn_graph &dbg,
													const std::unordered_set<node> strts,
													uint64_t dpth);
					edge_iterator<Q> begin(void) const; // may visit some edges twice
					edge_iterator<Q> end(void) const;
			};

			edges_container<std::stack> reachable_edges(void) const;
			edges_container<std::stack> connected_component_edges(node n) const; 
			edges_container<std::queue> nearby_edges(node n, uint64_t depth) const; 

			std::set<edge> left_edges(node n) const;
			std::set<edge> right_edges(node n) const;
			std::set<edge> duplex_edges(node n) const;
			std::set<edge> edges(node n) const;

			std::set<node> left_neighbors(node n) const;
			std::set<node> right_neighbors(node n) const;
			std::set<node> duplex_neighbors(node n) const;
			std::set<node> neighbors(node n) const;

			uint64_t left_degree(node n) const;
			uint64_t right_degree(node n) const;
			uint64_t duplex_degree(node n) const;
			uint64_t degree(node n) const;

			uint64_t left_abundance(node n) const;
			uint64_t right_abundance(node n) const;
			uint64_t duplex_abundance(node n) const;

			bool is_leaf(node n) const;
			bool is_leaf_edge(edge e) const;

			template<template<class, class> class Q>
				void print_graph(std::string name,
												 nodes_container<Q> nodes,
												 edges_container<Q> edges,
												 std::ostream &os) const;
			void print_reachable_graph(std::string name, std::ostream &os) const;
			void print_nearby_graph(std::string name, std::ostream &os, node n, uint64_t dist) const;

		private:
			void record_occurence(edge e);
			void record_left_end(node n);
			void record_right_end(node n);
			void record_end(node n, edge e);

			uint64_t get_occurences(edge e) const;
			uint64_t get_left_ends(node n) const;
			uint64_t get_right_ends(node n) const;

			std::unordered_set<node> get_all_left_ends(void) const;
			std::unordered_set<node> get_all_right_ends(void) const;
			std::unordered_set<node> get_all_ends(void) const;

			int64_t left_invariant_value(node n) const;
			int64_t right_invariant_value(node n) const;
			int64_t invariant_discrepancy(node n) const;
			int64_t invariant_implied_abundance_correction(node n,
																										 edge e) const;
			/* bool abundance_is_correct_whp_helper(node origin, uint8_t level, */
			/* 					 node n) const; */
			bool abundance_is_correct_whp(edge e) const;

			class unresolved_edge_set : public std::unordered_set<edge> {
			public:
				uint64_t unresolved_count;
				unresolved_edge_set(void);
			};
			
			typedef
				std::unordered_map<uint64_t, unresolved_edge_set>
				unresolved_edges_map;

			bool abundance_is_definitely_correct(edge e,
																					 unresolved_edges_map &mbi) const;
			bool left_side_is_definitely_correct(node n,
																					 unresolved_edges_map &mbi) const;
			bool right_side_is_definitely_correct(node n,
																						unresolved_edges_map &mbi) const;
			edge find_single_error_edge(node n, unresolved_edges_map &mbi) const;
			void add_correction(edge e, int64_t correction);
			uint64_t get_corrections(edge e) const;
			void mark_edge_as_correct(edge e, unresolved_edges_map &mbi,
																std::unordered_set<node> &work_queue,
																bool can_use_hashing_rules);
			void do_work(edge e, unresolved_edges_map &mbi,
									 std::unordered_set<node> &work_queue,
									 bool can_use_hashing_rules);
			void fix_leaf_edge(edge e, unresolved_edges_map &mbi,
												 std::unordered_set<node> &work_queue);

			template<template<class, class> class Q>
				void update_corrections(nodes_container<Q> &, bool, bool);
	};

	debruijn_graph::node left_node(debruijn_graph::edge e);
	debruijn_graph::node right_node(debruijn_graph::edge e);
	std::set<debruijn_graph::node> nodes(debruijn_graph::edge e);

	bool is_left_edge(debruijn_graph::node n, debruijn_graph::edge e);
	bool is_right_edge(debruijn_graph::node n, debruijn_graph::edge e);
	bool is_duplex_edge(debruijn_graph::node n, debruijn_graph::edge e);

	bool is_left_node(debruijn_graph::edge e, debruijn_graph::node n);
	bool is_right_node(debruijn_graph::edge e, debruijn_graph::node n);

	bool is_self_revcomp(debruijn_graph::node n);
	bool is_constant(debruijn_graph::node n);
	bool can_have_duplex_edges(debruijn_graph::node n);
	bool is_loop(debruijn_graph::edge e);

};
