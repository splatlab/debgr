/*
 * =====================================================================================
 *
 *       Filename:  kmer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2016 05:06:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * =====================================================================================
 */

#ifndef _KMER_H_
#define _KMER_H_

#include <stdio.h>
#include <string>

#define K 28

enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)

using namespace std;

namespace kmercounting {

	class kmer {
		public:
			static inline char map_int(uint8_t base);
			static inline uint8_t map_base(char base);
			static uint64_t str_to_int(string str);
			static string int_to_str(uint64_t kmer);
			static inline int reverse_complement_base(int x);
			static uint64_t reverse_complement(uint64_t kmer);
			static uint64_t reverse_complement_node(uint64_t kmer);
			static inline bool compare_kmers(uint64_t kmer, uint64_t kmer_rev);
			static inline unsigned __int128 word_reverse_complement(unsigned __int128 w);
			static inline int64_t word_reverse_complement(uint64_t w);
			static inline uint32_t word_reverse_complement(uint32_t w);
			static inline uint64_t prefix(uint64_t k, int len) { return k >> (2*(K - len)); }
			static inline uint64_t suffix(uint64_t k, int len) { return k & BITMASK(2*len); }

		private:
			kmer();
	};

	/*return the integer representation of the base */
	inline char kmer::map_int(uint8_t base)
	{
		switch(base) {
			case DNA_MAP::A: { return 'A'; }
			case DNA_MAP::T: { return 'T'; }
			case DNA_MAP::C: { return 'C'; }
			case DNA_MAP::G: { return 'G'; }
			default:  { return DNA_MAP::G+1; }
		}
	}

	/*return the integer representation of the base */
	inline uint8_t kmer::map_base(char base)
	{
		switch(base) {
			case 'A': { return DNA_MAP::A; }
			case 'T': { return DNA_MAP::T; }
			case 'C': { return DNA_MAP::C; }
			case 'G': { return DNA_MAP::G; }
			default:  { return DNA_MAP::G+1; }
		}
	}

	/**
	 * Converts a string of "ATCG" to a uint64_t
	 * where each character is represented by using only two bits
	 */
	uint64_t str_to_int(string str)
	{
		uint64_t strint = 0;
		for (auto it = str.begin(); it != str.end(); it++) {
			uint8_t curr = 0;
			switch (*it) {
				case 'A': { curr = DNA_MAP::A; break; }
				case 'T': { curr = DNA_MAP::T; break; }
				case 'C': { curr = DNA_MAP::C; break; }
				case 'G': { curr = DNA_MAP::G; break; }
			}
			strint = strint | curr;
			strint = strint << 2;
		}
		return strint >> 2;
	}

	/**
	 * Converts a uint64_t to a string of "ACTG"
	 * where each character is represented by using only two bits
	 */
	string int_to_str(uint64_t kmer)
	{
		uint8_t base;
		string str;
		for (int i=K; i>0; i--) {
			base = (kmer >> (i*2-2)) & 3ULL;
			char chr = kmer::map_int(base);
			str.push_back(chr);
		}
		return str;
	}

	string int_to_str_node(uint64_t kmer)
	{
		uint8_t base;
		string str;
		for (int i=K-1; i>0; i--) {
			base = (kmer >> (i*2-2)) & 3ULL;
			char chr = kmer::map_int(base);
			str.push_back(chr);
		}
		return str;
	}

	/* Return the reverse complement of a base */
	inline int kmer::reverse_complement_base(int x) { return 3 - x; }

	/* Calculate the revsese complement of a kmer */
	uint64_t kmer::reverse_complement(uint64_t kmer)
	{
		uint64_t rc = 0;
		uint8_t base = 0;
		for (int i=0; i<K; i++) {
			base = kmer & 3ULL;
			base = reverse_complement_base(base);
			kmer >>= 2;
			rc |= base;
			rc <<= 2;
		}
		rc >>=2;
		return rc;
	}

	uint64_t kmer::reverse_complement_node(uint64_t kmer)
	{
		uint64_t rc = 0;
		uint8_t base = 0;
		for (int i=0; i<K-1; i++) {
			base = kmer & 3ULL;
			base = reverse_complement_base(base);
			kmer >>= 2;
			rc |= base;
			rc <<= 2;
		}
		rc >>=2;
		return rc;
	}

	/* Compare the kmer and its reverse complement and return the result 
	 * Return true if the kmer is greater than or equal to its
	 * reverse complement. 
	 * */
	inline bool kmer::compare_kmers(uint64_t kmer, uint64_t kmer_rev)
	{
		return kmer >= kmer_rev;
	}

	/* This code is taken from Jellyfish 2.0
	 * git@github.com:gmarcais/Jellyfish.git 
	 * */

	// Checkered mask. cmask<uint16_t, 1> is every other bit on
	// (0x55). cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
	template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
		struct cmask {
			static const U v =
				(cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
		};
	template<typename U, int len>
		struct cmask<U, len, 0> {
			static const U v = 0;
		};

	// Fast reverse complement of one word through bit tweedling.
	inline uint32_t kmer::word_reverse_complement(uint32_t w) {
		typedef uint64_t U;
		w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
		w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
		w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
		w = ( w >> 16                   ) | ( w                    << 16);
		return ((U)-1) - w;
	}

	inline int64_t kmer::word_reverse_complement(uint64_t w) {
		typedef uint64_t U;
		w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
		w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
		w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
		w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
		w = ( w >> 32                   ) | ( w                    << 32);
		return ((U)-1) - w;
	}

#ifdef HAVE_INT128
	inline static unsigned __int128 kmer::word_reverse_complement(unsigned __int128 w) {
		typedef unsigned __int128 U;
		w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
		w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
		w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
		w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
		w = ((w >> 32) & cmask<U, 32>::v) | ((w & cmask<U, 32>::v) << 32);
		w = ( w >> 64                   ) | ( w                    << 64);
		return ((U)-1) - w;
	}
#endif
}

#endif
