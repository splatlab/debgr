/*
 * ============================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/06/2016 09:56:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/atomic.hpp>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "threadsafe-gqf/gqf.h"
#include "hashutil.h"
#include "chunk.h"
#include "kmer.h"
#include "reader.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define QBITS_LOCAL_QF 16
#define SPARE_EMPTY_LOCAL_QFS 16

using namespace std;
using namespace kmercounting;

typedef struct {
	QF *first_qf;
	QF *last_qf;
	QF *local_qf;
	QF *main_qf;
	QF *exact_qf;
	uint32_t count {0};
	ofstream *kmerlog;
}flush_object;

struct file_pointer {
	std::unique_ptr<reader> freader{nullptr};
	char* part{nullptr};
	char* part_buffer{nullptr};
	int mode{0};
	uint64_t size{0};
	uint64_t part_filled{0};
};

/*create a multi-prod multi-cons queue for storing the chunk of fastq file.*/
boost::lockfree::queue<file_pointer*, boost::lockfree::fixed_sized<true> > ip_files(64);
boost::atomic<int> num_files {0};

/* Count distinct items in a sorted list */
uint64_t count_distinct_kmers(multiset<uint64_t> kmers)
{
	uint64_t cnt = 0;
	uint64_t curr_kmer = 0;

	for(uint64_t kmer: kmers) {
		if (kmer != curr_kmer) {
			curr_kmer = kmer;
			cnt++;
		}
	}
	return cnt;
}

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

/* check if it's the end of the file. */
inline bool is_eof(reader &file_reader, int mode)
{
	if (mode == 0)
		return feof(file_reader.in) != 0;
	else if (mode == 1)
		return gzeof(file_reader.in_gzip) != 0;
	else if (mode == 2)
		return file_reader.bzerror == BZ_STREAM_END;

	return true;
}

/* move the pointer to the end of the next newline. */
bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos)
{
	int64_t i;
	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' ||
																								 part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;
	pos = i+1;

	return true;
}

/* read a part of the fastq file. */
static bool fastq_read_parts(int mode, file_pointer *fp)
{
	char *& _part = (fp->part);
	uint64_t& _size = fp->size;
	char*& part_buffer = (fp->part_buffer);
	uint64_t& part_filled = fp->part_filled;
	reader& file_reader = *(fp->freader.get());

	uint32_t OVERHEAD_SIZE = 65535;
	uint64_t part_size = 1ULL << 23;
	char *part = (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
	memcpy(part, part_buffer, part_filled);

	if(is_eof(file_reader, mode))
		return false;

	uint64_t readed = 0;

	if (mode == 0)
		readed = fread(part+part_filled, 1, part_size, file_reader.in);
	else if (mode == 1)
		readed = gzread(file_reader.in_gzip, part+part_filled, (int) part_size);
	else if (mode == 2)
		readed = BZ2_bzRead(&file_reader.bzerror, file_reader.in_bzip2,
												part+part_filled, (int) part_size);
	else 
		readed = 0;

	int64_t total_filled = part_filled + readed;
	int64_t i;
	if(part_filled >= OVERHEAD_SIZE)
	{
		cout << "Error: Wrong input file!\n";
		exit(EXIT_FAILURE);
	}
	if(is_eof(file_reader, mode))
	{
		_part = part;
		_size = total_filled;
		part = NULL;
		return true;
	}
	// Looking for a FASTQ record at the end of the area
	{
		int64_t line_start[9];
		int32_t j;
		i = total_filled - OVERHEAD_SIZE / 2;
		for(j = 0; j < 9; ++j)
		{
			if(!skip_next_eol(part, i, total_filled))
				break;
			line_start[j] = i;
		}
		_part = part;
		if(j < 9)
			_size = 0;
		else
		{
			int k;
			for(k = 0; k < 4; ++k)
			{
				if(part[line_start[k]+0] == '@' && part[line_start[k+2]+0] == '+')
				{
					if(part[line_start[k+2]+1] == '\n' || part[line_start[k+2]+1] == '\r')
						break;
					if(line_start[k+1]-line_start[k] == line_start[k+3]-line_start[k+2] &&
						 memcmp(part+line_start[k]+1, part+line_start[k+2]+1,
										line_start[k+3]-line_start[k+2]-1) == 0)
						break;
				}
			}
			if(k == 4)
				_size = 0;
			else
				_size = line_start[k];
		}
	}

	copy(_part+_size, _part+total_filled, part_buffer);
	part_filled = total_filled - _size;

	return true;
}

/* dump the contents of a local QF into the main QF */
static void dump_local_qf_to_main(flush_object *obj)
{
	QFi local_cfi;

	if (qf_iterator(obj->local_qf, &local_cfi, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&local_cfi, &key, &value, &count);
			qf_insert(obj->main_qf, key, 0, count, true, true);
		} while (!qfi_next(&local_cfi));
		qf_reset(obj->local_qf);
	}
}

const DNA_MAP bases[4] = {C, A, T, G};

uint64_t canonicalize(uint64_t e)
{
	return e > kmer::reverse_complement(e) ? e : kmer::reverse_complement(e);
}

uint64_t canonicalize_node(uint64_t n)
{
	return n > kmer::reverse_complement_node(n) ? n : kmer::reverse_complement_node(n);
}

bool is_self_revcomp(uint64_t n) {
	return n == kmer::reverse_complement_node(n);
}
	
bool is_constant(uint64_t n) {
	uint64_t b = n & 3ULL;
	uint64_t val = b * 0x5555555555555555;
	uint64_t n2 = val & BITMASK(2*K);
	return n == n2;
}
	
bool can_have_duplex_edges(uint64_t n) {
	return is_self_revcomp(n) || is_constant(n);
}

uint64_t node_r_edge(uint64_t node, DNA_MAP base)
{
	return (node << 2) | base;
}

uint64_t node_l_edge(uint64_t node, DNA_MAP base)
{
	uint64_t int_base = base;
	return node | (int_base << (K*2-2));
}

bool is_duplex_edge(uint64_t n, uint64_t e) {
	bool lefty = false, righty = false;
	for (const auto b : bases) {
		righty |= (e == node_r_edge(n, b));
		lefty |= (e == node_l_edge(n, b));
	}
	return righty && lefty;
}

bool is_left_edge(uint64_t n, uint64_t e) {
	for (const auto b : bases)
		if (e == canonicalize(node_l_edge(n, b)))
			return !is_duplex_edge(n, e);
	return false;
}

bool is_right_edge(uint64_t n, uint64_t e) {
	for (const auto b : bases)
		if (e == canonicalize(node_r_edge(n, b)))
			return !is_duplex_edge(n, e);
	return false;
}

void record_end(uint64_t n, uint64_t e, flush_object *obj) {
	uint64_t hash = HashUtil::hash_64(n, BITMASK(2*(K-1)));
	if (is_left_edge(n, e)) {
		qf_insert(obj->first_qf, hash, 0, 1, true, true);
		//cout << "Begin node: " << int_to_str_node(n) << endl;
		//cout << "Begin edge: " << int_to_str(e) << endl;
	}
	else if (is_right_edge(n, e)) {
		qf_insert(obj->last_qf, hash, 0, 1, true, true);
		//cout << "End node: " << int_to_str_node(n) << endl;
		//cout << "End edge: " << int_to_str(e) << endl;
	}
}

/* convert a chunk of the fastq file into kmers */
void reads_to_kmers(chunk &c, flush_object *obj)
{
	auto fs = c.get_reads();
	auto fe = c.get_reads();
	auto end = fs + c.get_size();
	while (fs && fs!=end) {
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
		fs++; // increment the pointer

		fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
		string read(fs, fe-fs);
		/*cout << read << endl;*/

start_read:
		if (read.length() < K) // start with the next read if length is smaller than K
			goto next_read;
		{
			uint64_t first = 0;
			uint64_t first_rev = 0;
			uint64_t first_can = 0;
			uint64_t item = 0;
			//cout << "K " << read.substr(0,K) << endl;
			for(int i=0; i<K; i++) { //First kmer
				uint8_t curr = kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					read = read.substr(i+1, read.length());
					goto start_read;
				}
				first = first | curr;
				first = first << 2;
			}
			first = first >> 2;
			first_rev = kmer::reverse_complement(first);

			//cout << "kmer: "; cout << int_to_str(first);
			//cout << " reverse-comp: "; cout << int_to_str(first_rev) << endl;

			if (kmer::compare_kmers(first, first_rev))
				first_can = first;
			else
				first_can = first_rev;

			// create a set of true k-mers
			uint64_t exact_hash = HashUtil::hash_64(first_can, BITMASK(2*K));
			qf_insert(obj->exact_qf, exact_hash, 0, 1, true, true);
			//true_kmers.insert(first_can);

			// this is the first kmer of the read, hash using invertible hash and insert
			// in the first QF
			// No need to use the range, because by default the QF will use the
			// lower order 56 (q+r) bits.
			uint64_t node = kmer::prefix(first, K-1);
			node = canonicalize_node(node);
			record_end(node, first_can, obj);

			// hash the kmer using murmurhash/xxHash before adding to the list
			item = HashUtil::MurmurHash64A(((void*)&first_can), sizeof(first_can),
																		 obj->local_qf->metadata->seed);
			/*
			 * first try and insert in the main QF.
			 * If lock can't be accuired in the first attempt then
			 * insert the item in the local QF.
			 */
			if (!qf_insert(obj->main_qf, item%obj->main_qf->metadata->range, 0, 1,
										 true, false)) {
				qf_insert(obj->local_qf, item%obj->local_qf->metadata->range, 0, 1,
									false, false);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					dump_local_qf_to_main(obj);
					obj->count = 0;
				}
			}
			//cout<< "X " << bitset<64>(first)<<endl;

			uint64_t next = first;
			uint64_t next_rev = first_rev;
			uint64_t next_can = first_can;

			for(uint32_t i=K; i<read.length(); i++) { //next kmers
				uint64_t nextn = kmer::suffix(next, K-1);
				nextn = canonicalize_node(nextn);
				if (can_have_duplex_edges(nextn))
					record_end(nextn, next_can, obj);
				//cout << "K: " << read.substr(i-K+1,K) << endl;
				uint8_t curr = kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					read = read.substr(i+1, read.length());
					// don't insert if i = K and node is a duplex node.
					if (i > K || !can_have_duplex_edges(nextn)) {
						nextn = kmer::suffix(next, K-1);
						nextn = canonicalize_node(nextn);
						record_end(nextn, next_can, obj);
					}
					goto start_read;
				}
				next = (next << 2) & BITMASK(2*K);
				next_rev = next_rev >> 2;
				
				next |= curr;
				uint64_t tmp = kmer::reverse_complement_base(curr);
				tmp <<= (K*2-2);
				next_rev = next_rev | tmp;
				if (kmer::compare_kmers(next, next_rev))
					next_can = next;
				else
					next_can = next_rev;

				uint64_t exact_hash = HashUtil::hash_64(next_can, BITMASK(2*K));
				qf_insert(obj->exact_qf, exact_hash, 0, 1, true, true);
				//true_kmers.insert(next_can);

				// hash the kmer using murmurhash/xxHash before adding to the list
				item = HashUtil::MurmurHash64A(((void*)&next_can), sizeof(next_can),
																			 obj->local_qf->metadata->seed);
				//item = XXH63 (((void*)&item), sizeof(item), seed);

				/*
				 * first try and insert in the main QF.
				 * If lock can't be accuired in the first attempt then
				 * insert the item in the local QF.
				 */
				if (!qf_insert(obj->main_qf, item%obj->main_qf->metadata->range, 0, 1, true,
											 false)) {
					qf_insert(obj->local_qf, item%obj->local_qf->metadata->range, 0, 1, false,
										false);
					obj->count++;
					// check of the load factor of the local QF is more than 50%
					if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
						dump_local_qf_to_main(obj);
						obj->count = 0;
					}
				}

				//cout<<bitset<64>(next)<<endl;
				//assert(next == str_to_int(read.substr(i-K+1,K)));

				if (can_have_duplex_edges(nextn))
					record_end(nextn, next_can, obj);
			}
			node = kmer::suffix(next, K-1);
			node = canonicalize_node(node);
			record_end(node, next_can, obj);
		}

next_read:
		fs = ++fe;		// increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
		fs++; // increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
		fs++; // increment the pointer
	}
	free(c.get_reads());
}

/* read a part of the fastq file, parse it, convert the reads to kmers, and
 * insert them in the CQF
 */
static bool fastq_to_uint64kmers_prod(flush_object* obj) 
{
	file_pointer* fp;

	while (num_files) {
		while (ip_files.pop(fp)) {
			if (fastq_read_parts(fp->mode, fp)) {
				ip_files.push(fp);
				chunk c(fp->part, fp->size);
				reads_to_kmers(c, obj);
			} else {
				/* close the file */
				if (fp->mode == 0)
					fclose(fp->freader->in);
				else if (fp->mode == 1)
					gzclose(fp->freader->in_gzip);
				else if (fp->mode == 2)
					if (fp->freader->in) {
						BZ2_bzReadClose(&(fp->freader->bzerror), fp->freader->in_bzip2);
						fclose(fp->freader->in);
					}
				delete[] fp->part_buffer;
				delete fp;
				num_files--;
			}
		}
	}
	if (obj->count) {
		dump_local_qf_to_main(obj);
		obj->count = 0;
	}

		return true;
}

bool getFileReader(int mode, char* fastq_file, reader* file_reader) 
{
	uint64_t gzip_buffer_size = 1ULL << 26;
	uint64_t bzip2_buffer_size = 1ULL << 26;

	if (mode == 0) {
		if ((file_reader->in = fopen(fastq_file, "rb")) == NULL)
			return false;
	} else if (mode == 1) {
		if ((file_reader->in_gzip = gzopen(fastq_file, "rb")) == NULL)
			return false;
		gzbuffer(file_reader->in_gzip, gzip_buffer_size);
	} else if (mode == 2) {
		file_reader->in = fopen(fastq_file, "rb");
		if (!file_reader->in)
			return false;
		setvbuf(file_reader->in, NULL, _IOFBF, bzip2_buffer_size);
		if ((file_reader->in_bzip2 = BZ2_bzReadOpen(&file_reader->bzerror,
																								file_reader->in, 0, 0, NULL,
																								0)) == NULL) {
			fclose(file_reader->in);
			return false;
		}
	}
	return true;
}

/* main method */
int main(int argc, char *argv[])
{
	QF cf, last_cf, first_cf, exact_cf;
	QFi cfi;
	QF local_qfs[50];
	int mode = atoi(argv[1]);
	int qbits = atoi(argv[2]);
	int qbits_aux = atoi(argv[3]);
	int numthreads = atoi(argv[4]);
	int num_hash_bits = qbits+8;	  // we use 8 bits for remainders in the main QF
	int num_hash_bits_exact = 2*K;	  // for exact k-mer counting
	int num_hash_bits_aux = 2*(K-1);  // for exact counting.
	string first_ext(".first");
	string last_ext(".last");
	string exact_ext(".exact");
	string ser_ext(".ser");
	string log_ext(".log");
	string cluster_ext(".cluster");
	string freq_ext(".freq");
	struct timeval start1, start2, end1, end2;
	struct timezone tzp;
	uint32_t OVERHEAD_SIZE = 65535;

	for (int i = 5; i < argc; i++) {
		auto* fr = new reader;
		if (getFileReader(mode, argv[i], fr)) {
			file_pointer* fp = new file_pointer;
			fp->mode = mode;
			fp->freader.reset(fr);
			fp->part_buffer = new char[OVERHEAD_SIZE];
			ip_files.push(fp);
			num_files++;
		} else {
			delete fr;
		}
	}

	string ds_file = string(argv[5]) + ser_ext;
	string ds_file_first = string(argv[5]) + first_ext;
	string ds_file_last = string(argv[5]) + last_ext;
	string ds_file_exact = string(argv[5]) + exact_ext;
	string log_file = string(argv[5]) + log_ext;
	string cluster_file = string(argv[5]) + cluster_ext;
	//string freq_file = string(argv[5]) + freq_ext;
	
	uint32_t seed = time(NULL);
	//Initialize the main  QF and aux QFs
	qf_init(&cf, (1ULL<<qbits), num_hash_bits, 0, true, "", seed);
	qf_init(&exact_cf, (1ULL<<qbits), num_hash_bits_exact, 0, true, "", seed);
	qf_init(&first_cf, (1ULL<<qbits_aux), num_hash_bits_aux, 0, true, "", seed);
	qf_init(&last_cf, (1ULL<<qbits_aux), num_hash_bits_aux, 0, true, "", seed);

	boost::thread_group prod_threads;

	for (int i = 0; i < numthreads; i++) {
		qf_init(&local_qfs[i], (1ULL << QBITS_LOCAL_QF), num_hash_bits, 0, true,
						"", seed);
		flush_object* obj = (flush_object*)malloc(sizeof(flush_object));
		obj->local_qf = &local_qfs[i];
		obj->main_qf = &cf;
		obj->first_qf = &first_cf;
		obj->last_qf = &last_cf;
		obj->exact_qf = &exact_cf;
		prod_threads.add_thread(new boost::thread(fastq_to_uint64kmers_prod,
																							obj));
	}

	cout << "Reading from the fastq file and inserting in the QF" << endl;
	gettimeofday(&start1, &tzp);
	prod_threads.join_all();
	qf_serialize(&cf, ds_file.c_str());
	qf_serialize(&exact_cf, ds_file_exact.c_str());
	qf_serialize(&first_cf, ds_file_first.c_str());
	qf_serialize(&last_cf, ds_file_last.c_str());
	gettimeofday(&end1, &tzp);
	print_time_elapsed("", &start1, &end1);

	cout << "Calc freq distribution: " << endl;
	//ofstream freq_file;
	//freq_file.open(freq_file.c_str());
	uint64_t max_cnt = 0;
	qf_iterator(&cf, &cfi, 0);
	gettimeofday(&start2, &tzp);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&cfi, &key, &value, &count);
		//freq_file << key << " " << count << endl;
		if (max_cnt < count)
			max_cnt = count;
	} while (!qfi_next(&cfi));
	gettimeofday(&end2, &tzp);
	print_time_elapsed("", &start2, &end2);

	cout << "Maximum freq: " << max_cnt << endl;
	//freq_file.close();

	cout << "Num distinct elem: " << cf.metadata->ndistinct_elts << endl;
	cout << "Total num elems: " << cf.metadata->nelts << endl;

	//to validate the invertible hashing and exact QF
	qf_iterator(&exact_cf, &cfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&cfi, &key, &value, &count);
		//cout << int_to_str_node(HashUtil::hash_64i(key, BITMASK(2*K))) << endl;
	} while (!qfi_next(&cfi));
	cout << "Exact: Num distinct elem: " << exact_cf.metadata->ndistinct_elts << endl;
	cout << "Exact: Total num elems: " << exact_cf.metadata->nelts << endl;

	//to validate the invertible hashing and exact QF
	qf_iterator(&first_cf, &cfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&cfi, &key, &value, &count);
		//cout << int_to_str_node(HashUtil::hash_64i(key, BITMASK(2*K))) << endl;
	} while (!qfi_next(&cfi));
	cout << "First: Num distinct elem: " << first_cf.metadata->ndistinct_elts << endl;
	cout << "First: Total num elems: " << first_cf.metadata->nelts << endl;

	//to validate the invertible hashing and exact QF
	qf_iterator(&last_cf, &cfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&cfi, &key, &value, &count);
		//cout << int_to_str_node(HashUtil::hash_64i(key, BITMASK(2*K))) << endl;
	} while (!qfi_next(&cfi));
	cout << "Last: Num distinct elem: " << last_cf.metadata->ndistinct_elts << endl;
	cout << "Last: Total num elems: " << last_cf.metadata->nelts << endl;

#ifdef LOG_WAIT_TIME
	ofstream wait_time_log;
	wait_time_log.open(log_file.c_str());
	wait_time_log << "Id\tTotalTimeSingle\tTotalTimeSpinning\tNumLocks\tNumSingleAttempt\tPercentageSpinningTimes"
		<< endl;
	for (uint32_t i=0; i<cf.num_locks; i++)
		wait_time_log << i << "\t" << cf.wait_times[i].total_time_single << "\t\t\t" 
			<< cf.wait_times[i].total_time_spinning << "\t\t\t" 
			<< cf.wait_times[i].locks_taken << "\t\t\t" 
			<< cf.wait_times[i].locks_acquired_single_attempt << "\t\t\t"
			<< ((double)(cf.wait_times[i].locks_taken
									 -cf.wait_times[i].locks_acquired_single_attempt)
					/(double)cf.wait_times[i].locks_taken)*100
			<< endl;
	wait_time_log.close();
#endif

#ifdef LOG_CLUSTER_LENGTH
	ofstream cluster_len_log;
	cluster_len_log.open(cluster_file.c_str());
	cluster_len_log << "StartingIndex\tLength" << endl;
	for (uint32_t i = 0; i < cfi.num_clusters; i++)
		cluster_len_log << cfi.c_info[i].start_index << "\t\t" << 
			cfi.c_info[i].length << endl;
	cluster_len_log.close();

#endif

	//destroy the QF and reclaim the memory
	qf_destroy(&cf, true);
	qf_destroy(&first_cf, true);
	qf_destroy(&last_cf, true);

	return 0;
}
