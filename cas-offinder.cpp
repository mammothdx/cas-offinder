#include "config.h"

#include "oclkernels.h"
#include "cas-offinder.h"

#ifndef MIN
#define MIN(a,b) ( ((a)<(b))?(a):(b) )
#endif
#ifndef MAX
#define MAX(a,b) ( ((a)>(b))?(a):(b) )
#endif
#include <sstream>
#include <iterator>

using namespace std;

#define add_compare(a, b, c)\
	m_compares[a].push_back(make_pair(b, c));

static inline string trim_left_bases(const string& seq, size_t trim) {
	return (trim >= seq.size()) ? string() : seq.substr(trim);
}

static inline string trim_right_bases(const string& seq, size_t trim) {
	return (trim >= seq.size()) ? string() : seq.substr(0, seq.size() - trim);
}

static inline string trim_search_padding(const string& seq, size_t trim, bool trim_right) {
	return trim_right ? trim_right_bases(seq, trim) : trim_left_bases(seq, trim);
}

static inline char complement_base(char base) {
	switch (base) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'G': return 'C';
		case 'C': return 'G';
		case 'R': return 'Y';
		case 'Y': return 'R';
		case 'M': return 'K';
		case 'K': return 'M';
		case 'H': return 'D';
		case 'D': return 'H';
		case 'B': return 'V';
		case 'V': return 'B';
		default: return base;
	}
}

static inline string reverse_complement_string(const string& seq) {
	string rc(seq.rbegin(), seq.rend());
	transform(rc.begin(), rc.end(), rc.begin(), complement_base);
	return rc;
}

// Report the forward-strand genome index of the first base in the emitted DNA
// slice after trimming hidden search padding in guide orientation.
static inline int output_alignment_start_offset(size_t trim, bool trim_right, char direction) {
	if (direction == '-')
		return (int)(trim_right ? trim : 0);
	return (int)(trim_right ? 0 : trim);
}

vector<string> split(string const &input) {
	istringstream sbuffer(input);
	vector<string> ret((istream_iterator<string>(sbuffer)), istream_iterator<string>());
	return ret;
}

vector<string> split(string const &input, char delim) {
	istringstream sbuffer(input);
	vector<string> ret;
	string item;
	while (getline(sbuffer, item, delim))
		ret.push_back(item);
	return ret;
}

static inline bool iupac_contains_base(char pattern, char base) {
	base = (char)toupper(base);
	switch ((char)toupper(pattern)) {
		case 'A': return base == 'A';
		case 'C': return base == 'C';
		case 'G': return base == 'G';
		case 'T': return base == 'T';
		case 'R': return base == 'A' || base == 'G';
		case 'Y': return base == 'C' || base == 'T';
		case 'K': return base == 'G' || base == 'T';
		case 'M': return base == 'A' || base == 'C';
		case 'W': return base == 'A' || base == 'T';
		case 'S': return base == 'C' || base == 'G';
		case 'H': return base == 'A' || base == 'C' || base == 'T';
		case 'B': return base == 'C' || base == 'G' || base == 'T';
		case 'V': return base == 'A' || base == 'C' || base == 'G';
		case 'D': return base == 'A' || base == 'G' || base == 'T';
		case 'N': return base == 'A' || base == 'C' || base == 'G' || base == 'T';
		default: return false;
	}
}

void Cas_OFFinder::set_complementary_sequence(cl_char* seq, size_t seqlen) {
	size_t i, l = 0;
	cl_char tmp;

	for (i = 0; i < seqlen; i++) {
		if (seq[i] == 'A') seq[i] = 'T';
		else if (seq[i] == 'T') seq[i] = 'A';
		else if (seq[i] == 'G') seq[i] = 'C';
		else if (seq[i] == 'C') seq[i] = 'G';
		else if (seq[i] == 'R') seq[i] = 'Y';
		else if (seq[i] == 'Y') seq[i] = 'R';
		else if (seq[i] == 'M') seq[i] = 'K';
		else if (seq[i] == 'K') seq[i] = 'M';
		else if (seq[i] == 'H') seq[i] = 'D';
		else if (seq[i] == 'D') seq[i] = 'H';
		else if (seq[i] == 'B') seq[i] = 'V';
		else if (seq[i] == 'V') seq[i] = 'B';
	}
	for (i = 0; i < seqlen / 2; i++) {
		tmp = seq[i];
		seq[i] = seq[seqlen - i - 1];
		seq[seqlen - i - 1] = tmp;
	}
}

void Cas_OFFinder::set_seq_flags(int* seq_flags, const cl_char* seq, size_t seqlen) {
	int i, n = 0;
	for (i = 0; i < seqlen; i++) {
		if (seq[i] != 'N') {
			seq_flags[n] = i;
			n++;
		}
	}
	if (i != n)
		seq_flags[n] = -1;
}

void Cas_OFFinder::initOpenCLPlatforms() {
	oclGetPlatformIDs(MAX_PLATFORM_NUM, platforms, &platform_cnt);
	if (platform_cnt == 0) {
		cerr << "No OpenCL platforms found. Check OpenCL installation!" << endl;
		exit(1);
	}
}

void Cas_OFFinder::initOpenCLDevices(vector<unsigned int> dev_ids) {
	unsigned int i, j;

	cl_device_id* found_devices = new cl_device_id[MAX_DEVICE_NUM];
	cl_uint device_cnt;

	unsigned int dev_id = 0;
	vector<cl_device_id> devices;

	for (i = 0; i < platform_cnt; i++) {
		oclGetDeviceIDs(platforms[i], m_devtype, MAX_DEVICE_NUM, found_devices, &device_cnt);

		for (j = 0; j < device_cnt; j++) {
			if (dev_ids.size() == 0 || (dev_ids.size() > 0 && find(dev_ids.begin(), dev_ids.end(), dev_id) != dev_ids.end()))
				devices.push_back(found_devices[j]);
			dev_id += 1;
		}
	}

	m_devnum = devices.size();

	if (m_devnum == 0) {
		cerr << "No OpenCL devices found." << endl;
		exit(1);
	}

	cl_context context;
	cl_program program;

	const size_t src_len = strlen(program_src);
	for (i = 0; i < m_devnum; i++) {
		// Create completely separate contexts per device to avoid unknown errors
		context = oclCreateContext(0, 1, &devices[i], 0, 0);
		m_contexts.push_back(context);
		program = oclCreateProgramWithSource(context, 1, &program_src, &src_len);
		oclBuildProgram(program, 1, &devices[i], "", 0, 0);
				if (m_devtype == CL_DEVICE_TYPE_CPU) {
		m_finderkernels.push_back(oclCreateKernel(program, "finder_cpu"));
		m_comparerkernels.push_back(oclCreateKernel(program, "comparer_cpu"));
	} else {
		m_finderkernels.push_back(oclCreateKernel(program, "finder"));
		m_comparerkernels.push_back(oclCreateKernel(program, "comparer"));
	}
		m_queues.push_back(oclCreateCommandQueue(m_contexts[i], devices[i], 0));
		MAX_ALLOC_MEMORY.push_back(0);
		oclGetDeviceInfo(devices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &MAX_ALLOC_MEMORY[i], 0);
	}
	delete[] found_devices;
	cerr << "Total " << m_devnum << " device(s) found." << endl;
}

Cas_OFFinder::Cas_OFFinder(cl_device_type devtype, string devarg) {
	unsigned int i, j;
	int step;
	vector<unsigned int> dev_ids;

	m_devtype = devtype;
	vector<string> id_args = split(devarg, ',');
	vector<string> id_indices;
	for (i = 0; i < id_args.size(); i++) {
		id_indices = split(id_args[i], ':');
		if (id_indices.size() == 1) {
			dev_ids.push_back(atoi(id_indices[0].c_str()));
		}
		else if (id_indices.size() == 2 || id_indices.size() == 3) {
			step = 1;
			if (id_indices.size() == 3) step = atoi(id_indices[2].c_str());
			for (j = (unsigned int)atoi(id_indices[0].c_str()); j < (unsigned int)atoi(id_indices[1].c_str()); j += step) {
				dev_ids.push_back(j);
			}
		}
		else {
			cerr << "Something wrong with the device ID argument. Use all available devices instead..." << endl;
		}
	}
	initOpenCLDevices(dev_ids);
}

Cas_OFFinder::~Cas_OFFinder() {
	unsigned int i;
	for (i = 0; i < m_finderkernels.size(); i++)
		oclReleaseKernel(m_finderkernels[i]);
	for (i = 0; i < m_comparerkernels.size(); i++)
		oclReleaseKernel(m_comparerkernels[i]);
	for (i = 0; i < m_devnum; i++) {
		oclReleaseCommandQueue(m_queues[i]);
		oclReleaseContext(m_contexts[i]);
	}
	clearbufvec(&m_patternbufs);
	clearbufvec(&m_patternflagbufs);
	clearbufvec(&m_comparebufs);
	clearbufvec(&m_compareflagbufs);
	clearbufvec(&m_entrycountbufs);
}

void Cas_OFFinder::setChrData() {
	unsigned int dev_index;

	m_chrdatasize = m_chrdata.size();
	m_totalanalyzedsize = 0;
	m_lasttotalanalyzedsize = 0;
	m_lastloci = 0;

	m_dicesizes.clear();
	clearbufvec(&m_chrdatabufs);
	clearbufvec(&m_flagbufs);
	clearbufvec(&m_locibufs);

	for (dev_index = 0; dev_index < m_devnum; dev_index++) {
		m_dicesizes.push_back(
			(size_t)MIN(
				(MAX_ALLOC_MEMORY[dev_index] - sizeof(cl_char) * (3 * m_patternlen - 1) - sizeof(cl_uint) * (2 * m_patternlen + 3) - sizeof(cl_ushort)) / (4 * sizeof(cl_char) + 3 * sizeof(cl_uint) + 2 * sizeof(cl_ushort)),
				((m_chrdatasize / m_devnum) + ((m_chrdatasize%m_devnum == 0) ? 0 : 1))
			)
		); // No more than maximum allocation per device
		// cerr << "Dicesize: " << m_dicesizes[dev_index] << endl;
		m_chrdatabufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char) * (m_dicesizes[dev_index] + m_patternlen - 1), 0));
		m_flagbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_char) * m_dicesizes[dev_index], 0));
		m_locibufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_uint) * m_dicesizes[dev_index], 0));

		oclSetKernelArg(m_finderkernels[dev_index], 0, sizeof(cl_mem), &m_chrdatabufs[dev_index]);
		oclSetKernelArg(m_finderkernels[dev_index], 4, sizeof(cl_mem), &m_flagbufs[dev_index]);
		oclSetKernelArg(m_finderkernels[dev_index], 6, sizeof(cl_mem), &m_locibufs[dev_index]);

		oclSetKernelArg(m_comparerkernels[dev_index], 0, sizeof(cl_mem), &m_chrdatabufs[dev_index]);
		oclSetKernelArg(m_comparerkernels[dev_index], 1, sizeof(cl_mem), &m_locibufs[dev_index]);
		oclSetKernelArg(m_comparerkernels[dev_index], 7, sizeof(cl_mem), &m_flagbufs[dev_index]);
	}
}

bool Cas_OFFinder::loadNextChunk() {
	if (m_totalanalyzedsize == m_chrdatasize)
		return false;

	unsigned int dev_index;
	unsigned long long tailsize;
	size_t overlap;

	m_activedevnum = 0;
	m_worksizes.clear();
	m_lasttotalanalyzedsize = m_totalanalyzedsize;

	for (dev_index = 0; dev_index < m_devnum; dev_index++) {
		tailsize = m_chrdatasize - m_totalanalyzedsize;
		m_activedevnum++;
		if (tailsize <= m_dicesizes[dev_index]) {
			if (tailsize < m_patternlen) {
				m_totalanalyzedsize += tailsize;
				m_worksizes.push_back(0);
				break;
			}
			oclEnqueueWriteBuffer(m_queues[dev_index], m_chrdatabufs[dev_index], CL_TRUE, 0, sizeof(cl_char) * tailsize, (cl_char *)m_chrdata.c_str() + m_totalanalyzedsize, 0, 0, 0);
			m_totalanalyzedsize += tailsize;
			m_worksizes.push_back(tailsize - m_patternlen + 1);
#ifdef DEBUG
			cerr << "Worksize: " << m_worksizes[dev_index] << ", Tailsize: " << tailsize << endl;
#endif
			break;
		}
		else {
			overlap = MIN(m_patternlen - 1, tailsize - m_dicesizes[dev_index]);
			oclEnqueueWriteBuffer(m_queues[dev_index], m_chrdatabufs[dev_index], CL_TRUE, 0, sizeof(cl_char) * (m_dicesizes[dev_index] + overlap), (cl_char *)m_chrdata.c_str() + m_totalanalyzedsize, 0, 0, 0);
			m_totalanalyzedsize += m_dicesizes[dev_index];
			m_worksizes.push_back(m_dicesizes[dev_index] - m_patternlen + overlap + 1);
#ifdef DEBUG
			cerr << "Worksize: " << m_worksizes[dev_index] << ", Tailsize: " << tailsize << endl;
#endif
		}
	}
	cerr << m_activedevnum << " devices selected to analyze..." << endl;

	return true;
}

void Cas_OFFinder::findPattern() {
	unsigned int dev_index;
	cl_uint zero = 0;
	for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
		const size_t worksize = (size_t)m_worksizes[dev_index];
		if (worksize == 0)
			continue;
		oclEnqueueWriteBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &zero, 0, 0, 0);
		oclEnqueueNDRangeKernel(m_queues[dev_index], m_finderkernels[dev_index], 1, 0, &worksize, 0, 0, 0, 0);
	}

	for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
		m_locicnts.push_back(0);
		if (m_worksizes[dev_index] > 0) {
			oclFinish(m_queues[dev_index]);
			oclEnqueueReadBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &m_locicnts[dev_index], 0, 0, 0);
		}
		if (m_locicnts[dev_index] > 0) {
			m_flags.push_back((cl_char *)malloc(sizeof(cl_char) * m_locicnts[dev_index]));
			oclEnqueueReadBuffer(m_queues[dev_index], m_flagbufs[dev_index], CL_TRUE, 0, sizeof(cl_char) * m_locicnts[dev_index], m_flags[dev_index], 0, 0, 0);

			m_mmcounts.push_back((cl_ushort *)malloc(sizeof(cl_ushort) * m_locicnts[dev_index] * 2)); // Maximum numbers of mismatch counts
			m_directions.push_back((cl_char *)malloc(sizeof(cl_char) * m_locicnts[dev_index] * 2));
			m_mmlocis.push_back((cl_uint *)malloc(sizeof(cl_uint) * m_locicnts[dev_index] * 2));

			m_mmlocibufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_uint) * m_locicnts[dev_index] * 2, 0));
			m_mmcountbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_ushort) * m_locicnts[dev_index] * 2, 0));
			m_directionbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_WRITE_ONLY, sizeof(cl_char) * m_locicnts[dev_index] * 2, 0));

			oclSetKernelArg(m_comparerkernels[dev_index], 2, sizeof(cl_mem), &m_mmlocibufs[dev_index]);
			oclSetKernelArg(m_comparerkernels[dev_index], 8, sizeof(cl_mem), &m_mmcountbufs[dev_index]);
			oclSetKernelArg(m_comparerkernels[dev_index], 9, sizeof(cl_mem), &m_directionbufs[dev_index]);
		}
		else {
			m_flags.push_back(0);
			m_mmcounts.push_back(0);
			m_directions.push_back(0);
			m_mmlocis.push_back(0);
			m_mmlocibufs.push_back(0);
			m_mmcountbufs.push_back(0);
			m_directionbufs.push_back(0);
		}
	}
}

void Cas_OFFinder::indicate_mismatches(cl_char* seq, cl_char* comp) {
	unsigned int k;
	for (k = 0; k < m_patternlen; k++)
		if ((comp[k] == 'R' && (seq[k] == 'C' || seq[k] == 'T')) ||
			(comp[k] == 'Y' && (seq[k] == 'A' || seq[k] == 'G')) ||
			(comp[k] == 'K' && (seq[k] == 'A' || seq[k] == 'C')) ||
			(comp[k] == 'M' && (seq[k] == 'G' || seq[k] == 'T')) ||
			(comp[k] == 'W' && (seq[k] == 'C' || seq[k] == 'G')) ||
			(comp[k] == 'S' && (seq[k] == 'A' || seq[k] == 'T')) ||
			(comp[k] == 'H' && (seq[k] == 'G')) ||
			(comp[k] == 'B' && (seq[k] == 'A')) ||
			(comp[k] == 'V' && (seq[k] == 'T')) ||
			(comp[k] == 'D' && (seq[k] == 'C')) ||
			(comp[k] == 'A' && (seq[k] != 'A')) ||
			(comp[k] == 'G' && (seq[k] != 'G')) ||
			(comp[k] == 'C' && (seq[k] != 'C')) ||
			(comp[k] == 'T' && (seq[k] != 'T')))
			seq[k] += 32;
}

void Cas_OFFinder::writeHeaders(const char* outfilename) {
	ostream *fo;
	if (strlen(outfilename) == 1 && outfilename[0] == '-') {
		fo = &cout;
	} else {
		fo = new ofstream(outfilename, ios::out | ios::trunc);
	}
	(*fo) << "##Generated by Cas-OFFinder " << CAS_OFFINDER_VERSION << endl;
	(*fo) << "#Id\tBulge Type\tcrRNA\tDNA\tChromosome\tLocation\tDirection\tMismatches\tBulge Size" << endl;
}

void Cas_OFFinder::compareAll(const char* outfilename, bool issummary) {
	unsigned int i, j, dev_index;
	unsigned int bulge_size, bulge_index;
	bool guide_reversed_pam;
	bool trim_right;
	cl_ushort threshold;
	cl_ushort max_threshold;
	string compare;
	string seq_rna;
	string seq_dna;
	string id;
	string bulge_type;
	int offset;
	unsigned long long loci;
	cl_uint zero = 0;

	cl_char *cl_compare = new cl_char[m_patternlen * 2];
	cl_int *cl_compare_flags = new cl_int[m_patternlen * 2];

	char *strbuf = new char[m_patternlen + 1]; strbuf[m_patternlen] = 0;
	bool isfile = false;

	ostream *fo;
	if (strlen(outfilename) == 1 && outfilename[0] == '-') {
		fo = &cout;
	} else {
		fo = new ofstream(outfilename, ios::out | ios::app);
		isfile = true;
	}
	for (auto &ci: m_compares) {
		compare = ci.first;
		max_threshold = 0;
		for (const auto& entry : ci.second)
			max_threshold = MAX(max_threshold, entry.first.second);
		memcpy(cl_compare, compare.c_str(), m_patternlen);
		memcpy(cl_compare + m_patternlen, compare.c_str(), m_patternlen);
		set_complementary_sequence(cl_compare + m_patternlen, m_patternlen);
		set_seq_flags(cl_compare_flags, cl_compare, m_patternlen);
		set_seq_flags(cl_compare_flags + m_patternlen, cl_compare + m_patternlen, m_patternlen);

		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			if (m_locicnts[dev_index] > 0) {
				oclEnqueueWriteBuffer(m_queues[dev_index], m_comparebufs[dev_index], CL_FALSE, 0, sizeof(cl_char) * m_patternlen * 2, cl_compare, 0, 0, 0);
				oclEnqueueWriteBuffer(m_queues[dev_index], m_compareflagbufs[dev_index], CL_FALSE, 0, sizeof(cl_int) * m_patternlen * 2, cl_compare_flags, 0, 0, 0);
				oclEnqueueWriteBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_FALSE, 0, sizeof(cl_uint), &zero, 0, 0, 0);
				oclFinish(m_queues[dev_index]);
				oclSetKernelArg(m_comparerkernels[dev_index], 6, sizeof(cl_ushort), &max_threshold);
				const size_t locicnts = m_locicnts[dev_index];
				oclEnqueueNDRangeKernel(m_queues[dev_index], m_comparerkernels[dev_index], 1, 0, &locicnts, 0, 0, 0, 0);
			}
		}

		unsigned long long localanalyzedsize = 0;
		unsigned int cnt = 0;
		unsigned int idx;
		char key[10000];
		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			if (m_locicnts[dev_index] > 0) {
				oclFinish(m_queues[dev_index]);
				oclEnqueueReadBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_TRUE, 0, sizeof(cl_uint), &cnt, 0, 0, 0);
				if (cnt > 0) {
					oclEnqueueReadBuffer(m_queues[dev_index], m_mmcountbufs[dev_index], CL_FALSE, 0, sizeof(cl_ushort) * cnt, m_mmcounts[dev_index], 0, 0, 0);
					oclEnqueueReadBuffer(m_queues[dev_index], m_directionbufs[dev_index], CL_FALSE, 0, sizeof(cl_char) * cnt, m_directions[dev_index], 0, 0, 0);
					oclEnqueueReadBuffer(m_queues[dev_index], m_mmlocibufs[dev_index], CL_FALSE, 0, sizeof(cl_uint) * cnt, m_mmlocis[dev_index], 0, 0, 0);
					oclFinish(m_queues[dev_index]);
					for (i = 0; i < cnt; i++) {
						loci = m_mmlocis[dev_index][i] + m_lasttotalanalyzedsize + localanalyzedsize;
						if (m_mmcounts[dev_index][i] <= max_threshold) {
							strncpy(strbuf, (char *)(m_chrdata.c_str() + loci), m_patternlen);
							if (m_directions[dev_index][i] == '-') set_complementary_sequence((cl_char *)strbuf, m_patternlen);
							indicate_mismatches((cl_char*)strbuf, (cl_char*)compare.c_str());
							for (j = 0; ((j < m_chrpos.size()) && (loci >= m_chrpos[j])); j++) idx = j;
							for (const auto& entry : ci.second) {
								const bulgeinfo& bi = entry.second;
								threshold = entry.first.second;
								if (m_mmcounts[dev_index][i] > threshold)
									continue;
								id = entry.first.first;
								seq_dna = string(strbuf);
								guide_reversed_pam = bulge_is_reversed_pam(bi.second);
								// The reverse-complement step already puts seq_dna into guide
								// orientation, so hidden search padding stays on the guide
								// side regardless of strand.
								trim_right = guide_reversed_pam;
								bulge_index = decode_bulge_index(bi.second);
								const string& guide = m_guides[id];
								if (isnumeric(bi.first)) {
									bulge_size = (unsigned int)stoi(bi.first);
									offset = output_alignment_start_offset(m_dnabulgesize - bulge_size, trim_right, m_directions[dev_index][i]);
									seq_rna = guide;
									if (bulge_size > 0)
										seq_rna = guide.substr(0, bulge_index) + string(bulge_size, '-') + guide.substr(bulge_index);
									seq_dna = trim_search_padding(seq_dna, m_dnabulgesize - bulge_size, trim_right);
									bulge_type = (bulge_size == 0) ? "X" : "DNA";
								} else {
									bulge_size = (unsigned int)bi.first.size();
									offset = output_alignment_start_offset(m_dnabulgesize + bulge_size, trim_right, m_directions[dev_index][i]);
									seq_rna = guide;
									seq_dna = trim_search_padding(seq_dna, m_dnabulgesize + bulge_size, trim_right);
									seq_dna = seq_dna.substr(0, bulge_index) + string(bulge_size, '-') + seq_dna.substr(bulge_index);
									bulge_type = "RNA";
								}
								(*fo) << id << "\t" << bulge_type << "\t" << seq_rna << "\t" << seq_dna << "\t" << m_chrnames[idx] << "\t" << loci - m_chrpos[idx] + offset << "\t" << m_directions[dev_index][i] << "\t" << m_mmcounts[dev_index][i] << "\t" << bulge_size << endl;
								if (issummary) {
									snprintf(key, 10000, "%s,%s,%d,%d", id.c_str(), bulge_type.c_str(), bulge_size, m_mmcounts[dev_index][i]);
									++m_summarytable[string(key)];
								}
							}
						}
					}
				}
			}

			localanalyzedsize += m_worksizes[dev_index];
		}
		fo->flush();
	}

	if (m_has_left_pam && m_has_right_pam && m_rnabulgesize > 0) {
		unsigned long long chunk_offset = 0;
		unsigned int idx;
		for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
			for (unsigned long long pos = 0; pos < m_worksizes[dev_index]; pos++) {
				loci = m_lasttotalanalyzedsize + chunk_offset + pos;
				strncpy(strbuf, (char *)(m_chrdata.c_str() + loci), m_patternlen);
				set_complementary_sequence((cl_char *)strbuf, m_patternlen);
				for (j = 0; ((j < m_chrpos.size()) && (loci >= m_chrpos[j])); j++) idx = j;
				for (const auto& guide_entry : m_guides) {
					const string& guide = guide_entry.second;
					const cl_ushort guide_threshold = m_thresholds[guide_entry.first];
					const unsigned int first_non_n = (unsigned int)guide.find_first_not_of('N');
					const unsigned int last_non_n = (unsigned int)guide.find_last_not_of('N');
					for (unsigned int bulge_size_fallback = 1; bulge_size_fallback <= m_rnabulgesize; bulge_size_fallback++) {
						const string dna_trimmed = trim_search_padding(string(strbuf), m_dnabulgesize + bulge_size_fallback, true);
						for (unsigned int bulge_index_fallback = first_non_n; bulge_index_fallback + bulge_size_fallback <= last_non_n + 1; bulge_index_fallback++) {
							string seq_dna_fallback = dna_trimmed.substr(0, bulge_index_fallback) + string(bulge_size_fallback, '-') + dna_trimmed.substr(bulge_index_fallback);
							cl_ushort mismatch_count_fallback = 0;
							bool valid = true;
							for (unsigned int col = 0; col < guide.size(); col++) {
								const char dna_char = seq_dna_fallback[col];
								if (dna_char == '-')
									continue;
								const char dna_base = (char)toupper(dna_char);
								if (!iupac_contains_base(m_input_pattern[col], dna_base)) {
									valid = false;
									break;
								}
								if (!iupac_contains_base(guide[col], dna_base)) {
									mismatch_count_fallback++;
									seq_dna_fallback[col] = (char)tolower(dna_char);
									if (mismatch_count_fallback > guide_threshold) {
										valid = false;
										break;
									}
								}
							}
							if (!valid)
								continue;
							const int offset_fallback = output_alignment_start_offset(m_dnabulgesize + bulge_size_fallback, true, '-');
							(*fo) << guide_entry.first << "\tRNA\t" << guide << "\t" << seq_dna_fallback << "\t"
							      << m_chrnames[idx] << "\t" << loci - m_chrpos[idx] + offset_fallback << "\t-\t"
							      << mismatch_count_fallback << "\t" << bulge_size_fallback << endl;
						}
					}
				}
			}
			chunk_offset += m_worksizes[dev_index];
		}
		fo->flush();
	}

	if (isfile)
		((ofstream *)fo)->close();
	delete [] strbuf;
	delete [] cl_compare;
	delete [] cl_compare_flags;
}

void Cas_OFFinder::writeSummaryTable(const char* summaryfilename) {
	bool isfile = false;
	ostream *fo;
	int cnt;
	size_t pos;
	string key, token;
	vector<string> summaryheaders = {
		"Id", "Bulge Type", "Bulge Size", "Mismatches"
	};
	if (strlen(summaryfilename) == 1 && summaryfilename[0] == '-') {
		fo = &cout;
	} else {
		fo = new ofstream(summaryfilename, ios::out);
		isfile = true;
	}
		for (const auto& kv : m_summarytable) {
		cnt = 0;
		key = string(kv.first);
		(*fo) << "##";
		while ((pos = key.find(",")) != string::npos) {
			token = key.substr(0, pos);
			(*fo) << summaryheaders[cnt++] << "=" << token << ";";
			key.erase(0, pos + 1);
		}
		(*fo) << summaryheaders[cnt++] << "=" << key << ";";
				(*fo) << "Number of Found Targets=" << kv.second << endl;
		}
	if (isfile)
		((ofstream *)fo)->close();
}

void Cas_OFFinder::releaseLociinfo() {
	unsigned int dev_index;

	for (dev_index = 0; dev_index < m_activedevnum; dev_index++) {
		free((void *)m_mmcounts[dev_index]);
		free((void *)m_flags[dev_index]);
		free((void *)m_directions[dev_index]);
		free((void *)m_mmlocis[dev_index]);
	}
	m_directions.clear();
	m_mmlocis.clear();
	m_mmcounts.clear();
	m_locicnts.clear();
	clearbufvec(&m_mmlocibufs);
	m_flags.clear();
	clearbufvec(&m_mmcountbufs);
	clearbufvec(&m_directionbufs);
}

void Cas_OFFinder::print_usage() {
	unsigned int i, j;
	cout << "Cas-OFFinder " << CAS_OFFINDER_VERSION << endl <<
		endl <<
		"Copyright (c) 2021 Jeongbin Park and Sangsu Bae" << endl <<
		"Website: " << CAS_OFFINDER_HOMEPAGE_URL << endl <<
		endl <<
		"Usage: cas-offinder [options] {input_filename|-} {C|G|A}[device_id(s)] {output_filename|-}" << endl <<
		"(C: using CPUs, G: using GPUs, A: using accelerators)" << endl
		<< endl <<
		"Options" << endl <<
		"  --summary <file>				 Print summary table to the specified file." << endl
		<< endl <<
		"Example input file (DNA bulge 2, RNA bulge 1):" << endl <<
		"/var/chromosomes/human_grch38" << endl <<
		"NNNNNNNNNNNNNNNNNNNNNRG 2 1" << endl <<
		"GGCCGACCTGTCGCTGACGCNNN 5" << endl <<
		"CGCCAGCGTCAGCGACAGGTNNN 5" << endl <<
		"ACGGCGCCAGCGTCAGCGACNNN 5" << endl <<
		"GTCGCTGACGCTGGCGCCGTNNN 5" << endl <<
		endl <<
		"Available device list:" << endl;

	cl_device_id devices_per_platform[MAX_DEVICE_NUM];
	cl_uint device_cnt;
	cl_char devname[255] = { 0, };
	cl_char platformname[255] = { 0, };

	unsigned int cpu_id = 0, gpu_id = 0, acc_id = 0;
	for (i = 0; i < platform_cnt; i++) {
		oclGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 255, &platformname, 0);
		oclGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_CPU, MAX_DEVICE_NUM, devices_per_platform, &device_cnt);
		for (j = 0; j < device_cnt; j++) {
			oclGetDeviceInfo(devices_per_platform[j], CL_DEVICE_NAME, 255, &devname, 0);
			cout << "Type: CPU, ID: " << cpu_id++ << ", <" << devname << "> on <" << platformname << ">" << endl;
		}
		oclGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, MAX_DEVICE_NUM, devices_per_platform, &device_cnt);
		for (j = 0; j < device_cnt; j++) {
			oclGetDeviceInfo(devices_per_platform[j], CL_DEVICE_NAME, 255, &devname, 0);
			cout << "Type: GPU, ID: " << gpu_id++ << ", <" << devname << "> on <" << platformname << ">" << endl;
		}
		oclGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ACCELERATOR, MAX_DEVICE_NUM, devices_per_platform, &device_cnt);
		for (j = 0; j < device_cnt; j++) {
			oclGetDeviceInfo(devices_per_platform[j], CL_DEVICE_NAME, 255, &devname, 0);
			cout << "Type: ACCELERATOR, ID: " << acc_id++ << ", <" << devname << "> on <" << platformname << ">" << endl;
		}
	}
}

static inline bool parseLine(istream& input, string& line, bool expect_line = true) {
	if (expect_line && input.eof()) {
		throw runtime_error("Unexpected end of file.");
	} else {
		if (!getline(input, line))
			return false;

		// Remove possible Windows-only '\r'.
		if (line[line.size()-1] == '\r')
			line = line.substr(0, line.size()-1);
	}
	return true;
}

/*
 * Reads the input file and configures this Cas_OFFinder instance.
 *
 * Sets public members: m_chrdir
 * Sets private members: m_dnabulgesize, m_rnabulgesize, m_pattern
 */
void Cas_OFFinder::parseInput(istream& input) {
	string line;
	vector<string> sline;
	unsigned int i, j, preNcnt;
	int threshold;
	string id = "";
	compareinfo ci;
	bulgeinfo bi;
	bool is_reversed_pam = false;
	bool has_right_pam = false;
	string compare, tmp;

	try {
		if (!input.good())
			throw runtime_error("Cannot read input file.");

		parseLine(input, m_chrdir);
		parseLine(input, line);
		sline = split(line);
		if (sline.size() != 1 && sline.size() != 3)
			throw runtime_error("Malformed input file.");
		
		m_dnabulgesize = 0;
		m_rnabulgesize = 0;
		m_has_left_pam = false;
		m_has_right_pam = false;

		if (sline.size() == 3) {
			m_dnabulgesize = atoi(sline[1].c_str());
			m_rnabulgesize = atoi(sline[2].c_str());
		}
		m_input_pattern = sline[0];
		has_right_pam = (sline[0][sline[0].size() - 1] != 'N');
		m_has_left_pam = (sline[0][0] != 'N');
		m_has_right_pam = has_right_pam;
		if (sline[0][0] != 'N') { // As of today (2021.07.29), there is no reversed PAM starts with 'N'
			is_reversed_pam = true;
			m_pattern = sline[0] + string(m_dnabulgesize, 'N');
		} else {
			m_pattern = string(m_dnabulgesize, 'N') + sline[0];
		}
		transform(m_pattern.begin(), m_pattern.end(), m_pattern.begin(), ::toupper);

		size_t linecnt = 0;
		size_t entrycnt = 0;
		while (parseLine(input, line, false)) {
			if (line.empty())
				break;
			sline = split(line);
			if (sline.size() != 2 && sline.size() != 3) {
				throw runtime_error("Malformed input file.");
				break;
			}
			if (sline[0].size() + m_dnabulgesize != m_pattern.size()) {
				throw runtime_error("The length of target sequences should match with the length of pattern sequence.");
			}
			transform(sline[0].begin(), sline[0].end(), sline[0].begin(), ::toupper);
			threshold = atoi(sline[1].c_str());
			if (sline.size() == 3)
				id = sline[2];
			else
				id = to_string(linecnt);
			m_guides[id] = sline[0];
			m_thresholds[id] = (cl_ushort)threshold;
			ci = make_pair(id, threshold);
			if (m_dnabulgesize == 0) {
				bi = make_pair("0", 0);
				add_compare(sline[0], ci, bi);
			} else {
				if (is_reversed_pam)
					tmp = sline[0] + string(m_dnabulgesize, 'N');
				else
					tmp = string(m_dnabulgesize, 'N') + sline[0];
				for (i = 1; i <= m_dnabulgesize; i++) {
					preNcnt = m_dnabulgesize - i;
					for (j = sline[0].find_first_not_of('N'); j < sline[0].find_last_not_of('N') + 2; j++) {
						if (is_reversed_pam) {
							compare = sline[0].substr(0, j) + string(i, 'N') + sline[0].substr(j) + string(preNcnt, 'N');
						} else {
							compare = string(preNcnt, 'N') + sline[0].substr(0, j) + string(i, 'N') + sline[0].substr(j);
						}
						if (compare == tmp)
							bi = make_pair("0", encode_bulge_index(j, is_reversed_pam));
						else
							bi = make_pair(to_string(i), encode_bulge_index(j, is_reversed_pam));
						add_compare(compare, ci, bi);
					}
				}
			}
			for (i = 1; i <= m_rnabulgesize; i++) {
				preNcnt = m_dnabulgesize + i;
				for (j = sline[0].find_first_not_of('N'); j < sline[0].find_last_not_of('N') + 1 - i; j++) {
					bi = make_pair(sline[0].substr(j, i), encode_bulge_index(j, is_reversed_pam));
					if (is_reversed_pam) {
						compare = sline[0].substr(0, j) + sline[0].substr(j + i) + string(preNcnt, 'N');
					} else {
						compare = string(preNcnt, 'N') + sline[0].substr(0, j) + sline[0].substr(j + i);
					}
					add_compare(compare, ci, bi);
					if (is_reversed_pam && has_right_pam) {
						const string guide_rc = reverse_complement_string(sline[0]);
						const unsigned int mirrored_index = (unsigned int)guide_rc.size() - j - i;
						const string compare_rc = guide_rc.substr(0, mirrored_index) + guide_rc.substr(mirrored_index + i) + string(preNcnt, 'N');
						add_compare(reverse_complement_string(compare_rc), ci, bi);
					}
				}
			}
			if (entrycnt == 0) {
				entrycnt = sline.size();
			} else if (entrycnt != sline.size()) {
				throw runtime_error("The number of entries below 2nd line should be consistent.");
			}
			linecnt++;
		}
	} catch (const exception& e) {
		cerr << "Critical error! " << e.what() << endl;
		exit(1);
	}
}

void Cas_OFFinder::readInputFile(const char* inputfile) {
	unsigned int dev_index;
	cl_uint zero = 0;

	if (strlen(inputfile) == 1 && inputfile[0] == '-') {
		parseInput(cin);
	} else {
		ifstream fi(inputfile, ios::in);
		parseInput(fi);
		fi.close();
	}

	m_patternlen = (cl_uint)(m_pattern.size());

	cl_char *cl_pattern = new cl_char[m_patternlen * 2];
	memcpy(cl_pattern, m_pattern.c_str(), m_patternlen);
	memcpy(cl_pattern + m_patternlen, m_pattern.c_str(), m_patternlen);
	set_complementary_sequence(cl_pattern+m_patternlen, m_patternlen);
	cl_int *cl_pattern_flags = new cl_int[m_patternlen * 2];
	set_seq_flags(cl_pattern_flags, cl_pattern, m_patternlen);
	set_seq_flags(cl_pattern_flags + m_patternlen, cl_pattern + m_patternlen, m_patternlen);

	for (dev_index = 0; dev_index < m_devnum; dev_index++) {
		m_patternbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char) * m_patternlen * 2, 0));
		m_patternflagbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_int) * m_patternlen * 2, 0));
		oclEnqueueWriteBuffer(m_queues[dev_index], m_patternbufs[dev_index], CL_FALSE, 0, sizeof(cl_char) * m_patternlen * 2, cl_pattern, 0, 0, 0);
		oclEnqueueWriteBuffer(m_queues[dev_index], m_patternflagbufs[dev_index], CL_FALSE, 0, sizeof(cl_int) * m_patternlen * 2, cl_pattern_flags, 0, 0, 0);
		
		m_comparebufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_char) * m_patternlen * 2, 0));
		m_compareflagbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_ONLY, sizeof(cl_uint) * m_patternlen * 2, 0));

		m_entrycountbufs.push_back(oclCreateBuffer(m_contexts[dev_index], CL_MEM_READ_WRITE, sizeof(cl_uint), 0));
		oclEnqueueWriteBuffer(m_queues[dev_index], m_entrycountbufs[dev_index], CL_FALSE, 0, sizeof(cl_uint), &zero, 0, 0, 0);
		oclFinish(m_queues[dev_index]);

		oclSetKernelArg(m_finderkernels[dev_index], 1, sizeof(cl_mem), &m_patternbufs[dev_index]);
		oclSetKernelArg(m_finderkernels[dev_index], 2, sizeof(cl_mem), &m_patternflagbufs[dev_index]);
		oclSetKernelArg(m_finderkernels[dev_index], 3, sizeof(cl_uint), &m_patternlen);
		oclSetKernelArg(m_finderkernels[dev_index], 5, sizeof(cl_mem), &m_entrycountbufs[dev_index]);

		oclSetKernelArg(m_comparerkernels[dev_index], 3, sizeof(cl_mem), &m_comparebufs[dev_index]);
		oclSetKernelArg(m_comparerkernels[dev_index], 4, sizeof(cl_mem), &m_compareflagbufs[dev_index]);
		oclSetKernelArg(m_comparerkernels[dev_index], 5, sizeof(cl_uint), &m_patternlen);
		oclSetKernelArg(m_comparerkernels[dev_index], 10, sizeof(cl_mem), &m_entrycountbufs[dev_index]);

		if (m_devtype != CL_DEVICE_TYPE_CPU) {
			oclSetKernelArg(m_finderkernels[dev_index], 7, sizeof(cl_char) * m_patternlen * 2, 0);
			oclSetKernelArg(m_finderkernels[dev_index], 8, sizeof(cl_int) * m_patternlen * 2, 0);
			oclSetKernelArg(m_comparerkernels[dev_index], 11, sizeof(cl_char) * m_patternlen * 2, 0);
			oclSetKernelArg(m_comparerkernels[dev_index], 12, sizeof(cl_int) * m_patternlen * 2, 0);
		}
	}

	delete[] cl_pattern;
	delete[] cl_pattern_flags;
}
