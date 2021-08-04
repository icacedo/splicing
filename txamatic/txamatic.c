#include "toolbox.h"
#include "sequence.h"
#include <assert.h>

// some of this may go into ik library

int kmer2dec(const char *kmer) {
	int k = strlen(kmer);
	int idx = 0;
	for (int i = 0; i < k; i++) {
		switch (kmer[i]) {
			case 'A': idx += pow(4, (k -i -1)) * 0; break;
			case 'C': idx += pow(4, (k -i -1)) * 1; break;
			case 'G': idx += pow(4, (k -i -1)) * 2; break;
			case 'T': idx += pow(4, (k -i -1)) * 3; break;
			default: ik_exit(1, "invalid nt character: %c\n", kmer[i]);
		}
	}
	return idx;
}

double prob2score(double p) {
	if (p == 0) return -100;
	return log(p) / log(2);
}

//

struct ik_PWM {
	char    *name;   // acceptor, donor
	int      size; // eg. 6
	double **score;  // score[pos][nt]
};
typedef struct ik_PWM * ik_pwm;

ik_pwm ik_read_pwm(const char *filename) {
	char    *line = NULL;
	size_t   len = 0;
	ssize_t  read;
	ik_pipe  io  = ik_pipe_open(filename, "r");
	char     blah[256];
	int      size;
	double **score = NULL;
	double   a, c, g, t;
	int      row = 0;
	
	while ((read = getline(&line, &len, io->stream)) != -1) {
		if (line[0] == '#') {
			assert(sscanf(line, "# PWM %s %d", blah, &size) == 2);
			score = malloc(sizeof(double*) * size);
			for (int i = 0; i < size; i++) {
				score[i] = malloc(sizeof(double) * 4);
			}
		} else if (sscanf(line, "%lf %lf %lf %lf", &a, &c, &g, &t) == 4) {
			score[row][0] = prob2score(a);
			score[row][1] = prob2score(c);
			score[row][2] = prob2score(g);
			score[row][3] = prob2score(t);
			row++;
		}
	}
	ik_pipe_close(io);
	if (line) free(line);
	
	ik_pwm model = malloc(sizeof(struct ik_PWM));
	model->name = malloc(strlen(filename)+1);
	strcpy(model->name, filename);
	model->size = size;
	model->score = score;
	
	return model;
}

double ik_score_pwm(const ik_pwm pwm, const char *seq, int pos) {
	double p = 0;
	for (int i = 0; i < pwm->size; i++) {
		switch (seq[i+pos]) {
			case 'A': p += pwm->score[i][0]; break;
			case 'C': p += pwm->score[i][1]; break;
			case 'G': p += pwm->score[i][2]; break;
			case 'T': p += pwm->score[i][3]; break;
			default: ik_exit(1, "invalid nt: %c\n", seq[i+pos]);
		}
	}
	return p;
}

struct ik_MM {
	char   *name;   // exon, intron
	int     k;      // kmer size
	int     size;   // size of array
	double *score;  // score[base-4 kmer] = value
};
typedef struct ik_MM * ik_mm;

ik_mm ik_read_mm(const char *filename) {
	char    *line = NULL;
	size_t   len = 0;
	ssize_t  read;
	ik_pipe  io  = ik_pipe_open(filename, "r");
	double  *score = NULL;
	char     kmer[16];
	char     blah[256];
	int      size;
	double   p;
	
	while ((read = getline(&line, &len, io->stream)) != -1) {
		if (line[0] == '#') {
			assert(sscanf(line, "# MM %s %d", blah, &size) == 2);
			score = malloc(sizeof(double) * size);
		} else if (sscanf(line, "%s %lf", kmer, &p) == 2) {
			int idx = kmer2dec(kmer);
			score[idx] = prob2score(p);
		}
	}
	ik_pipe_close(io);
	if (line) free(line);
	
	ik_mm model = malloc(sizeof(struct ik_MM));
	model->name = malloc(strlen(filename)+1);
	strcpy(model->name, filename);
	model->k = strlen(kmer);
	model->size = size;
	model->score = score;
	
	return model;
}

double ik_score_mm(const ik_mm mm, const char *seq, int pos, int end) {
	char kmer[16];
	double p = 0;
	
	if (pos < mm->k) pos = mm->k;
	for (int i = pos; i <= end; i++) {
		strncpy(kmer, seq+i, mm->k);
		kmer[mm->k] = '\0';
		int idx = kmer2dec(kmer);
		double val = mm->score[idx];
		p += val;
	}
	
	return p;
}

struct ik_LEN {
	char   *name;   // exon, intron
	ik_fvec score;  // the values for defined region, ik_fvec
	double  tail;
};
typedef struct ik_LEN * ik_len;

static double find_tail(double final, int x) {
	double lo = 0;
	double hi = 1000;
	double m;
	
	while (hi - lo > 1) {
		m = (hi + lo) / 2;
		double p = 1 / m;
		double f = pow(1-p, x-1) * p;
		//printf("%f %f %f %f\n", lo, hi, m, f);
		if (f < final) lo += (m - lo) / 2;
		else           hi -= (hi - m) / 2;
	}

	return m;
}

ik_len ik_read_len(const char *filename) {
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	ik_pipe io = ik_pipe_open(filename, "r");
	ik_fvec vals = ik_fvec_new();
	double p;
	
	while ((read = getline(&line, &len, io->stream)) != -1) {
		if (line[0] == '#') continue;
		sscanf(line, "%lf", &p);
		ik_fvec_push(vals, prob2score(p));
	}
	ik_pipe_close(io);
	if (line) free(line);
	
	ik_len model = malloc(sizeof(struct ik_LEN));
	model->name = malloc(strlen(filename)+1);
	strcpy(model->name, filename);
	model->score = vals;
	model->tail = find_tail(vals->elem[vals->size-1], vals->size);
	
	return model;
}

double ik_score_len(const ik_len len, int x) {
	assert(x > 0);
	if (x >= len->score->size) {
		double f;
		double p = 1 / len->tail;
		f = pow(1-p, x-1) * p;
		return f;
	} else {
		return len->score->elem[x];
	}
}

static void combo(ik_vec ans, ik_ivec tmp, int n, int left, int k) {
	if (k == 0) {
		ik_ivec keep = ik_ivec_new();
		for (int i = 0; i < tmp->size; i++) {
			ik_ivec_push(keep, tmp->elem[i]);
		}
		ik_vec_push(ans, (void*)keep);
		return;
	}
	
	for (int i = left; i <= n; i++) {
		ik_ivec_push(tmp, i - 1);
		combo(ans, tmp, n, i + 1, k - 1);
		ik_ivec_pop(tmp);
	}
}

static ik_vec get_combinations(const ik_ivec sites, int k) {
	int idx, val;
	ik_ivec iv;
	ik_vec  indexes = ik_vec_new();
	ik_ivec tmp = ik_ivec_new();
	combo(indexes, tmp, sites->size, 1, k);
	
	for (int i = 0; i < indexes->size; i++) {
		iv = indexes->elem[i];
		for (int j = 0; j < iv->size; j++) {
			idx = iv->elem[j];
			val = sites->elem[idx];
			iv->elem[j] = val;
		}
	}
	
	return indexes;
}

static int short_intron(ik_ivec dons, ik_ivec accs, int min_intron) {
	for (int i = 0; i < dons->size; i++) {
		int len = accs->elem[i] - dons->elem[i];
		if (len < min_intron) return 1;
	}
	return 0;
}

static int short_exon(ik_ivec dons, ik_ivec accs, int min_exon) {
	for (int i = 1; i < dons->size; i++) {
		int beg = accs->elem[i-1] + 1;
		int end = dons->elem[i] -1;
		int len = end - beg;
		if (len < min_exon) return 1;
	}
	return 0;
}

static void all_possible(const char *seq,
		int min_intron, int min_exon,
		int max, int flank,
		ik_pwm apwm, ik_pwm dpwm,
		ik_mm emm, ik_mm imm,
		ik_len elen, ik_len ilen) {
	int len = strlen(seq);
	int nsites;
	ik_ivec dons = ik_ivec_new();
	ik_ivec accs = ik_ivec_new();
	
	for (int i = flank; i < len - flank; i++) {
		if (seq[i]   == 'G' && seq[i+1] == 'T') ik_ivec_push(dons, i);
		if (seq[i-1] == 'A' && seq[i]   == 'G') ik_ivec_push(accs, i);
	}
	
	nsites = dons->size < accs->size ? dons->size : accs->size;
	if (nsites > max) nsites = max;
	
	int trials = 0;
	int ishort = 0;
	int eshort = 0;
	int passed = 0;
	for (int k = 1; k <= nsites; k++) {
		ik_vec dcs = get_combinations(dons, k);
		ik_vec acs = get_combinations(accs, k);
		
		for (int i = 0; i < dcs->size; i++) {
			for (int j = 0; j < acs->size; j++) {
				ik_ivec dv = dcs->elem[i];
				ik_ivec av = acs->elem[j];
				
				trials += 1;
				if (short_intron(dv, av, min_intron)) {
					ishort++;
					continue;
				}
				if (short_exon(dv, av, min_exon)) {
					eshort++;
					continue;	
				}
				passed++;
				
				// output
				// need to adjust probability depending on models
				//double p = 1.0;
				//if (apwm) p *= ik_score_pwm(apwm, seq, pos);
				
				/*
				printf("%g ", p);
				for (int a = 0; a < k; a++) {
					printf("%d..%d ", dv->elem[a], av->elem[a]);
				}
				printf("\n");
				*/
			}
		}
		
		for (int i = 0; i < dcs->size; i++) ik_ivec_free(dcs->elem[i]);
		ik_vec_free(dcs);
		for (int i = 0; i < acs->size; i++) ik_ivec_free(acs->elem[i]);
		ik_vec_free(acs);
	}
	
	fprintf(stderr, "don:%d acc:%d n:%d in:%d ex:%d ok:%d\n",
		dons->size, accs->size,
		trials, ishort, eshort, passed);
	
	ik_ivec_free(dons);
	ik_ivec_free(accs);
	
}

static char *usage = "\
txamatic - generate all possible transcripts from sequences\n\n\
usage: isonum <fasta file> [options]\n\
options:\n\
  -intron <int>  minimum intron size [35]\n\
  -exon   <int>  minimum exon size   [25]\n\
  -max    <int>  maximum # introns   [3]\n\
  -flank  <int>  flank to ignore     [100]\n\
  -apwm   <file> use acceptor pwm\n\
  -dpwm   <file> use donor pwm\n\
  -emm    <file> use exon Markov model\n\
  -imm    <file> use intron Markov model\n\
  -elen   <file> use exon length model\n\
  -ilen   <file> use intron length model\n\
";

int main(int argc, char **argv) {
	char *file   = NULL; // path to fasta file
	int   intron = 35;   // minimum intron size
	int   exon   = 25;   // minimum exon size
	int   max    = 3;    // maximum number of introns
	int   flank  = 100;  // promoter and such
	
	ik_pwm apwm = NULL; // acceptor pwm
	ik_pwm dpwm = NULL; // donor pwm
	ik_mm  emm  = NULL; // exon Markov model
	ik_mm  imm  = NULL; // intron Markov model
	ik_len elen = NULL; // exon length model
	ik_len ilen = NULL; // intron length model
	
	ik_pipe  io = NULL; // for reading fasta files
	ik_fasta ff = NULL; // for individual fasta entries
	
	// CLI - setup
	ik_set_program_name(argv[0]);
	ik_register_option("-intron", 1);
	ik_register_option("-exon", 1);
	ik_register_option("-max", 1);
	ik_register_option("-flank", 1);
	ik_register_option("-apwm", 1);
	ik_register_option("-dpwm", 1);
	ik_register_option("-emm", 1);
	ik_register_option("-imm", 1);
	ik_register_option("-elen", 1);
	ik_register_option("-ilen", 1);
	ik_parse_options(&argc, argv);
	if (argc == 1) ik_exit(1, "%s", usage);
	
	// CLI - harvest 
	file = argv[1];
	if (ik_option("-intron")) intron = atoi(ik_option("-intron"));
	if (ik_option("-exon"))   exon   = atoi(ik_option("-exon"));
	if (ik_option("-max"))    max    = atoi(ik_option("-max"));
	if (ik_option("-flank"))  flank  = atoi(ik_option("-flank"));
	if (ik_option("-apwm"))   apwm   = ik_read_pwm(ik_option("-apwm"));
	if (ik_option("-dpwm"))   dpwm   = ik_read_pwm(ik_option("-dpwm"));
	if (ik_option("-emm"))    emm    = ik_read_mm(ik_option("-emm"));
	if (ik_option("-imm"))    imm    = ik_read_mm(ik_option("-imm"));
	if (ik_option("-elen"))   elen   = ik_read_len(ik_option("-elen"));
	if (ik_option("-ilen"))   ilen   = ik_read_len(ik_option("-ilen"));
	
	// main loop
	io = ik_pipe_open(file, "r");
	while ((ff = ik_fasta_read(io->stream)) != NULL) {
		all_possible(ff->seq, intron, exon, max, flank,
			apwm, dpwm, emm, imm, elen, ilen
		);
		ik_fasta_free(ff);
	}

	return 0;
}
