#include "toolbox.h"
#include "sequence.h"
#include "model.h"

static void combo(ik_vec ans, ik_ivec tmp, int n, int left, int k) {
	if (k == 0) {
		ik_ivec keep = ik_ivec_new();
		for (int i = 0; i < tmp->size; i++) ik_ivec_push(keep, tmp->elem[i]);
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
	int     idx, val;
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

static double score_apwm(const char *seq, const ik_pwm pwm,
		const ik_ivec accs) {
	double score = 0;
	for (int i = 0; i < accs->size; i++) {
		int pos = accs->elem[i] - pwm->size +1;
		score += ik_score_pwm(pwm, seq, pos);
	}
	return score - (-2 * pwm->size);
}

static double score_dpwm(const char *seq, const ik_pwm pwm,
		const ik_ivec dons) {
	double score = 0;
	for (int i = 0; i < dons->size; i++) {
		int pos = dons->elem[i];
		score += ik_score_pwm(pwm, seq, pos);
	}
	return score - (-2 * pwm->size);
}

static double score_elen(const ik_len len,
		const ik_ivec dons, const ik_ivec accs) {
	return 0;
}

static double score_ilen(const ik_len len,
		const ik_ivec dons, const ik_ivec accs) {
	return 0;
}

static double score_emm(const char *seq, const ik_mm mm,
		const ik_ivec dons, const ik_ivec accs) {
	return 0;
}

static double score_imm(const char *seq, const ik_mm mm,
		const ik_ivec dons, const ik_ivec accs) {
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
		if (seq[i] == 'G' && seq[i+1] == 'T') ik_ivec_push(dons, i);
		if (seq[i] == 'A' && seq[i+1] == 'G') ik_ivec_push(accs, i+1);
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
				double score = 0;
				if (apwm) score += score_apwm(seq, apwm, av);
				if (dpwm) score += score_dpwm(seq, dpwm, dv);
				if (elen) score += score_elen(elen, dv, av);
				if (ilen) score += score_ilen(ilen, dv, av);
				if (emm)  score += score_emm(seq, emm, dv, av);
				if (imm)  score += score_imm(seq, imm, dv, av);
				
				printf("%g", score);
				for (int a = 0; a < k; a++) {
					printf(" %d..%d", dv->elem[a], av->elem[a]);
				}
				printf("\n");
				
				//ik_exit("testing");
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
	if (argc == 1) ik_exit("%s", usage);
	
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
