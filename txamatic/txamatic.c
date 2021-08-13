#include "toolbox.h"
#include "sequence.h"
#include "model.h"
#include "feature.h"

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

	ik_ivec_free(tmp);

	return indexes;
}

static int short_intron(const ik_ivec dons, const ik_ivec accs, int min) {
	for (int i = 0; i < dons->size; i++) {
		int len = accs->elem[i] - dons->elem[i] + 1;
		if (len < min) return 1;
	}
	return 0;
}

static int short_exon(const ik_ivec dons, const ik_ivec accs,
		int seqlen, int flank, int min) {
	
	// first exon
	int exon_beg = flank + 1;
	int exon_end = dons->elem[0] -1;
	int exon_len = exon_end - exon_beg + 1;
	if (exon_len < min) return 1;
	
	// last exon
	exon_beg = accs->elem[accs->size -1] + 1;
	exon_end = seqlen - flank -1;
	exon_len = exon_end - exon_beg + 1;
	if (exon_len < min) return 1;
	
	// internal exons
	for (int i = 1; i < dons->size; i++) {
		int beg = accs->elem[i-1] + 1;
		int end = dons->elem[i] -1;
		int len = end - beg;
		if (len < min) return 1;
	}
	return 0;
}

static double score_apwm(const ik_pwm pwm, const ik_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		ik_feat f = tx->introns->elem[i];
		double s = ik_pwm_score(pwm, f->seq, f->end -pwm->size +1);
		score += s;
	}
	return score;
}

static double score_dpwm(const ik_pwm pwm, const ik_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		ik_feat f = tx->introns->elem[i];
		double s = ik_pwm_score(pwm, f->seq, f->beg);
		score += s;
	}
	return score;
}

static double score_elen(const ik_len model, const ik_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->exons->size; i++) {
		ik_feat f = tx->exons->elem[i];
		int len = f->end - f->beg + 1;
		double s = ik_len_score(model, len);
		score += s;
	}
	return score;
}

static double score_ilen(const ik_len model, const ik_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		ik_feat f = tx->introns->elem[i];
		int len = f->end - f->beg + 1;
		double s = ik_len_score(model, len);
		score += s;
	}
	return score;
}

static double score_emm(const ik_mm mm, const ik_mRNA tx) {
	double score = 0;
	for (int i = 0; i < tx->exons->size; i++) {
		ik_feat f = tx->exons->elem[i];
		double s = ik_mm_score(mm, f->seq, f->beg, f->end);
		score += s;
	}
	return score;
}

static double score_imm(const ik_mm mm, const ik_mRNA tx,
		const ik_pwm dpwm, const ik_pwm apwm) {
	double score = 0;
	for (int i = 0; i < tx->introns->size; i++) {
		ik_feat f = tx->introns->elem[i];
		double s = ik_mm_score(mm, f->seq, f->beg + dpwm->size,
			f->end - apwm->size);
		score += s;
	}
	return score;
}

static int cmptx(const ik_mRNA a, const ik_mRNA b) {
	if      (a->score < b->score) return -1;
	else if (a->score > b->score) return  1;
	else                          return  0;
}

static int mysort(const void *a, const void *b) {
	return cmptx( *(ik_mRNA *)b, *(ik_mRNA *)a );
}

struct isoret {
	int     dons;
	int     accs;
	int     trials;
	int     ishort;
	int     eshort;
	ik_vec  mRNAs;
};

static struct isoret isoforms(const char *seq,
		int min_intron, int min_exon,
		int max, int flank,
		ik_pwm apwm, ik_pwm dpwm,
		ik_mm emm, ik_mm imm,
		ik_len elen, ik_len ilen) {
	int len = strlen(seq);
	int nsites;
	ik_ivec dons = ik_ivec_new();
	ik_ivec accs = ik_ivec_new();
	struct isoret ret;

	/* init */
	for (int i = flank + min_exon; i < len - flank - min_exon; i++) {
		if (seq[i] == 'G' && seq[i+1] == 'T') ik_ivec_push(dons, i);
		if (seq[i] == 'A' && seq[i+1] == 'G') ik_ivec_push(accs, i+1);
	}

	nsites = dons->size < accs->size ? dons->size : accs->size;
	if (nsites > max) nsites = max;

	int trials = 0;
	int ishort = 0;
	int eshort = 0;

	// main loop
	ik_vec txs = ik_vec_new();
	for (int n = 1; n <= nsites; n++) {
		ik_vec dcombos = get_combinations(dons, n);
		ik_vec acombos = get_combinations(accs, n);

		for (int i = 0; i < dcombos->size; i++) {
			for (int j = 0; j < acombos->size; j++) {
				ik_ivec dsites = dcombos->elem[i];
				ik_ivec asites = acombos->elem[j];
				assert(dsites->size == asites->size);

				trials += 1;
				if (short_intron(dsites, asites, min_intron)) {
					ishort++;
					continue;
				}
				if (short_exon(dsites, asites, len, flank, min_exon)) {
					eshort++; // what about short flanks?
					continue;
				}

				ik_mRNA tx = ik_mRNA_new(seq, flank, len -flank,
						dsites, asites);
				double score = 0;
				if (apwm) score += score_apwm(apwm, tx);
				if (dpwm) score += score_dpwm(dpwm, tx);
				if (elen) score += score_elen(elen, tx);
				if (ilen) score += score_ilen(ilen, tx);
				if (emm)  score += score_emm(emm, tx);
				if (imm)  score += score_imm(emm, tx, dpwm, apwm);
				tx->score = score;
				ik_vec_push(txs, tx);
			}
		}
		qsort(txs->elem, txs->size, sizeof(ik_mRNA), mysort);

		ret.dons   = dons->size;
		ret.accs   = accs->size;
		ret.trials = trials;
		ret.ishort = ishort;
		ret.eshort = eshort;
		ret.mRNAs  = txs;

		// free combos
		for (int i = 0; i < dcombos->size; i++) {
			ik_ivec v = dcombos->elem[i];
			ik_ivec_free(v);
		}
		ik_vec_free(dcombos);
		for (int i = 0; i < acombos->size; i++) {
			ik_ivec v = acombos->elem[i];
			ik_ivec_free(v);
		}
		ik_vec_free(acombos);
	}

	// clean up
	ik_ivec_free(dons);
	ik_ivec_free(accs);

	return ret;
}

static char *usage = "\
txamatic - generate all possible isoforms from sequences\n\n\
usage: txamatic <fasta file> [options]\n\
options:\n\
  -intron <int>  minimum intron size [35]\n\
  -exon   <int>  minimum exon size   [25]\n\
  -max    <int>  maximum # introns   [3]\n\
  -flank  <int>  flank to ignore     [99]\n\
  -apwm   <file> use acceptor pwm\n\
  -dpwm   <file> use donor pwm\n\
  -emm    <file> use exon Markov model\n\
  -imm    <file> use intron Markov model (requires -dpwm & -apwm)\n\
  -elen   <file> use exon length model\n\
  -ilen   <file> use intron length model\n\
  -full          full report\n\
";

int main(int argc, char **argv) {
	char *file   = NULL; // path to fasta file
	int   intron = 35;   // minimum intron size
	int   exon   = 25;   // minimum exon size
	int   max    = 3;    // maximum number of introns
	int   flank  = 99;  // promoter and such

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
	ik_register_option("-full", 0);
	ik_parse_options(&argc, argv);
	if (argc == 1) ik_exit("%s", usage);

	// CLI - harvest
	file = argv[1];
	if (ik_option("-intron")) intron = atoi(ik_option("-intron"));
	if (ik_option("-exon"))   exon   = atoi(ik_option("-exon"));
	if (ik_option("-max"))    max    = atoi(ik_option("-max"));
	if (ik_option("-flank"))  flank  = atoi(ik_option("-flank"));
	if (ik_option("-apwm"))   apwm   = ik_pwm_read(ik_option("-apwm"));
	if (ik_option("-dpwm"))   dpwm   = ik_pwm_read(ik_option("-dpwm"));
	if (ik_option("-emm"))    emm    = ik_mm_read(ik_option("-emm"));
	if (ik_option("-imm"))    imm    = ik_mm_read(ik_option("-imm"));
	if (ik_option("-elen"))   elen   = ik_len_read(ik_option("-elen"));
	if (ik_option("-ilen"))   ilen   = ik_len_read(ik_option("-ilen"));

	// main loop
	io = ik_pipe_open(file, "r");
	while ((ff = ik_fasta_read(io->stream)) != NULL) {
		struct isoret iso = isoforms(ff->seq, intron, exon, max, flank,
			apwm, dpwm, emm, imm, elen, ilen);

		// output
		ik_vec txs = iso.mRNAs;
		printf("seq: %s\n", ff->def);
		printf("len: %lu\n", strlen(ff->seq));
		printf("dons: %d\n", iso.dons);
		printf("accs: %d\n", iso.accs);
		printf("trials: %d\n", iso.trials);
		printf("ishort: %d\n", iso.ishort);
		printf("eshort: %d\n", iso.eshort);
		printf("mRNAs: %d\n", iso.mRNAs->size);
		if (ik_option("-full")) {
			for (int i = 0; i < txs->size; i++) {
				ik_mRNA tx = txs->elem[i];
				printf("%f", tx->score);
				for (int j = 0; j < tx->exons->size; j++) {
					ik_feat f = tx->exons->elem[j];
					printf(" %d..%d", f->beg+1, f->end+1);
				}
				printf("\n");
			}
		}

		// clean up
		ik_fasta_free(ff);
		for (int i = 0; i < txs->size; i++) {
			ik_mRNA m = txs->elem[i];
			ik_mRNA_free(m);
		}
		ik_vec_free(txs);
	}

	return 0;
}
