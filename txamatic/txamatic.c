#include "toolbox.h"
#include "sequence.h"


void make_combi_util(ik_vec ans, ik_ivec tmp, int n, int left, int k) {
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
		make_combi_util(ans, tmp, n, i + 1, k - 1);
		ik_ivec_pop(tmp);
	}
}

ik_vec make_combinations(const ik_ivec sites, int k) {
	int idx, val;
	ik_ivec iv;
	ik_vec  indexes = ik_vec_new();
	ik_ivec tmp = ik_ivec_new();
	make_combi_util(indexes, tmp, sites->size, 1, k);
	
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

int short_intron(ik_ivec dons, ik_ivec accs, int min_intron) {
	for (int i = 0; i < dons->size; i++) {
		int len = accs->elem[i] - dons->elem[i];
		if (len < min_intron) return 1;
	}
	return 0;
}

int short_exon(ik_ivec dons, ik_ivec accs, int min_exon) {
	for (int i = 1; i < dons->size; i++) {
		int beg = accs->elem[i-1] + 1;
		int end = dons->elem[i] -1;
		int len = end - beg;
		if (len < min_exon) return 1;
	}
	return 0;
}

static void all_possible(const char *seq, int min_intron, int min_exon, int max, int flank) {
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
	//printf("min %d, don %d, acc %d\n", nsites, dons->size, accs->size);
	
	int trials = 0;
	int ishort = 0;
	int eshort = 0;
	int passed = 0;
	for (int k = 1; k <= nsites; k++) {
		ik_vec dcs = make_combinations(dons, k);
		ik_vec acs = make_combinations(accs, k);
		
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
				
				for (int a = 0; a < k; a++) {
					printf("%d..%d ", dv->elem[a], av->elem[a]);
				}
				printf("\n");
			}
		}
	}
	
	printf("%d %d %d %d\n", trials, ishort, eshort, passed);
	
	// THERE IS NO GARBAGE COLLECTION! //
}

static char *usage = "\
txamatic - generate all possible transcripts from sequences\n\n\
usage: isonum <fasta file> [options]\n\
options:\n\
  -intron <int>  minimum intron size [35]\n\
  -exon   <int>  minimum exon size   [25]\n\
  -max    <int>  maximum # introns   [3]\n\
  -flank  <int>  flank to ignore     [100]\n\
  -full          generate full report\n\
";

int main(int argc, char **argv) {
	char *file = NULL;  // path to fasta file
	int   intron = 35;  // minimum intron size
	int   exon   = 25;  // minimum exon size
	int   max    = 3;   // maximum number of introns
	int   flank  = 100; // promoter and such
	int   full   = 0;   // full report?
	ik_pipe io;
	ik_fasta ff;
	
	// Command Line Interface
	ik_set_program_name(argv[0]);
	ik_register_option("-intron", 1);
	ik_register_option("-exon", 1);
	ik_register_option("-max", 1);
	ik_register_option("-flank", 1);
	ik_register_option("-full", 0);
	ik_parse_options(&argc, argv);
	
	if (argc == 1) ik_exit(1, "%s", usage);
	
	file = argv[1];
	if (ik_option("-intron")) intron = atoi(ik_option("-intron"));
	if (ik_option("-exon"))   exon   = atoi(ik_option("-exon"));
	if (ik_option("-max"))    max    = atoi(ik_option("-max"));
	if (ik_option("-flank"))  flank  = atoi(ik_option("-flank"));
	if (ik_option("-full"))   full = 1;
	
	// main loop
	io = ik_pipe_open(file, "r");
	while ((ff = ik_fasta_read(io->stream)) != NULL) {
		all_possible(ff->seq, intron, exon, max, flank);
		ik_fasta_free(ff);
	}

	return 0;
}
