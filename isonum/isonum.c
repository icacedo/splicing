#include "toolbox.h"
#include "sequence.h"

/*
static go(int offset, int k, const ik_vec people, ik_ivec combos) {
	if (k == 0) return;
	for (int i = offset; i <= people->size -k ; i++) {
		ik_vec_push(combos, people->elem[i]);
		go(i+1, k-1, people, combos)
		int trash = ik_vec_pop(combos);
	}
}
*/


static void all_possible(const char *seq, int min_intron, int min_exon) {
	int len = strlen(seq);
	int nsites;
	ik_ivec dons = ik_ivec_new();
	ik_ivec accs = ik_ivec_new();
	
	for (int i = min_exon; i < len - min_exon; i++) {
		if      (seq[i  ] == 'G' && seq[i+1] == 'T') ik_ivec_push(dons, i);
		else if (seq[i-1] == 'A' && seq[i  ] == 'G') ik_ivec_push(accs, i);
	}
	
	nsites = dons->size < accs->size ? dons->size : accs->size;
	
	


	printf("%d sites %d %d\n", nsites, dons->size, accs->size);
	
}

static char *usage = "\
isonum - generate all possible isoforms of sequences\n\n\
usage: isonum <fasta file> [options]\n\
options:\n\
  -intron <int>  minimum intron size [35]\n\
  -exon   <int>  minimum exon size   [15]\n\
  -full          generate full report\n\
";

int main(int argc, char **argv) {
	char *file = NULL;  // path to fasta file
	int   intron = 35;  // minimum intron size
	int   exon   = 15;  // minimum exon size
	int   full   = 0;   // full report?
	ik_pipe io;
	ik_fasta ff;
	
	// Command Line Interface
	ik_set_program_name(argv[0]);
	ik_register_option("-intron", 1);
	ik_register_option("-exon", 1);
	ik_register_option("-full", 0);
	ik_parse_options(&argc, argv);
	
	if (argc == 1) ik_exit(1, "%s", usage);
	
	file = argv[1];
	if (ik_option("-intron")) intron = atoi(ik_option("-intron"));
	if (ik_option("-exon"))   exon   = atoi(ik_option("-exon"));
	if (ik_option("-full"))   full = 1;
	
	// main loop
	io = ik_pipe_open(file, "r");
	while ((ff = ik_fasta_read(io->stream)) != NULL) {
		all_possible(ff->seq, intron, exon);
		ik_fasta_free(ff);
	}
	
	// testing new functions
	ik_ivec iv = ik_ivec_new();
	for (int i = 1; i <= 5; i++) ik_ivec_push(iv, i);
	
	for (int i = 0; i < iv->size; i++) printf("%d ", iv->elem[i]);
	printf("\n");

	return 0;
}
