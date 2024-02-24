import isosort_lib as isl
import argparse
import modelib as ml
# use ln -s for modelib in gff_analysis/

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
	help='fasta file')
parser.add_argument('wb_gff', type=str, metavar='<file>',
	help='gff file')

parser.add_argument('--exon_len', type=str, metavar='<file>',
	help='exon length model .tsv')
parser.add_argument('--intron_len', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--exon_mm', type=str, metavar='<file>',
	help='exon markov model .tsv')
parser.add_argument('--intron_mm', type=str, metavar='<file>',
	help='intron markov model .tsv')
parser.add_argument('--donor_pwm', type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--acceptor_pwm', type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--icost', required=False, type=float, default=0.0, 
	metavar='<float>', help='intron cost %(default).2d')

args = parser.parse_args()

seq = None
for seqid, seq in ml.read_fastas(args.fasta):
    seqid = seqid
    seq = seq

wbginfo = isl.get_wbgene_info(args.wb_gff, seq)


print(wbginfo)

re_elen_pdf, re_elen_log2 = ml.read_exin_len(args.exon_len)
ea, eb, eg = ml.read_len_params(args.exon_len)

re_ilen_pdf, re_ilen_log2 = ml.read_exin_len(args.intron_len)
ia, ib, ig = ml.read_len_params(args.intron_len)

re_emm_prob, re_emm_log2 = ml.read_exin_mm(args.exon_mm)
re_imm_prob, re_imm_log2 = ml.read_exin_mm(args.intron_mm)

re_dppm, re_dpwm = ml.read_pwm(args.donor_pwm)
re_appm, re_apwm = ml.read_pwm(args.acceptor_pwm)

# THE REASON SCORES DONT MATCH #
# apc_isogen.py and modelib.py uses 0 indexing
# wbginfo has 1 indexing

for gene in wbginfo:
    total_score = 0
    for exon in wbginfo[gene]['exons']:
        print(exon)
        exon = (exon[0]-1, exon[1]-1)
        elen_score = ml.get_exin_len_score(exon, re_elen_log2, ea, eb, eg)
        print(elen_score)
        emm_score = ml.get_exin_mm_score(exon, seq, re_emm_log2)
        print(emm_score)
        escore = elen_score + emm_score
        print(escore)        
    print('#####')
    for intron in wbginfo[gene]['introns']:
        print(intron)
        intron = (intron[0]-1, intron[1]-1)
        ilen_score = ml.get_exin_len_score(intron, re_ilen_log2, ia, ib, ig)
        print(ilen_score)
        imm_score = ml.get_exin_mm_score(intron, seq, re_imm_log2, 'GT', 'AG')
        print(imm_score)
        dseq, aseq = ml.get_donacc_seq(intron, seq)
        dpwm_score = ml.get_donacc_pwm_score(dseq, re_dpwm)
        apwm_score = ml.get_donacc_pwm_score(aseq, re_apwm)
        print(dpwm_score)
        print(apwm_score)



# need to include icost in wb gene score

