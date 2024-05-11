import isomod as im
import test_isoform as isoform
import sys

seq = (
    'TTGTTTTTATTTTGCTGAATTCTTTTGTTTTGGTTGTATACCAATTTTTTAAACTTATCCTTGTTTTCATTTAGTTCAAA'
    'TCTTACAATTTTTTTTCAAGAACCATGTTCAAAAAATTCGACGAGAAGGAAGATGTGACAGGCGCTACTCAGCTCAAGTC'
)

d, a = im.get_daseq((26, 131), seq)
print(d)
print(a)

seeq = 'CCCCGTCCCCCAGCCCCC'
d, a = im.get_daseq((4, 12), seeq)
print(d)
print(a)

pwm = sys.argv[1]



print(isoform.read_pwm(pwm))