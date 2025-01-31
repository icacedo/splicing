# problem: geniso2 on the command line returns different probabilites
# when compared to locus.isoforms?

import isoform
from isoform import Locus
import subprocess
import json
import sys

with open('config_optiso2.json', 'r') as file:
    data = json.load(file)

minin = data['params']['minin']
minex = data['params']['minex']
flank = data['params']['flank']
limit = data['params']['limit']

dpwm = isoform.read_pwm(data['models']['dpwm'])
apwm = isoform.read_pwm(data['models']['apwm'])
elen = isoform.read_len(data['models']['elen'])
ilen = isoform.read_len(data['models']['ilen'])
emm = isoform.read_markov(data['models']['emm'])
imm = isoform.read_markov(data['models']['imm'])

fasta = '../../datacore2024/project_splicing/smallgenes/ch.2_1.fa'
gff3 = '../../datacore2024/project_splicing/smallgenes/ch.2_1.gff3'
name, seq = next(isoform.read_fasta(fasta))
models = (dpwm, apwm, emm, imm, elen, ilen)
weights = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

# the problem was you needed to transform icost before passing it to Locus
icost = isoform.prob2score(0.0)
locus = Locus(name, seq, minin, minex, flank, models, weights, 
                  icost, limit=2)
print('###############################')
print(locus.gff(sys.stdout))
print('###############################')
# are introns offset by 1 base? internal isoform index starts at 0

weight = []
total = 0
isos = locus.isoforms
print(len(isos))
for i in isos:
    print(i, '$$$')
print('###########')
print(isos, '1')
for tx in isos:
    s = 0 
    s += isoform.score_apwm(apwm, tx) * 1.0
    s += isoform.score_apwm(apwm, tx) * 1.0
    s += isoform.score_apwm(apwm, tx) * 1.0
    s += isoform.score_apwm(apwm, tx) * 1.0
    s += isoform.score_apwm(apwm, tx) * 1.0
    s += isoform.score_apwm(apwm, tx) * 1.0
    s += len(tx['introns']) * isoform.prob2score(0.0)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!
    # big problem: isos is not getting updated after exiting the loop
    tx['score'] = 100
    #print(tx, '&&&&')
print(isos, '2')
'''
    for tx in isos:
        #print(tx)
        for intron in tx['introns']:
            print(intron, tx['prob'], '@wow')
    weight = []
    total = 0
    for tx in isos:
        w = 2 ** s
        weight.append(w)
        total += w
    prob = [w / total for w in weight]       
    for p, tx in zip(prob, isos):
        tx['prob'] = p
    for tx in isos:
        for intron in tx['introns']:
            print(intron, tx['prob'], '@@@')
'''

 
subprocess.run([
    'geniso2', '../../datacore2024/project_splicing/smallgenes/ch.2_1.fa',
    '--dpwm', '../../isoforms/models/don.pwm',
    '--apwm', '../../isoforms/models/acc.pwm',
    '--emm', '../../isoforms/models/exon.mm', 
    '--imm', '../../isoforms/models/intron.mm', 
    '--elen', '../../isoforms/models/exon.len',
    '--ilen', '../../isoforms/models/intron.len',
    '--wdpwm', '1.0',
    '--wapwm', '1.0',
    '--wemm', '1.0',
    '--wimm', '1.0',
    '--welen', '1.0',
    '--wilen', '1.0',
    '--icost', '0.0',
    '--limit', '2'
], stdout=open('test.ch.2_1.geniso2.gff', 'w'))

