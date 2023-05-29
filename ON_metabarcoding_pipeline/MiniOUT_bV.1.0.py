
from itertools import combinations



import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import numpy as np
import scipy
from scipy import stats
from collections import Counter
import os

def makefig1(aa,b,a):
    fig, ax = plt.subplots(figsize=(8, 3))
    ax00 = plt.subplot2grid((1, 3), (0, 0), colspan=1, rowspan=1)
    ax0 = plt.subplot2grid((1, 3), (0, 1), colspan=1, rowspan=1)
    ax1 = plt.subplot2grid((1, 3), (0, 2), colspan=1, rowspan=1)
    ax00.hist(aa, color='grey', bins=100, )
    ax0.hist(b, color='grey', bins=100, )
    ax1.hist(a, color='grey', bins=100, )
    ax00.set_ylabel('Number of reads')
    ax00.set_xlabel('Reads lenght')
    ax0.set_xlabel('Reads lenght')
    ax1.set_xlabel('Reads lenght')

    ax00.set_title('All reads')
    ax0.set_title('Clipped reads')
    ax1.set_title(str(min) + '-' + str(max) + ' reads')
    plt.tight_layout()
    plt.savefig('Figure_1.png', dpi=800, format='png')


def makefig2(a):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax0 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)

    ax0.hist(a, color='grey', bins=100, )
    ax0.set_xticks([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])

    ax0.set_ylabel('Number of occurence')
    ax0.set_xlabel('Distance')
    ax0.set_xlim(0,1)

    plt.tight_layout()
    plt.savefig('Figure_2.png', dpi=800, format='png')



def dist(s1, s2):
    n1, n2 = s1.split('\n',1)[0], s2.split('\n',1)[0]
    s1, s2 = s1.split('\n',1)[1].replace('\n',''), s2.split('\n',1)[1].replace('\n','')
    tag, nm = 0, 0
    for x, y in zip(s1,s2):
        if x == '-' or y == '-':
            if tag == 0:
                nm = nm + 1
        else:
            tag = 1




def sortLenght(inp):
    #useless function
    n = 0
    bef, aft, seqL = [], [], []
    for x in inp[1:]:
        lg = len(x.split('\n', 1)[1])
        bef.append(lg)
        if lg >= min and lg <= max:
            aft.append(lg)
            seqL.append(x)
            n = n + 1
    fout = open('temp1_1.fas', 'w')
    ftp = open('temp1_2.fas', 'w')
    tpl = []
    for S in random.sample(seqL,100):  #random select 100 sequences
        fout.write('>' + S)
        tpl.append(S)
    fout.close()
    for S in seqL:
        if S not in tpl:
            ftp.write('>' + S)
    makefig1(bef,aft)
    return n, seqL


def mean_confidence_interval(data, confidence):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


def cleanT(dw):
    os.system('blastn -db ref.fas -query temp0.fas -max_target_seqs 1 -num_threads 24 -outfmt 6 > vectRes')
    print('blast done')
    fx = open('vectRes')
    bl = fx.readlines()
    fx.close()
    dicS = {}
    for x in bl:
        dicS[x.split('\t')[0]] = x.split('\t')[7]
    fx = open("temp0.fas")
    ct = fx.read().split('>')
    fx.close()
    n, t = 0, 0
    bef, aft,seqL,uncliped = [], [], [],[]
    fclip = open('cliped_seq.fas','w')
    for x in ct[1:]:
        try:
            clip = int(dicS[x.split(' ')[0]])
        except:
            #clip = 140      ########## if vector was removed during sequencing, this should be 0
            clip = 0
            pass
        seq = x.split('\n',1)[1]
        uncliped.append(len(seq))
        seq = seq[clip:]
        #print('>' + x.split(' ')[0] + '\n' + seq + '\n')
        fclip.write('>' + x.split(' ')[0] + '\n')  #fclip.write('>' + x.split(' ')[0] + '\n' + seq + '\n')
        fd = 0
        li_w = ['A' * 30, 'T' * 30, 'C' * 30, 'G' * 30, 'GA' * 20, 'GT' * 20, 'GC' * 20, 'TA' * 20, 'TC' * 20,
                    'TG' * 20, 'CA' * 20, 'CT' * 20, 'CG' * 20]
        for mot in li_w:
            if mot in seq:
                fd = 1
        if fd == 0:
            n = n + 1
            #aft.append(len(seq))
            if len(seq) > 199:     #249
                bef.append(len(seq))
            seqL.append((x.split(' ')[0], seq))

    #mn, dw, upp = mean_confidence_interval(bef,0.9)
    fout = open('temp1.fas', 'w')
    seqF = []
    for x in seqL:
        if len(x[1]) >= dw and len(x[1]) <= 1000:
            aft.append(len(x[1]))
            #fout.write('>' + x[0] + '\n' + x[1] + '\n')
            
            fout.write('>' + x[0] + '\n')
            seqF.append(x[0])# + '\n' + x[1] + '\n')

    fout.close()
    print(len(aft))
    makefig1(uncliped, bef, aft)
    return n, len(ct),seqF


def buildclusters(ssq):
    final = open('Clusters_1st_round.fasta','w')
    fx = open('temp2.opti_mcc.list')
    line = fx.readlines()[1].split('\t')
    fx.close()
    dico = {}
    for x in ssq:
        dico[x.split('\n')[0]] = x.split('\n')[1]
    n = 1
    for clust in line[2:]:
        fn = open('cluster_' + str(n) + '.fas','w')
        k = 0
        for y in clust.split(','):
            fn.write('>' + y.replace('\n','') + '\n' + dico[y.replace('\n','')] + '\n')
            k = k + 1
        fn.close()
        if k >= 3:
            os.system('mafft --auto cluster_' + str(n) + '.fas > cluster_aln_' + str(n) + '.fas')
            cs = consensus('cluster_aln_' + str(n) + '.fas')
            final.write('>cluster_aln_' + str(n) + '\n' + cs.replace('-','') + '\n')
        elif k == 2:
            fx = open('cluster_' + str(n) + '.fas')
            ct = fx.read().split('>')
            fx.close()
            final.write('>cluster_' + str(n) + '\n' + ct[1].split('\n',1)[1].replace('-', '').replace('\n', '') + '\n')
        elif k == 1:
            fx = open('cluster_' + str(n) + '.fas')
            ct = fx.read()
            fx.close()
            final.write('>cluster_' + str(n) + '\n' + ct.split('\n', 1)[1].replace('-', '').replace('\n', '') + '\n')

        n = n + 1
    final.close()

def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

def consensus(alig):
    fx= open(alig)
    ct = fx.read().split('>')
    fx.close()
    ls = []
    for x in ct[1:]:
        ls.append(x.split('\n',1)[1].replace('\n',''))
    lss = [list(i) for i in zip(*ls)]
    cons = ''
    for ssls in lss:

        cons = cons + Most_Common(ssls)
    fx = open(alig,'a')
    fx.write('>consensus\n' + cons + '\n')
    fx.close()
    return cons



def get5quantile(fichier):
    fx = open(fichier)
    ct = fx.readlines()
    fx.close()
    ls = []
    for x in ct:
        ls.append(float(x.split(' ')[2].replace('\n', '')))
    makefig2(ls)
    return np.quantile(ls, .25)


def getIdent(aln):
    fx = open(aln)
    op = fx.read().split('>')
    fx.close()
    tag = ''
    s1, s2 = op[1].split('\n',1)[1].replace('\n',''),op[2].split('\n',1)[1].replace('\n','')

    for x,y in zip(s1,s2):
        if x == '-' or y == '-':
            tag = tag + '*'
        else:
            tag = tag + 'N'
    i = 0
    while i <= len(tag):
        if tag[i] == 'N' and tag[i+1] == 'N':
            deb = i
            i = len(tag)
        i = i + 1
    i = len(tag) - 1
    while i > 0:
        if tag[i] == 'N' and tag[i-1] == 'N':
            fin = i + 1
            i = 0
        i = i - 1
    seq1,seq2 = s1[deb:fin], s2[deb:fin]

    n = 0.0
    for x,y in zip(seq1,seq2):
        if x == y:
            n = n + 1
    print(1.0 - (n / len(seq1)))
    return 1.0 - (n / len(seq1))

def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


def checklocus(ref,fastq):
    print('running checklocus')
    os.system('samtools faidx ' + ref )
    os.system('bwa index ' + ref)
    os.system('./minimap/minimap2 -ax map-ont ' + ref + ' ' + fastq + ' > Rs.sam')
    os.system('samtools view -b -F 4 Rs.sam > only.sam')
    os.system('samtools bam2fq only.sam | seqtk seq -a > temp0.fas')



def parse_blast(blastres):
    fx = open(blastres)
    bl = fx.readlines()
    fx.close()
    uniq = []
    for x in bl:
        if x.split('\t')[0].replace('_aln','') not in uniq:
            uniq.append(x.split('\t')[0].replace('_aln',''))
    dicbl = {}
    for x in uniq:
        i = 0
        while i <= len(bl):
            if x == bl[i].split('\t')[0].replace('_aln',''):
                locus = 0

                os.system('../edirect/efetch -db nucleotide -id ' + bl[i].split('\t')[1] + ' -format gp > outnt')
                fx = open('outnt').readlines()
                os.system('rm outnt')

                for y in fx:
                    if 'ORGA' in y:
                        sp = y.replace('ORGANISM  ','')
                        index = fx.index(y)
                        taxo = fx[index + 1].replace('\t','')
                        #print('sp and taxo', sp, taxo)
                        dicbl[x] = sp.replace('\n','')  + '\t' + taxo.replace('\n','') + '\t' + bl[i].split('\t')[1] + '\t' + bl[i].split('\t')[2]
                i = len(bl)
            i = i + 1
    return dicbl


def fin_out(dicCon, dicPol, dITS,d18S,d28S, dEF1):

    fr = open('racon.fasta')
    cr = fr.read().split('>')
    fr.close()
    dicr = {}
    for x in cr[1:]:
        dicr[x.split(' ')[0].replace('_aln', '')] = x.split('\n', 1)[1].replace('\n', '').upper()

    fx = open('Clusters_1st_round.fasta')
    ct = fx.read().split('>')
    fx.close()
    blst = []
    for x in ct[1:]:
        line = ""
        name = x.split('\n', 1)[0].replace('_aln', '')
        ft = open(name + '.fas')
        tp = ft.read().split('>')
        tp = len(tp) - 1
        ft.close()

        line = line + name + '\t' + str(tp) + '\t' + x.split('\n', 1)[1].replace('\n', '').upper() + '\t' + str(len(x.split('\n', 1)[1].replace('\n', '')))
        try:
            line = line + '\t' + dicCon[name]
        except:
            line = line + '\t-\t-'
        try:
            line = line + '\t' + name + '\t' + dicr[name] + '\t' + str(len(dicr[name])) + '\t' + dicPol[name] + '\t'
        except:
            line = line + '\t-\t-\t-\t-\t'

        if dITS != False:
            try:
                line = line + dITS[name]
            except:
                line = line + '\t-'

        elif dEF1 != False:
            try:
                line = line + '\t' + dEF1[name]
            except:
                line = line + '\t-'
        blst.append(line)
    return blst


def deci_out():
    decPH = []
    fb = open("decipher_out")
    dc = fb.read()
    fb.close()
    block1 = dc.split('"x"')[1].split('\n')
    block2 = dc.split('"x"')[2].split('\n')
    for x,y in zip(block1[1:-1],block2[1:-1]):
        pb,tgg = [],[]
        for k,l in zip(y.split('" "')[1].split(";"),x.split('" "')[1].split(";")     ):
            if "unidentified" not in k:
                tgg.append(k.replace('\n',''))
                pb.append(l.replace('\n',''))

        decPH.append(x.split('" "')[0].replace('"','') + '\t' + tgg[-1] + '\t' + pb[-1])
    return decPH





# MAIN #

#usage ==> python Ba.py 409_40790_1/barcode01 ITS  #ITS or EF1 will trigger the checklocus function
cmdline = sys.argv[0] + ' ' + sys.argv[1] + ' ' + sys.argv[2] + ' ' + sys.argv[3]
# sys.argv[3] is the minimum size of reads



logfile = open('_log.txt','w')
logfile.write(cmdline + '\n\n')
print(cmdline)
try:
    os.system('gunzip ' + sys.argv[1] + '/*.gz')
except:
    pass
os.system('rm ' + sys.argv[1] + '/' +  sys.argv[1].split('/')[-1] + '.fastq')
os.system('cat ' + sys.argv[1] + '/*.fastq > ' + sys.argv[1] + '/' +  sys.argv[1].split('/')[-1] + '.fastq')
inpt = sys.argv[1] + '/' + sys.argv[1].split('/')[-1] + '.fastq'



if sys.argv[2] == 'ITS':
    #os.system("cat ITS_RefSeq.fasta 18S_RefSeq.fasta 28S_RefSeq.fasta > rDNA_fungi.fas")
    #checklocus('rDNA_fungi.fas',inpt)
    checklocus('ITS_RefSeq.fasta',inpt)
elif sys.argv[2] == 'EF1':
    checklocus('EF1_fungi.fas',inpt)


logfile.write('Starting CleanT function\n')
logfile.write('Clipping vector and selected reads of expected size\n')

nb, tot, AllS = cleanT(int(sys.argv[3]))



logfile.write('Starting mafft\n')
os.system('mafft --auto temp1.fas > temp2.fas')
logfile.write('mafft alignment done\n')
logfile.write('Alignment is in temp2.fas\n\n')



logfile.write('Starting mothur - pairwise distances\n')
os.system("mothur/mothur '#dist.seqs(fasta= temp2.fas, calc=onegap, countends=F)'")
qt = get5quantile('temp2.dist')
logfile.write('Pairwise distances computation done\n')
logfile.write('Suggested cut-off = ' + str(qt) + '\n')
logfile.write('Pairwise distances are in temp2.dist\n')
logfile.write('Built Figure2.png with distribution of pairwise distances\n\n')

#qt = input('Please, enter a cutoff value:')  #0.42 tested with 400 to 1000 bp on barcode01
if sys.argv[2] == 'ITS':
    qt = 0.44
elif sys.argv[2] == 'EF1':
    qt = 0.35  #try 0.35 for fusarium 0.42-45 for ITS

logfile.write('Starting mothur - clustering\n')
os.system("mothur/mothur '#unique.seqs(fasta=temp1.fas)'")
os.system("mothur/mothur '#cluster(column=temp2.dist, count=temp1.count_table, cutoff=" + str(qt) + ")'")

logfile.write('Clustering done\n')
logfile.write('Building clusters alignments and consensus sequences - buildclusters function\n')
buildclusters(AllS)
logfile.write('buildclusters done\n')
logfile.write('Consensus sequences are in Clusters_1st_round.fasta\n')

logfile.write('Starting bwa - remapping reads to clusters\n')
os.system('bwa index Clusters_1st_round.fasta')
os.system('bwa mem -t 14 Clusters_1st_round.fasta ' + inpt + ' > mapping.sam')
logfile.write('bwa alignments done\n')


logfile.write('Starting polishing with racon\n')
os.system('./racon/build/bin/racon -m 8 -x -6 -g -8 -w 500 -t 14 ' + inpt + ' mapping.sam Clusters_1st_round.fasta > racon.fasta')
logfile.write('Polishing done\n')
logfile.write('Polished consensus sequences are in racon.fasta\n\n')

logfile.write('Querying NCBI nt with consensus sequences - blastn\n')
os.system('blastn -db nt -query Clusters_1st_round.fasta -num_threads 12 -outfmt 6 -max_target_seqs 5 > consensus_res.txt')
os.system('blastn -db nt -query racon.fasta -outfmt 6 -num_threads 12 -max_target_seqs 5 > polishing_res.txt')

if sys.argv[2] == 'ITS':
    os.system("Rscript decipher.R")
    os.system("makeblastdb -in ITS_RefSeq.fasta -dbtype nucl")
    #os.system("makeblastdb -in 18S_RefSeq.fasta -dbtype nucl")
    #os.system("makeblastdb -in 28S_RefSeq.fasta -dbtype nucl")
    os.system('blastn -db ITS_RefSeq.fasta -query racon.fasta -outfmt 6 -num_threads 12 -max_target_seqs 5 > ITS_RefSeq_res.txt')
    #os.system('blastn -db 18S_RefSeq.fasta -query racon.fasta -outfmt 6 -num_threads 12 -max_target_seqs 5 > 18S_RefSeq_res.txt')
    #os.system('blastn -db 28S_RefSeq.fasta -query racon.fasta -outfmt 6 -num_threads 12 -max_target_seqs 5 > 28S_RefSeq_res.txt')
    logfile.write('blastn done\n')
elif sys.argv[2] == 'EF1':
    os.system("makeblastdb -in EF1_RefSeq.fasta -dbtype nucl")
    os.system('blastn -db EF1_RefSeq.fasta -query racon.fasta -outfmt 6 -num_threads 12 -max_target_seqs 5 > EF1_RefSeq_res.txt')


logfile.write('Parsing blast results with efetch\n')
dicCon = parse_blast('consensus_res.txt')
dicPol = parse_blast('polishing_res.txt')


if sys.argv[2] == 'ITS':
    dicITS = parse_blast('ITS_RefSeq_res.txt')
elif sys.argv[2] == 'EF1':
    dicEF1 = parse_blast('EF1_RefSeq_res.txt')



logfile.write('efetch done\n\n')
logfile.write('Taxo results in consensus_res.txt and polishing_res.txt\n\n')

if sys.argv[2] == 'ITS':
    lsBlast = fin_out(dicCon, dicPol, dicITS, False, False, False)


    assign = deci_out()

elif sys.argv[2] == 'EF1':
    lsBlast = fin_out(dicCon, dicPol, False, False, False, dicEF1)

hits = open(sys.argv[1].split('/')[0] + '_' + sys.argv[1].split('/')[-1] + '_assignment_repport.xls', 'w')
hits.write(cmdline + '\tcutoff value=\t' + str(qt) + '\n')
hits.write("Consensus name\tNb. Of reads\tConsensus sequence\tLength\tFirst hit\tTaxonomy\tAccession number\t% similarity\tPolished consensus name\tPolished consensus sequence\tLenght\tFirst hit\tTaxonomy\tAccession number\t% similarity\n")


if sys.argv[2] == 'ITS':

    for m, k in zip(lsBlast,assign):
        hits.write(m + '\t' + k + '\n')
else:
    for m in lsBlast:
        hits.write(m + '\n')

hits.close()



os.system("rm racon.fas Cluster* cluster* mothur* temp* vectRes Rs.sam only.sam mapped* checklocus.txt decipher_out")

logfile.write('Output in Assignement_repport.xls\n')


logfile.close()
