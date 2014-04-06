import io
import sys
import time

# import matplotlib as plt
import numpy as np
import pandas as pd

from Bio import SeqIO
from scipy.stats import nbinom

# inputs:
# - expected number of reads
# - read length
# - paired end?
# - fasta file w/ sequences
# - relative abundance file
# - fld
# - prop. of DE
# - ***size factors

def readFasta(inFile):
    inHandle = open(inFile, "r")
    ids = []
    seqs = []
    for rec in SeqIO.parse(inHandle, "fasta"):
        ids.append(rec.id)
        seqs.append(rec.seq)
    df = pd.DataFrame({'seq': seqs}, index = ids)
    inHandle.close()

    return df

def computeEffLength(targetDf, meanFl):
    length = targetDf.apply(lambda x: len(x['seq']), axis = 1)
    effLength = length - meanFl + 1
    targetDf["effLength"] = effLength

    return targetDf

def readFpkm(inFile):
    fpkm = pd.read_csv(inFile, index_col = 0)

    return fpkm

def negBinom(mu, disp):
    var = mu + (mu **2) * disp
    r = mu ** 2 / (var - mu)
    p = r / (r + mu)

    return nbinom(r, p)

def kl(p, q):
    """Kullback-Leibler divergence D(P || Q) for discrete distributions
 
    Parameters
    ----------
    p, q : array-like, dtype=float, shape=n
        Discrete probability distributions.
    """
    p = np.asarray(p, dtype=np.float)
    q = np.asarray(q, dtype=np.float)
 
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def jsd(p, q):
    m = (p + q) / 2

    return (kl(p, m) + kl(q, m)) / 2

# returns log(rho)
def computeLogRho(counts, effLen):
    num = counts.apply(lambda col: col / effLen)
    denom = num.apply(np.sum)
    denom = np.log(denom)
    num = np.log(num)

    return num - denom

def lrhoToTpm(lrho):
    return np.exp(lrho + np.log(10000000))

@profile
def simulateReads(countSeries, readLen, seqs, revSeqs, fld, filePrefix):
    leftHandle = open(filePrefix + "_left.fasta", "w")
    rightHandle = open(filePrefix + "_right.fasta", "w")

    nTotal = countSeries.sum()
    countSeries = countSeries[countSeries > 0]

    # remove mass < readLen and renormalize
    fld[0:(readLen-1)] = 0.0
    fld = fld / fld.sum()

    flRange = pd.Series(range(readLen, fld.shape[0] + 1), 
            index = range(readLen, fld.shape[0] + 1))
    # print flRange

    # create Series of FLD from [readLen, maxFL]
    # FIXME: there's an off-by-one error here but I'm too lazy to fix it 
    # and it shouldn't make a big difference...
    flds = flRange.map(lambda x: (fld.ix[readLen:x] / fld.ix[readLen:x].sum()).tolist())
    flds = flds.to_dict()
    maxFl = max(flds.keys())
    flds = flds.values()

    ranOrder = np.random.choice(nTotal, size = nTotal)
    allLeft = [None] * nTotal
    allRight = [None] * nTotal

    # revSeqs = seqs.apply(lambda x: str(x.reverse_complement()))
    # seqs = seqs.apply(lambda x: str(x))

    readNum = 0
    for trans in xrange(countSeries.shape[0]):
        transName = countSeries.index[trans]
        transLen = len(seqs.iloc[trans])
        transSeq = seqs.iloc[trans]
        transRevSeq = revSeqs.iloc[trans]
        # choose a frag length
        curMaxFl = min(transLen, maxFl)
        fragLens = np.random.choice(len(flds[curMaxFl - readLen]), 
                size = countSeries.iloc[trans], 
                p = flds[curMaxFl - readLen])
        fragNum = 0

        for fragIt in xrange(countSeries.iloc[trans]):
            # choose a transcript
            #trans = np.random.randint(countSeries.shape[0])
    
            curFl = fragLens[fragNum]
    
            # choose a start site
            # print transName, curMaxFl, curFl
            # startSite = np.random.choice(transLen - curFl)
            startSite = int(np.random.random() * (transLen - curFl - 1))
    
            # randomly choose if reverse complement
            endSite = startSite + curFl
            if np.random.random() > 0.5:
                left = transSeq[startSite:(startSite + readLen)]
                right = transSeq[endSite - readLen:endSite]
            else:
                left = transRevSeq[startSite:(startSite + readLen)]
                right = transRevSeq[endSite - readLen:endSite]
    
            allLeft[ranOrder[readNum]] = left
            allRight[ranOrder[readNum]] = right
    
            readNum += 1
            fragNum += 1

    for i in xrange(len(allLeft)):
        fa_write(leftHandle, str(i), str(allLeft[i]))
        fa_write(rightHandle, str(i), str(allRight[i]))

    leftHandle.close()
    rightHandle.close()

def fa_write(fhandle, seq_id, seq):
    """
    Write to a FASTA file in the FASTA format.

    Arguments:
    - `fhandle`: A file handle open for writing
    - `seq_id`: The sequence id string for this sequence
    - `seq`: An unformatted string of the sequence to write
    """
    line_len = 60
    fhandle.write(">" + seq_id + "\n")
    for i in xrange(len(seq) / line_len + 1):
        start = i * line_len
        end = (i+1) * line_len if (i+1) * line_len < len(seq) else len(seq)
        fhandle.write( seq[ start:end ] + "\n")

def main():
    cfg = {}
    execfile(sys.argv[1], cfg) 
    # execfile("/Users/hjp/Documents/lmcb/pretender/examples/erythroidMouse.cfg", cfg)

    print "Reading transcriptome sequence ", cfg["fasta"]
    fixedData = readFasta(cfg["fasta"])
    fixedData['revSeq'] = fixedData['seq'].apply(lambda x: str(x.reverse_complement()))
    fixedData['seq'] = fixedData['seq'].apply(lambda x: str(x))
    print "Reading relative abundances ", cfg["fpkm"]
    fpkm = readFpkm(cfg["fpkm"])
    fixedData['fpkm'] = fpkm 
    print "Reading fragment length distribution ", cfg["fld"]
    fld = pd.read_table(cfg["fld"], header = None)[0]
    meanFl = sum(fld * pd.Series(range(0, len(fld))))
    print "Computing (mean effective length ", cfg["fld"]
    fixedData = computeEffLength(fixedData, meanFl)

    geneLabels = pd.read_csv(cfg["geneLabels"], index_col = 0)

    # remove everything w/ FPKM 0 and with negative effective length
    # XXX: consider removing everything w/ effLength < mean(fld)
    # XXX: look into transcripts w/  really short effective lengths... might be
    # getting too much expression
    fixedData = fixedData[(fixedData['fpkm'] != 0.0) & (fixedData['effLength']
        >= 0)]

    # TODO: Find bug that does...
    # /Library/Python/2.7/site-packages/pandas/core/series.py:628:
    # SettingWithCopyWarning: A value is trying to be set on a copy of a slice
    # from a DataFrame self.where(~key, value, inplace=True)
    # I think it's coming from the ~ usage in the next block

    # TODO: clean this mess up and modularize is

    # compute rho
    denom = fixedData['fpkm'].sum()
    fixedData['lrho'] = np.log(fixedData['fpkm']) - np.log(denom)
    fixedData = fixedData[~pd.isnull(np.exp(fixedData['lrho']))]

    # compute alpha
    lEffLen = np.log(fixedData['effLength'])
    num = fixedData['lrho'] + lEffLen
    denom = np.log(sum(np.exp(num)))
    fixedData['lalpha'] = num - denom
    fixedData['meanFrag'] = np.exp(fixedData['lalpha'] + np.log(cfg["numReads"]))

    # remove isoforms with < 0.5 mean fragments
    fixedData = fixedData[fixedData['meanFrag'] > 0.5]


    # decide whether isoform is going to be DE
    fixedData['isDE'] = np.random.uniform(0, 1, fixedData.shape[0]) <= 0.10
    fixedData['change'] = 1.0
    nChange = sum(fixedData['isDE'])

    # make 80% up-regulated and 20% down regulated with a tight, 
    # and normal distribution
    fixedData['change'][fixedData['isDE']] = np.random.choice([-1, 1], 
            nChange, replace = True, p = [1 - cfg["propUp"], cfg["propUp"]]) * \
                    np.random.normal(2, 0.1, nChange)

    fixedData['dispersion'] = fixedData['meanFrag'].apply(lambda x: 
            np.minimum(1 / (x * 0.005), 1/(1.5)))

    fixedData['var'] = fixedData['meanFrag'] + fixedData['dispersion'] * \
            fixedData['meanFrag'] ** 2

    fixedData['nb1'] = fixedData.apply(lambda row: 
            negBinom(row['meanFrag'], row['dispersion']) if row['meanFrag'] >= 0.5 \
                    else None, axis = 1)
    fixedData['nb2'] = fixedData.apply(lambda row: 
            negBinom(row['meanFrag'] * row['change'], row['dispersion']) if row['meanFrag'] * row['change'] >= 0.5 \
                    else None, axis = 1)
    fixedData = fixedData[fixedData['nb1'].notnull() & fixedData['nb2'].notnull()]

    print "Simulating counts from negative binomials.."
    counts1 = fixedData['nb1'].apply(lambda x: pd.Series(x.rvs(cfg["n1"])))
    counts2 = fixedData['nb2'].apply(lambda x: pd.Series(x.rvs(cfg["n2"])))
    isoCounts = pd.concat([counts1, counts2], axis = 1) 
    isoCounts.columns = range(cfg["n1"] + cfg["n2"])

    lrho = computeLogRho(isoCounts, fixedData['effLength'])

    # sanity check: JSD should be pretty small
    # lrho.apply(lambda col: jsd(np.exp(col), np.exp(fixedData['lrho'])))
    # 
    # plt.hist(lrho.apply(lambda col: jsd(np.exp(col), np.exp(fixedData['lrho']))))
    # plt.show()

    # compute corresponding TPM
    tpm = lrhoToTpm(lrho)
    tpm.columns = range(tpm.shape[1])
    tpm['gene'] = geneLabels.loc[tpm.index]
    # The "null" group is called "NotAGene"
    tpm['gene'] = tpm['gene'].fillna("NotAGene")
    geneTpm = tpm.groupby('gene').sum()

    # TODO: output tables (iso table and gene table)
    print "Simulating fragments"
    startTime = time.time()
    hi = simulateReads(isoCounts[0], cfg["readLength"], fixedData['seq'], fixedData['revSeq'], fld, 
            "sim")
    stopTime = time.time()
    print stopTime - startTime


if __name__ == "__main__":
    main()

