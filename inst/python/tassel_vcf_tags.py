#!/usr/bin/env python3

import argparse
import bisect
import re
import sys

if sys.version_info[0] < 3:
  raise Exception("Python 3 required.")

parser = argparse.ArgumentParser(description = "For polyRAD allele names imported from a TASSEL VCF, import tag sequences from a SAM file.")
parser.add_argument('-a', '--alleles', help = "File containing allele names.",
                    required = True)
parser.add_argument('-s', '--sam',
                    help = "One or more SAM files in which to search for tags.",
                    required = True, nargs = '+')
parser.add_argument('-o', '--out', help = "File name for output.",
                    required = True)

args = parser.parse_args()

alfile = args.alleles
samfiles = args.sam
outfile = args.out

# Import alleles
print("Importing alleles...")
with open(alfile, mode = 'r') as mycon:
  alnames = mycon.readlines()

nal = len(alnames)
chrom = [a.split("_")[0] for a in alnames]
pos = [int(a.split("_")[1]) for a in alnames]
alnuc = [a.split("_")[2].strip() for a in alnames]
chrom, pos, alnuc = zip(*sorted(zip(chrom, pos, alnuc)))

# allele dictionary - by chromosome, have position, alleles, and tags
alleleDict = {c: [(), [], []] for c in set(chrom)}

print("Building allele dictionary...")
for c in alleleDict.keys():
  chromindex = range(bisect.bisect_left(chrom, c), bisect.bisect_right(chrom, c))
  thispos, thisalnuc = zip(*[(pos[i], alnuc[i]) for i in chromindex])
  alleleDict[c][0] = tuple(sorted(set(thispos)))
  npos = len(alleleDict[c][0])
  albypos = [[] for i in range(npos)]
  tagsbypos = [[] for i in range(npos)]
  m = 0
  for i in range(len(thispos)):
    assert thispos[i] == alleleDict[c][0][m]
    albypos[m].append(thisalnuc[i])
    tagsbypos[m].append([])
    if i < len(thispos) - 1 and thispos[i] != thispos[i+1]:
      m += 1
  alleleDict[c][1] = albypos
  alleleDict[c][2] = tagsbypos

# helper functions for SAM files
def parseCigar(cigarstr):
  '''Split up a CIGAR string into a tuple of lengths and a tuple of codes
  indicating the meaning of those lenghts.'''
  nums = re.split('[MIDNSHP=X]', cigarstr)
  nums = [int(n) for n in nums if n != '']
  codes = [c for c in re.split('[0-9]', cigarstr) if c != '']
  return (tuple(nums), tuple(codes))

def matchAllele(tag, alleles, cigar, tagpos, alpos):
  '''For a tag (which may have already been reverse complemented) and a set
  of possible alleles, return the index of the matching allele, or -1 if there
  are none.'''
  # get reference position for every nucleotide in the tag
  tagref = []
  nextpos = tagpos
  for i in range(len(cigar[0])):
    if cigar[1][i] in {'M', '=', 'X'}:
      tagref.extend(list(range(nextpos, nextpos + cigar[0][i])))
    if cigar[1][i] in {'M', '=', 'X', 'D', 'N'}:
      nextpos += cigar[0][i]
    if cigar[1][i] in {'I', 'S'}:
      tagref.extend([-1 for j in range(cigar[0][i])])
  assert len(tagref) == len(tag)
  # get reference position for each nucleotide in the allele
  allen = len(alleles[0])
  assert all([len(a) == allen for a in alleles])
  assert alpos + allen - 1 <= max(tagref)
  alref = range(alpos, alpos + allen)
  # locations of variable sites within alleles
  varpos = [i for i in range(allen) if len(set([a[i] for a in alleles])) > 1]
  # determine which alleles could match tag
  possibleAl = list(range(len(alleles)))
  for i in varpos:
    if alref[i] not in tagref:
      return -1 # deletions in tag sequence, can't ID allele
    thisnuc = tag[tagref.index(alref[i])]
    possibleAl = [a for a in possibleAl if alleles[a][i] in {thisnuc, 'N'}]
    if len(possibleAl) == 0:
      return -1
  if len(possibleAl) != 1:
    possAlNuc = [alleles[a] for a in possibleAl]
    if len(set(possAlNuc)) < len(possAlNuc):
      print(tag)
      print(possAlNuc)
      raise Exception("Identical alleles found.  Be sure to only list each polyRAD allele once in the input file.")
    assert not all([p in tagref for p in alref])
    return -1 # otherwise, there were probably multiple matches due to a deletion
  return possibleAl[0]

def reverseComplement(sequence):
    '''Make the reverse complement of a nucleotide sequence.'''
    x = {ord('A'): 'T', ord('C'): 'G', ord('G'): 'C', ord('T'): 'A'}
    return(sequence.translate(x)[::-1])

# read and process SAM files
print("Reading SAM files...")
matchcount = 0
for sf in samfiles:
  linecount = 0
  with open(sf, mode = 'r') as samcon:
    for line in samcon:
      linecount += 1
      if linecount % 100000 == 0:
        print("{} {}".format(sf, linecount))
      if line[0] == '@':
        continue # skip header
      row = line.split('\t')
      flag = int(row[1])
      if flag & 4 == 4:
        continue # no alignment
      # chromosome name, formatted for TASSEL
      c = 'S' + row[2].upper().replace('CHR', '').replace('CHROMOSOME', '')
      if c not in alleleDict.keys():
        continue # no alleles on this chromosome
      leftpos = int(row[3]) # leftmost mapping position
      cigar = parseCigar(row[5])
      ciglen = len(cigar[0])
      # rightmost mapping position
      rightpos = leftpos + sum([cigar[0][i] for i in range(ciglen) if cigar[1][i] in {'M', 'D', 'N', '=', 'X'}]) - 1
      # determine if this tag overlaps an allele of interest
      posindex = bisect.bisect_left(alleleDict[c][0], leftpos)
      while posindex < len(alleleDict[c][0]) and \
      alleleDict[c][0][posindex] + len(alleleDict[c][1][posindex][0]) - 1 <= rightpos:
        # search for matches
        thismatch = matchAllele(row[9], alleleDict[c][1][posindex], cigar, leftpos, alleleDict[c][0][posindex])
        if thismatch != -1:
          # store tag
          if flag & 16 == 16:
            alleleDict[c][2][posindex][thismatch].append(reverseComplement(row[9]))
          else:
            alleleDict[c][2][posindex][thismatch].append(row[9])
          matchcount += 1
          if matchcount % 100000 == 0:
            print("{} tags matched to alleles".format(matchcount))
        posindex += 1 # crawl along to see if any other alleles overlap

# export alleles and tags
print("Exporting results...")
with open(outfile, mode = 'w') as outcon:
  for c in sorted(alleleDict.keys()):
    for m in range(len(alleleDict[c][0])):
      for a in range(len(alleleDict[c][1][m])):
        outcon.write("{}_{}_{}\t{}\n".format(c, alleleDict[c][0][m], alleleDict[c][1][m][a], ";".join(alleleDict[c][2][m][a])))
