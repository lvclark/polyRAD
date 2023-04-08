#!/usr/bin/env python3
import isoloci_fun
import csv
import math
import argparse

## Script to process tag depth and sort haplotypes into isoloci ##
parser = argparse.ArgumentParser(description =
'''Process the output of one chunk from process_sam_multi.py.  Estimate
Hind/He both to filter markers and to determine which groups of tags could be
adjusted in terms of assignment of tags to isoloci.  For groups of tags needing
adjustment, analyze allele correlations to identify groups that putatively
belong to the same isolocus, then perform tabu search to find groups of tags
within Hind/He expections that minimize number of mutations from the reference
sequence.  Output read depth and assignment of alleles to loci, for import by
polyRAD.''',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("alignfile", nargs = '?',
                    help = "Path to alignment file output by process_sam_multi.py.")
parser.add_argument("depthfile", nargs = '?',
                    help = "Path to read depth file output by process_sam_multi.py.")
parser.add_argument("--out", "-o", nargs = '?', default = "",
                    help = "File name for output.  Generated from input files if not provided.")
parser.add_argument("--ploidy", "-p", nargs = '?', type = int, default = 2,
                    help = "Expected ploidy after splitting isoloci.")
parser.add_argument("--inbreeding", "-f", nargs = '?', type = float, default = 0.0,
                    help = "Inbreeding coefficient, ranging from 0 to 1.")
parser.add_argument("--expHindHe", "-e", nargs = '?', type = float,
                    help = "Expected value for Hind/He. Overrides ploidy and inbreeding if provided.")
parser.add_argument("--maxHindHe", "-m", nargs = '?', type = float,
                    help = "Maximum allowable value for Hind/He. Overrides ploidy and inbreeding if provided.")
parser.add_argument("--logfile", "-l", nargs = '?', default = "",
                    help = "Optional path to file where log should be written.")
parser.add_argument("--samples", "-s", nargs = '?', default = "",
                    help = "File listing names of samples to analyze and retain (one name per line).")
parser.add_argument("--maxalleles", "-a", nargs = '?', type = int, default = 500,
                    help = "Maximum number of alleles per locus.  Loci with more alleles are discarded.")

args = parser.parse_args()
alignfile = args.alignfile
depthfile = args.depthfile
outfile = args.out
ploidy = args.ploidy
inbreeding = args.inbreeding
logfile = args.logfile
samples_file = args.samples
max_alleles = args.maxalleles

# generate output file name if not provided
if outfile == "":
  outfile = alignfile.replace("align", "sorted")

# expected value for Hind/He
if args.expHindHe == None:
  expHindHe = (ploidy - 1)/ploidy * (1 - inbreeding)
else:
  expHindHe = args.expHindHe

if args.maxHindHe == None:
  # maximum tolerable Hind/He: halfway between this and the next ploidy, on a log scale
  p2 = ploidy * 2
  maxHindHe = math.exp((math.log((ploidy - 1)/ploidy) + math.log((p2 - 1)/p2) + 2 * math.log(1 - inbreeding))/2)
else:
  maxHindHe = args.maxHindHe

if expHindHe < 0 or expHindHe > 1:
  raise Exception("Expected Hind/He needs to be between zero and one.")
if expHindHe > maxHindHe:
  raise Exception("Max Hind/He needs to be greater than expected Hind/He.")

def ProcessRowGroup(alignrows, depthrows, nisoloci, thresh, expHindHe,
                    outwriter, logcon):
  '''Process two matching groups of rows showing alignment and depth for a
  group of tags corresponding to one set of alignment locations.'''
  # write marker being analyzed to log
  if logcon != None:
    logcon.write(" ".join(alignrows[0][:nisoloci]) + "\n")
  assert len(depthrows) > 1

  depths = [[int(d) for d in row[1:]] for row in depthrows] # integer depths
  # filter out any alleles without depth
  packed = [(dep, ar) for dep, ar in zip(depths, alignrows) if not all([d < 2 for d in dep])]
  if len(packed) < 2:
    if logcon != None:
      logcon.write("Insufficient read depth.\n")
    return None
  if len(packed) > max_alleles:
    if logcon != None:
      logcon.write("Too many alleles.\n")
    return None
  depths, alignrows = zip(*packed)

  # detect if fewer alignments than max isoloci
  nalign = sum([a != '' for a in alignrows[0][:nisoloci]])
  if nalign == 1: # only one isolocus
    hapAssign = [list(range(len(depths)))]
  else:
    # number of mutations from reference locations
    NM = [[int(row[nisoloci + 1 + i]) for row in alignrows] for i in range(nalign)]
    # workhorse function for making assignments
    hapAssign = isoloci_fun.TabuLocus(depths, NM, expHindHe, logcon = logcon)
  hindhe = isoloci_fun.HindHeByIsolocus(depths, hapAssign)

  # filter down to isoloci that don't exceed max Hind/He
  packed2 = [(hapAssign[i], alignrows[0][i]) for i in range(nalign)
             if hindhe[i] != None and hindhe[i] <= thresh]
  if len(packed2) == 0: # if no isoloci should be kept
    return None
  hapAssign, aligns = zip(*packed2)

  # get strings and starting positions for variable regions
  alNucs = [[] for i in range(len(hapAssign))]
  varpos = [0 for i in range(len(hapAssign))]
  for i in range(len(hapAssign)):
    # index the alignment position for this set of tags
    locind = alignrows[0][:nisoloci].index(aligns[i])
    cigars = [alignrows[j][nisoloci * 2 + 1 + locind] for j in hapAssign[i]]
    MDs = [alignrows[j][nisoloci * 3 + 1 + locind] for j in hapAssign[i]]
    tags = [alignrows[j][nisoloci] for j in hapAssign[i]]
    pos = int(aligns[i].split('-')[1])
    strand = aligns[i].split('-')[2]
    alNucs[i], varpos[i] = isoloci_fun.MakeAlleleStrings(tags, cigars, MDs, pos, strand)

  # write to file
  for i in range(len(aligns)):
    depthsOut = [depths[j] for j in hapAssign[i]]
    tagsOut = [alignrows[j][nisoloci] for j in hapAssign[i]]
    refNucs = alNucs[i].pop()
    # skip cases where tags varied in length and some became identical
    if any([an == '' for an in alNucs[i]]) or len(alNucs[i]) != len(set(alNucs[i])):
      continue
    [outwriter.writerow([aligns[i], varpos[i], alNucs[i][j], refNucs, tagsOut[j]] + depthsOut[j]) \
     for j in range(len(hapAssign[i]))]

  return None

# files need to already have tags in the same order, grouped by alignment location
# loop through markers
try:
  depthcon = open(depthfile, newline = '', mode = 'r')
  aligncon = open(alignfile, newline = '', mode = 'r')
  outcon = open(outfile, newline = '', mode = 'w')
  if logfile == "":
    logcon = None
  else:
    logcon = open(logfile, mode = 'w')
  depthreader = csv.reader(depthcon)
  alignreader = csv.reader(aligncon)
  outwriter = csv.writer(outcon)

  # header info from files
  depthheader = next(depthreader)
  ttd_samples = depthheader[1:]
  alignheader = next(alignreader)
  maxisoloci = sum([h.startswith("Alignment") for h in alignheader])
  if maxisoloci == 0:
    raise Exception("Columns starting with 'Alignment' not found in " + alignfile)
    
  # import samples if provided
  if samples_file == "":
    samples = ttd_samples
    sample_index = [i for i in range(len(samples))]
  else:
    with open(samples_file, mode = 'r') as mycon:
      samples = mycon.read().splitlines()
    if any([s not in ttd_samples for s in samples]):
      print("Names of samples not found in TagTaxaDist:")
      [print(s) for s in samples if s not in ttd_samples]
      raise Exception("Samples from {} not found in depth file.".format(samples_file))
    sample_index = [ttd_samples.index(s) for s in samples]
  
  # write header to output
  outwriter.writerow(["Marker", "Variable site", "Allele string", "Reference"] + \
  [depthheader[0]] + samples)

  newdepthrow = next(depthreader)
  newdepthrow = [newdepthrow[0]] + [newdepthrow[si + 1] for si in sample_index]
  currdepthrows = [newdepthrow]
  curralignrows = [next(alignreader)]

  rowcount = 0 # for testing, to prevent going thru whole file

  for row in alignreader:
    newalignrow = row
    newdepthrow = next(depthreader)
    newdepthrow = [newdepthrow[0]] + [newdepthrow[si + 1] for si in sample_index]
    tag = newdepthrow[0] # tag sequence
    assert newalignrow[maxisoloci] == tag
    if newalignrow[:maxisoloci] == curralignrows[0][:maxisoloci]: # same alignment
      currdepthrows.append(newdepthrow)
      curralignrows.append(newalignrow)
    else: # new alignment; process last one and start new
      ProcessRowGroup(curralignrows, currdepthrows, maxisoloci, maxHindHe,
                      expHindHe, outwriter, logcon)
      currdepthrows = [newdepthrow]
      curralignrows = [newalignrow]
    rowcount += 1
    # if rowcount > 5000: # for testing only
    #    break
    if rowcount % 1000 == 0:
      print(rowcount)
  ProcessRowGroup(curralignrows, currdepthrows, maxisoloci, maxHindHe,
                  expHindHe, outwriter, logcon)
finally:
  depthcon.close()
  aligncon.close()
  outcon.close()
  if logcon != None:
    logcon.close()
