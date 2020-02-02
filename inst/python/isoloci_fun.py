from statistics import mean, StatisticsError
from numpy.random import choice
from scipy.stats import kendalltau
from copy import deepcopy
import math
import sys
import re

if sys.version_info.major < 3:
  raise Exeption("Python 3 required.")

# functions for sorting out isoloci

def GiniSimpson(counts, N = None):
  "Estimate the Gini-Simpson index, not corrected for size, on a list of integers."
  if N == None:
    N = sum(counts)
  if N == 0:
    return None
  for c in counts: # check to see if it should be 0, before doing calc.
    if c == N:
      return 0.0
    elif c != 0:
      break
  freq = [c/N for c in counts]
  return 1.0 - sum([f**2 for f in freq])

def HindHe(countsmat):
  '''For one marker, estimate Hind/He.
  countsmat is a list of lists, with the first dimension representing alleles
  and the second dimension representing taxa. Each cell is read depth.
  Hind/He will be corrected for size on a per-taxon basis by multiplying by
  N/(N-1), then averaged across taxa.'''
  if len(countsmat) < 2:
    return None

  countsmatT = list(zip(*countsmat)) # transpose the matrix
  depthByInd = [sum(x) for x in countsmatT]
  nind = len(depthByInd)
  if all([d == 0 for d in depthByInd]):
    return None

  depthRatios = [[c / depthByInd[i] for c in countsmatT[i]] for i in range(nind) if depthByInd[i] > 0]
  meanDepthRatios = [sum(x)/len(x) for x in zip(*depthRatios)] # mean by allele
  He = GiniSimpson(meanDepthRatios, N = 1)
  assert He > 0

  HindHeByInd = [GiniSimpson(countsmatT[i], N = depthByInd[i]) * \
                 depthByInd[i] / (depthByInd[i] - 1)/ He for i in range(nind) if depthByInd[i] > 1]

  m = sum(HindHeByInd)/len(HindHeByInd)
  return m

def InitHapAssign(NMmat):
  '''Generate initial assignments of haplotypes to loci based on number of
  mutations between haplotype and the locus.'''
  nloc = len(NMmat)
  nhap = len(NMmat[0])
  hapAssign = [[] for i in range(nloc)] # to store output

  for h in range(nhap):
    theseNM = [nm[h] for nm in NMmat]
    minNM = min(theseNM)
    bestLoc = [l for l in range(nloc) if theseNM[l] == minNM]
    if len(bestLoc) > 1:
      bestLoc = choice(bestLoc, 1)
    bestLoc = bestLoc[0] # convert list to number
    hapAssign[bestLoc].append(h) # add haplotype to locus

  return hapAssign

def HindHeByIsolocus(countsmat, hapAssign):
  '''For a given set of assignments of haplotypes to isoloci, estimate Hind/He
  for each isolocus.  countsmat and hapAssign are as defined above.'''
  splitcounts = [[countsmat[h] for h in isolocus] for isolocus in hapAssign]
  return [HindHe(c) for c in splitcounts]

def MeanNMperLoc(NMmat, hapAssign):
  '''Get mean number of mutations per locus across all haplotypes versus their
  assigned locus.'''
  nLoc = len(NMmat)
  nHap = len(NMmat[0])
  out = sum([NMmat[L][h] for L in range(nLoc) for h in hapAssign[L]])/nHap
  return out

def IndexHapAssign(hapAssign):
  '''Generate an index for a given set of assignments of haplotypes to isoloci,
  so that we can track solutions that have already been examined.'''
  nLoc = len(hapAssign)
  return sum([i * nLoc ** h for i in range(nLoc) for h in hapAssign[i]])

def AlleleAssociations(countsmat):
  '''Generate a square matrix of p-values for alleles being negatively
  associated with each other.'''
  nHap = len(countsmat)
  nInd = len(countsmat[0])
  outP = [[1.0 for i in range(nHap)] for j in range(nHap)]

  for h1 in range(nHap - 1):
    # tot1 and tot2 are per individual total read depth, omitting h1 and h2, respectively
    tot1 = [sum([countsmat[h][i] for h in range(nHap) if h != h1]) for i in range(nInd)]
    for h2 in range(h1 + 1, nHap):
      tot2 = [sum([countsmat[h][i] for h in range(nHap) if h != h2]) for i in range(nInd)]
      # depth ratios that ignore the other allele being considered
      rat1 = [countsmat[h1][i]/tot2[i] if tot2[i] > 0 else 0.0 for i in range(nInd)]
      rat2 = [countsmat[h2][i]/tot1[i] if tot1[i] > 0 else 0.0 for i in range(nInd)]
      # omit ones that are both zero (missing data)
      rat1, rat2 = zip(*[(r1, r2) for r1, r2, t1, t2 in zip(rat1, rat2, tot1, tot2) if t1 > 0 or t2 > 0])
      # skip if there is too much missing data
      if len(rat1) < 10:
        continue
      # skip if there is no variance in either allele
      if len(set(rat1)) == 1 or len(set(rat2)) == 1:
        continue
      # perform test for association
      kout = kendalltau(rat1, rat2, nan_policy = 'raise', method = 'asymptotic')
      # convert p-value to one-tailed
      if kout[0] <= 0:
        pout = kout[1]/2
      else:
        pout = 1 - kout[1]/2
      # add to matrix
      outP[h1][h2] = pout
      outP[h2][h1] = pout
  return(outP)

def GrpAllowedAligns(grp, NMmat):
  '''Because some haplotypes may not have actually aligned to some locations,
  find allowable locations for a group, where every haplotype aligned.'''
  nalign = len(NMmat)
  dummy_NM = 999 # NM if there was really no alignment
  return [all([NMmat[a][h] != dummy_NM for h in grp]) for a in range(nalign)]

def GroupByAlAssociations(countsmat, NMmat, expHindHe, startP = 0.1):
  '''Find groups of alleles that are significantly negatively associated, and
  don't exceed the expected value of Hind/He.'''
  nHap = len(countsmat)
  pvals = AlleleAssociations(countsmat)
  currP = startP

  while True:
    grps = []
    for h1 in range(nHap - 1):
      for h2 in range(h1 + 1, nHap):
        if pvals[h1][h2] > currP:
          continue # doesn't meet threshold, don't add to group
        if not any(GrpAllowedAligns({h1, h2}, NMmat)):
          continue # there are no subgenomes to which both aligned
        if len(grps) == 0: # first group
          grps.append({h1, h2})
          continue
        # determine if it is already in any groups
        h1grp = [i for i in range(len(grps)) if h1 in grps[i]]
        h2grp = [i for i in range(len(grps)) if h2 in grps[i]]
        assert len(h1grp) < 2
        assert len(h2grp) < 2
        if len(h1grp) == 0 and len(h2grp) == 0:
          grps.append({h1, h2}) # new group
        elif len(h1grp) == 1 and len(h2grp) == 0 and \
        any(GrpAllowedAligns(grps[h1grp[0]] | {h2}, NMmat)):
          grps[h1grp[0]].add(h2) # add to existing group
        elif len(h1grp) == 0 and len(h2grp) == 1 and \
        any(GrpAllowedAligns(grps[h2grp[0]] | {h1}, NMmat)):
          grps[h2grp[0]].add(h1) # add to existing group
        elif len(h1grp) == 1 and len(h2grp) == 1 and h1grp[0] != h2grp[0] and \
        any(GrpAllowedAligns(grps[h1grp[0]] | grps[h2grp[0]], NMmat)):
          # merge groups
          grps[h1grp[0]].update(grps[h2grp[0]])
          grps.pop(h2grp[0])
    if(len(grps) == 0):
      break # no groups could be made
    # Test Hind/He for these groups
    hindheOK = [HindHe([countsmat[h] for h in g]) <= expHindHe for g in grps]
    if all(hindheOK):
      break
    currP = currP / 10
  return([grps, currP])

def AdjustHapAssignByAlAssociations(grps, hapAssign, NMmat):
  '''Adjust hapAssign if necessary so that each group that was made based on
  negative associations between alleles is in just one haplotype group.'''
  nLoc = len(hapAssign)
  for grp in grps:
    numPerHA = [sum([g in ha for g in grp]) for ha in hapAssign]
    assert sum(numPerHA) == len(grp)
    haInGrp = [i for i in range(nLoc) if numPerHA[i] > 0]
    if len(haInGrp) == 1:
      continue # no rearrangement needed
    # get alignment locations allowable for this group
    allowed = GrpAllowedAligns(grp, NMmat)
    assert any(allowed)
    if any([allowed[i] for i in haInGrp]):
      # go to the isolocus where most of these are, or a random one.
      maxPerHA = max([numPerHA[i] for i in haInGrp if allowed[i]])
      matchmax = [i for i in haInGrp if numPerHA[i] == maxPerHA and allowed[i]]
      if len(matchmax) == 1:
        targetLoc = matchmax[0]
      else:
        targetLoc = choice(matchmax, size = 1)[0]
    else:
      # go to the allowed isolocus with the best sequence match
      haOK = [i for i in range(nLoc) if allowed[i]]
      NMtot = [sum([NMmat[i][h] for h in grp]) for i in haOK]
      minNM = min(NMtot)
      matchmin = [haOK[i] for i in range(len(haOK)) if NMtot[i] == minNM]
      if len(matchmin) == 1:
        targetLoc = matchmin[0]
      else:
        targetLoc = choic(matchmin, size = 1)[0]

    # rearrange haplotypes
    [hapAssign[i].remove(h) for i in haInGrp for h in grp if h in hapAssign[i]]
    hapAssign[targetLoc].extend(grp)

  return hapAssign

def HindHeExcess(hindhe, expHindHe):
  '''Return the mean amount by which per-isolocus Hind/He exceeds the expected
  value.'''
  excess = [max([0, h - expHindHe]) for h in hindhe if h != None]
  if len(excess) == 0:
    return 0.0 # all isoloci fixed; want to keep this solution
  else:
    return mean(excess)

def FindNeighbors(hapAssign, corrgrps, NMmat, tabu):
  '''For the Tabu Search algorithm, find neighboring solutions to the current
  one that have not recently been examined.  corrgrps is a list of groups of
  haplotypes that can be swapped together, including all individual haplotypes
  that were not placed in groups.'''
  nLoc = len(hapAssign)
  assert sum([len(grp) for grp in corrgrps]) == sum([len(ha) for ha in hapAssign])
  haList = []
  for grp in corrgrps:
    for loc in range(nLoc):
      if not all([g in hapAssign[loc] for g in grp]):
        hapAssign_new = deepcopy(hapAssign)
        [[ha.remove(g) for g in grp if g in ha] for ha in hapAssign_new]
        hapAssign_new[loc].extend(grp)
        haInd = IndexHapAssign(hapAssign_new)
        allowed = [GrpAllowedAligns(ha, NMmat) for ha in hapAssign_new]
        if haInd not in tabu and \
        all([allowed[i][i] for i in range(nLoc)]):
          haList.append(hapAssign_new)
  return haList

def TabuLocus(countsmat, NMmat, expHindHe, \
reps = 25, maxTabu = 5, corrstartP = 0.01, logcon = None):
  '''A Tabu Search to try to optimize first Hind/He, then NM from reference,
  while keeping together groups of alleles that are negatively associated.'''
  nHap = len(countsmat)
  nLoc = len(NMmat)
  hapAssign = InitHapAssign(NMmat) # initial assignment of haplotypes to isoloci
  hindhe = HindHeByIsolocus(countsmat, hapAssign) # doesn't get updated
  NM_mean = MeanNMperLoc(NMmat, hapAssign)        # doesn't get updated
  if logcon != None:
    logcon.write("Initial Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
    logcon.write("Initial average NM: {}\n".format(NM_mean))
  # if already fixed or within Hind/He expectations at each isolocus, don't do search
  if all([h == None or h < expHindHe for h in hindhe]):
    return hapAssign

  # get groups based on allele correlations, and adjust hapAssign if needed
  corrgrps, corrP = GroupByAlAssociations(countsmat, NMmat, expHindHe, startP = corrstartP)
  hapAssign = AdjustHapAssignByAlAssociations(corrgrps, hapAssign, NMmat)
  hindhe = HindHeByIsolocus(countsmat, hapAssign)
  # if everything ok after adjusting by corrgrps, don't do search
  if all([h == None or h < expHindHe for h in hindhe]):
    return hapAssign
  # expand corrgrps to add individual alleles
  if len(corrgrps) == 0:
    corrgrps = [{h} for h in range(nHap)]
  else:
    corrgrps.extend([{h} for h in range(nHap) if all([h not in grp for grp in corrgrps])])
  nGrp = len(corrgrps)

  # set up tabu list
  tabu = [0 for i in range(maxTabu)]
  tabu[0] = IndexHapAssign(hapAssign)
  tabuIndex = 1

  # set aside objects to hold best solution
  hapAssign_best = deepcopy(hapAssign)
  hindhe_excess_best = HindHeExcess(hindhe, expHindHe)
  NM_mean_best = MeanNMperLoc(NMmat, hapAssign)
  best_rep = 0 # rep where best solution was found

  # tabu search algorithm
  for rep in range(reps):
    # get all non-tabu neighbors and choose the best
    neighbors = FindNeighbors(hapAssign, corrgrps, NMmat, tabu)
    if len(neighbors) == 0:
      break
    hindhe_neighbors = [HindHeByIsolocus(countsmat, ha) for ha in neighbors]
    hindhe_excess_neighbors = [HindHeExcess(hh, expHindHe) for hh in hindhe_neighbors]
    min_excess = min(hindhe_excess_neighbors)
    min_neighbors = [neighbors[i] for i in range(len(neighbors)) \
                     if hindhe_excess_neighbors[i] == min_excess]
    if len(min_neighbors) == 1:
      hapAssign = min_neighbors[0]
    else:
      NM_neighbors = [MeanNMperLoc(NMmat, ha) for ha in min_neighbors]
      min_NM = min(NM_neighbors)
      min_neighbors = [min_neighbors[i] for i in range(len(min_neighbors)) \
                       if NM_neighbors[i] == min_NM]
      if len(min_neighbors) == 1:
        hapAssign = min_neighbors[0]
      else:
        hapAssign = min_neighbors[int(choice(range(len(min_neighbors)), size = 1))]
    # update the tabu list
    tabu[tabuIndex % maxTabu] = IndexHapAssign(hapAssign)
    tabuIndex += 1
    # update the best solution if appropriate
    min_NM = MeanNMperLoc(NMmat, hapAssign)
    if min_excess < hindhe_excess_best or (min_excess == hindhe_excess_best and min_NM < NM_mean_best):
      hapAssign_best = deepcopy(hapAssign)
      hindhe_excess_best = min_excess
      NM_mean_best = min_NM
      best_rep = rep

  if logcon != None:
    logcon.write("Final Hind/He: {}\n".format(" ".join([str(h) for h in HindHeByIsolocus(countsmat, hapAssign_best)])))
    logcon.write("Final average NM: {}\n".format(NM_mean_best))
    logcon.write("Rep where best solution found: {}\n".format(best_rep))

  return hapAssign_best

def SplitCigar(cigar):
  '''Split a CIGAR string into its components to process one at a time.'''
  return re.findall("\d+[MIDNSHP=X]", cigar)

def SplitMD(MD):
  '''Split an MD string into its components to process one at a time.
  Goal: get numbers, letters, or ^ plus letters.'''
  out = re.findall("(\d+|[ACGTN]+|\^[ACGTN]+)", MD)
  out = [sub if re.match("\d+", sub) == None else int(sub) for sub in out]
  return out

def CigarsToNucPos(tags, cigars, pos, strand):
  '''Take a set of CIGAR strings, along with the position of the leftmost
  nucleotide of a tag and to what strand that tag aligned, and get the genomic
  position of each nucleotide in each tag.'''
  nucpos = [[0 for p in range(len(t))] for t in tags]
  for ti in range(len(tags)):
    # split cigar string into values of interest
    cigsplit = SplitCigar(cigars[ti])
    currpos = 0
    currnuc = 0
    for cs in cigsplit:
      n = int(cs[:-1])
      let = cs[-1]
      if let in 'M=X':
        for i in range(n):
          nucpos[ti][currnuc] = currpos
          currpos += 1
          currnuc += 1
      if let in 'DN':
        currpos += n
      if let in 'IS':
        for i in range(n):
          nucpos[ti][currnuc] = currpos
          currnuc += 1
  assert all([len(tags[ti]) == len(nucpos[ti]) for ti in range(len(tags))])
  # number from reference rather than starting at zero
  if strand == 'top':
    nucpos = [[p + pos for p in np] for np in nucpos]
  else:
    ends = [np[-1] for np in nucpos]
    nucpos = [[p - ends[ti] + pos for p in nucpos[ti]] for ti in range(len(nucpos))]
  return nucpos

def PadPosition(tagnucs):
  '''Taking a set of nucleotides matching one position, pad them to represent
  insertions or deletions.'''
  # use - for deletions, or missing sequence trimmed from read
  tagnucs = ['-' if tn == '' else tn for tn in tagnucs]
  # find where insertions happened, and represent with .
  nuclens = [len(tn) for tn in tagnucs]
  maxlen = max(nuclens)
  if maxlen > 1:
    tagnucs = ['.' * (maxlen - nuclens[ti]) + tagnucs[ti] for ti in range(len(tagnucs))]
  return tagnucs

def RecreateReference(tags, cigars, MDs):
  '''Recreate the reference tag from CIGAR and MD information from the SAM
  file.'''
  # Determine the length of the reference spanned by each tag, and pick the longest.
  cigsplits = [SplitCigar(cig) for cig in cigars]
  cigsplits = [[s for s in cs if s[-1] not in set('ISHP')] for cs in cigsplits] # remove insertions
  align_lengths = [sum([int(s[:-1]) for s in cs]) for cs in cigsplits] # get total
  longest = max(align_lengths)
  keep = [i for i in range(len(tags)) if align_lengths[i] == longest]
  
  tags = [tags[i] for i in keep]
  cigars = [cigars[i] for i in keep]
  MDs = [MDs[i] for i in keep]
  
  # See if there are any tags that are already known to be the reference
  refs = [i for i in range(len(tags)) if re.match("\d+M$", cigars[i]) != None and \
  re.match("\d$", MDs[i]) != None]

  if len(refs) > 1:
    raise Exception("Multiple tags seem to be the reference.")
  if len(refs) == 1:
    i = refs[0]
    if MDs[i] != cigars[i][:-1]:
      raise Exception("MD and CIGAR don't match.")
    return (tags[i], cigars[i], MDs[i])

  # Otherwise, use the first tag to recreate the reference.
  # Use CIGAR string to remove any insertions.
  cigsplit = SplitCigar(cigars[0])
  thistag = ""
  currpos = 0
  for cig in cigsplit:
    n_nuc = int(cig[:-1])
    if cig[-1] in {'M', '=', 'X', 'I', 'S'}:
      if cig[-1] in {'M', '=', 'X'}:
        thistag = thistag + tags[0][currpos:(currpos + n_nuc)]
      currpos += n_nuc
  # Use MD string to fix any mismatches and deletions
  MDsplit = SplitMD(MDs[0])
  reftag = ""
  currpos = 0
  for md in MDsplit:
    if isinstance(md, int):  # add matching part of tag
      reftag = reftag + thistag[currpos:(currpos + md)]
      currpos += md
    elif md.startswith("^"): # add deletion from reference
      reftag = reftag + md[1:]
    else:                    # add reference nucleotide from mismatch
      reftag = reftag + md
      currpos += len(md)
  cigout = "{}M".format(len(reftag))
  mdout = "{}".format(len(reftag))
  return (reftag, cigout, mdout)

def MakeAlleleStrings(tags, cigars, MDs, pos, strand):
  '''Taking tag sequences, CIGAR strings, a position for the alignment
  starting at the cut site, and a strand (top or bot), make a set of
  strings just showing the variable portion of the tags, and also return a
  position for the beginning of those strings in the reference genome.'''
  # Note: if ALL tags have an insertion with respect to the reference, it
  # won't be obvious from the output of this function. ### Fixing

  assert strand == 'top' or strand == 'bot'
  if strand == 'bot':
    # Get reverse complement for bottom strand
    trans = {ord('A'): 'T', ord('C'): 'G', ord('G'): 'C', ord('T'): 'A'}
    tags = [t.translate(trans)[::-1] for t in tags]

  # Add the reference tag to the list
  ref = RecreateReference(tags, cigars, MDs)
  tags.append(ref[0])
  cigars.append(ref[1])
  MDs.append(ref[2])

  # Get position for each nucleotide
  nucpos = CigarsToNucPos(tags, cigars, pos, strand)
  # set up allele strings
  alstrings = ['' for t in tags]
  # get starting and stopping position to examine
  start = min([min(np) for np in nucpos])
  end = max([max(np) for np in nucpos])
  # starting and stopping positions of the variable region (to update)
  foundvar = False
  endvar = 0
  outpos = start
  # go through and build allele strings
  for p in range(start, end + 1):
    nucind = [[pi for pi in range(len(np)) if np[pi] == p] for np in nucpos]
    tagnucs = [tags[ti][nucind[ti][0]:(nucind[ti][-1] + 1)] \
               if len(nucind[ti]) > 0 else '' for ti in range(len(tags))]
    if foundvar: # add nucleotides to existing strings
      tagnucs = PadPosition(tagnucs)
      alstrings = [alstrings[ti] + tagnucs[ti] for ti in range(len(tags))]
      if len(set(tagnucs)) > 1:
        endvar = len(alstrings[0])
    elif len(set(tagnucs)) > 1: # found first variable site
      foundvar = True
      tagnucs = PadPosition(tagnucs)
      alstrings = tagnucs
      outpos = p
      endvar = len(alstrings[0])
  alstrings = [st[:endvar] for st in alstrings]

  return alstrings, outpos
