#!/usr/bin/env python3

import csv
import re
import argparse
import bisect
import sys

if sys.version_info.major < 3:
    raise Exeption("Python 3 required.")

parser = argparse.ArgumentParser(description =
'''Process a SAM file reporting multiple alignments, and output one or more CSV
files that list tag sequences organized by sets of aligment locations.  A
TagTaxaDist file output by TASSEL is also processed to extract read depth for
the same tags.  Filtering by missing data rate is performed.  The output files
can then be used for sorting tags into isoloci.''',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("sam", nargs = '?', help = "Path to SAM file.")
parser.add_argument("ttd", nargs = '?', help = "Path to TagTaxaDist file.")
parser.add_argument("out", nargs = '?', help = "Base file name for output.")
parser.add_argument("--subgenomes", "-g", nargs = '?', type = int, default = 2,
                    help = "Number of subgenomes, i.e. maximum number of alignments expected per tag.")
parser.add_argument("--chunks", "-c", nargs = '?', type = int, default = 1,
                    help = "Number of files to split the output into.")
parser.add_argument("--min_ind_with_reads", "-m", nargs = '?', type = int,
                    default = 100,
                    help = "Minimum number of samples with sequencing reads needed to retain a tag set.")
parser.add_argument("--samples", "-s", nargs = '?', default = "",
                    help = "File listing names of samples to retain (one name per line).")
args = parser.parse_args()
mysam = args.sam
myttd = args.ttd
outbase = args.out
maxalign = args.subgenomes
nchunks = args.chunks
min_ind_with_reads = args.min_ind_with_reads
samples_file = args.samples

# Read in samples to keep, if provided, and confirm that they are
# present in TTD file.
with open(myttd, mode = 'r') as mycon:
  ttd_samples = next(mycon).split()[1:]
if samples_file == "":
  samples = ttd_samples
  sample_index = [i for i in range(len(samples))]
else:
  with open(samples_file, mode = 'r') as mycon:
    samples = mycon.read().splitlines()
  if any([s not in ttd_samples for s in samples]):
    print("Names of samples not found in TagTaxaDist:")
    [print(s) for s in samples if s not in ttd_samples]
    raise Exception("Samples from {} not found in TagTaxaDist.".format(samples_file))
  sample_index = [ttd_samples.index(s) for s in samples]

#mysam = "D:/TASSELGBS_Msa/190517aligned_tags_multi.sam"
# dictionary to store alignments
aligndict = dict()
count = 0
locsfound = dict() # set of all alignment locations, for lookup

# lists to store marker data for the current tag
these_mnames = [] # marker names
these_NM = [] # number of mutational steps between tag and reference
lasttagseq = "" # tag sequence

def pad_marker(old_mnames, new_mnames, alignlist):
  '''Function for padding out marker information, if the tag aligned to fewer
  locations than the tags it is being grouped with.
  old_mnames and new_mnames are tuples of alignment locations.
  alignlist is formatted as an element of aligndict.'''
  dummy_NM = 999 # what NM should be if there is not an actual alignment
  m_index = [new_mnames.index(m) for m in old_mnames]
  nalign = len(new_mnames)
  new_alignlist = []
  for aligninfo in alignlist:
    new_NM = [dummy_NM for i in range(nalign)]
    new_CIGAR = ['' for i in range(nalign)] ### Add stuff for MD
    for i in range(len(m_index)):
      mi = m_index[i]
      new_NM[mi] = aligninfo[1][i]
      new_CIGAR[mi] = aligninfo[2][i]
    new_NM = tuple(new_NM)
    new_CIGAR = tuple(new_CIGAR)
    new_alignlist.append((aligninfo[0], new_NM, new_CIGAR))
  return new_alignlist

# subroutine for updating the alignment dictionary
def update_aligndict(these_mnames, these_NM, these_CIGAR, these_MD, lasttagseq):
  these_mnames, these_NM, these_CIGAR, these_MD = \
  zip(*sorted(zip(these_mnames, these_NM, these_CIGAR, these_MD)))
  doupdate = True
  if these_mnames not in aligndict.keys():
    # see if any of these alignment locations have been found before
    mnames_match = [m in locsfound.keys() for m in these_mnames]
    anyoverlap = any(mnames_match)

    if anyoverlap:
      existing_mnames = set([locsfound[these_mnames[mi]] for mi in range(len(these_mnames)) if mnames_match[mi]])
      # if marker name is listed as None, we should ignore all of these alignment locations
      if any([m == None for m in existing_mnames]):
        othernames = [m for m in existing_mnames if m != None]
        for m in these_mnames:
          locsfound[m] = None
        for mset in othernames:
          for m in mset:
            locsfound[m] = None
          del aligndict[mset]
        doupdate = False
      # Merge this set and any matching alignment sets,
      # or delete if there are too many alignment locations
      else:
        allalign = set(these_mnames)
        for mnames in existing_mnames:
          allalign.update(mnames)
        if len(allalign) > maxalign: # too many alignment locations, burn
          for m in allalign:
            locsfound[m] = None
          for mset in existing_mnames:
            del aligndict[mset] 
          doupdate = False
        else: # merge alignment groups
          new_mnames = tuple(sorted(allalign))
          # copy aligndict to avoid overwriting
          aligndict_old = {old_mnames:aligndict[old_mnames] for old_mnames in existing_mnames}
          aligndict[new_mnames] = [] # start new marker group
          for old_mnames in existing_mnames:
            aligndict[new_mnames].extend(pad_marker(old_mnames, new_mnames, aligndict_old[old_mnames]))
            if old_mnames != new_mnames:
              del aligndict[old_mnames]
          these_updated = pad_marker(these_mnames, new_mnames, [(lasttagseq, these_NM, these_CIGAR, these_MD)])
          these_mnames = new_mnames
          these_NM = these_updated[0][1]
          these_CIGAR = these_updated[0][2]
          these_MD = these_updated[0][3]
          for m in these_mnames:
            locsfound[m] = these_mnames 
    else:
      aligndict[these_mnames] = [] # new marker
      for m in these_mnames:
        locsfound[m] = these_mnames
  if doupdate:
    aligndict[these_mnames].append((lasttagseq, these_NM, these_CIGAR, these_MD))

# read the sam file and process to aligndict
print("Processing SAM...")
with open(mysam, mode = "r") as samcon:
  for line in samcon:
    if line[0] == "@": # header lines
      continue
    row = line.split()

    flag = int(row[1])
    if flag & 4 == 4: # no alignment
        continue

    # extract various information about the alignment
    tagseq = row[0][7:]
    chrom = row[2]
    pos = int(row[3])
    cigar = row[5]
    if flag & 16 == 16:
      strand = "bot"
      deletions = sum([int(x[:-1]) for x in re.findall('\d+D', cigar)])
      insertions = sum([int(x[:-1]) for x in re.findall('\d+I', cigar)])
      pos = pos + len(tagseq) - insertions + deletions - 1
    else:
      strand = "top"
    mname = "{}-{:0>{width}}-{}".format(chrom, pos, strand, width=9)
    if row[17].startswith("NM:i:"):
      assert row[18].startswith("MD:Z:")
      NM = int(row[17][5:])
      MD = row[18][5:]
    else:
      assert row[16].startswith("NM:i:")
      assert row[17].startswith("MD:Z:")
      NM = int(row[16][5:])
      MD = row[17][5:]

    secondary = flag & 256 == 256
    if secondary: # add alignments to the list
      these_mnames.append(mname)
      these_NM.append(NM)
      these_CIGAR.append(cigar)
      these_MD.append(MD)
    else: # start a new tag
      # put the last tag into alignment dict if applicable
      if(len(these_mnames) <= maxalign and lasttagseq != ""):
        update_aligndict(these_mnames, these_NM, these_CIGAR, these_MD, lasttagseq)
        count += 1
        if count % 10000 == 0:
          print("{} tags aligning".format(count))
      # reset variables for new tag
      lasttagseq = tagseq
      these_mnames = [mname]
      these_NM = [NM]
      these_CIGAR = [cigar]
      these_MD = [MD]

# write the last marker to dict if necessary
if(len(these_mnames) <= maxalign):
  update_aligndict(these_mnames, these_NM, these_CIGAR, these_MD, lasttagseq)
  count += 1
print("{} tags aligning".format(count))

# filter to polymorphic sites
aligndict = {k:v for k, v in aligndict.items() if len(v) > 1}

## Part two: reading TagTaxaDist, filtering, and export

def extractTTD(tags, ttdfile, sample_index):
  "Extract a given set of tags from the TTD file"
  # set up output matrix for TagTaxaDist
  ttd_mat = [[tag] for tag in tags]
  ntags = len(tags) # number of tags
  # sort tags for binary search
  sorted_tags, tagindex = zip(*sorted(zip(tags, range(ntags))))
  with open(ttdfile, mode = 'r') as mycon:
    header = next(mycon).split()[1:]
    header = [header[i] for i in sample_index]
    for line in mycon:
      row = line.split()
      tag = row[0]
      ti = bisect.bisect_left(sorted_tags, tag)
      if ti == ntags or sorted_tags[ti] != tag:
        continue # skip if this is not a tag we wanted to keep
      ti2 = tagindex[ti]
      ttd_mat[ti2].extend([int(row[i+1]) for i in sample_index])
  if any([len(ttd) == 1 for ttd in ttd_mat]):
    tagsnotfound = [ttd[0] for ttd in ttd_mat if len(ttd) == 1]
    for i in range(min([10, len(tagsnotfound)])):
      print(tagsnotfound[i])
    raise Exception("Tag not found in TagTaxaDist file. Do SAM and TagTaxaDist match?")
  return header, ttd_mat

def keepMarker(depths, min_ind_with_reads):
  "Determine whether to keep a marker, based on missing data rate."
  nind = len(depths[0]) - 1 # first item is tag sequence
  ind_with_reads = [any([d[i + 1] > 0 for d in depths]) for i in range(nind)]
  return sum(ind_with_reads) >= min_ind_with_reads

def makeChunks(all_markers, nchunks):
  '''Get indices for starting and ending points for chunks, putting the
  breakpoints between chromosomes if possible.'''
  # chromosome name for first alignment
  all_chrom = [m[0][:m[0].find('-')] for m in all_markers]
  m_per_chunk = len(all_markers) // nchunks
  chunk_ranges = [[0,0] for i in range(nchunks)]
  for chnk in range(nchunks):
    if chnk > 0:
      # next chunk starts where last one ended
      chunk_ranges[chnk][0] = chunk_ranges[chnk - 1][1]
      # cut chunks short if we have reached end
      if chunk_ranges[chnk][0] == len(all_markers):
        chunk_ranges = chunk_ranges[:chnk]
        break
    if chnk == nchunks - 1:
      # end of all markers
      chunk_ranges[chnk][1] = len(all_markers)
    else:
      # find a good breakpoint
      endm = chunk_ranges[chnk][0] + m_per_chunk
      if endm > len(all_markers): # don't go past end
        endm = len(all_markers)
      thischr = all_chrom[endm - 1]
      for i in range(m_per_chunk // 2): # search to left and right
        if all_chrom[endm - 1 - i] != thischr:
          endm = endm - i
          break
        if endm + i < len(all_markers) and all_chrom[endm + i] != thischr:
          endm = endm + i
          break
      chunk_ranges[chnk][1] = endm
  return chunk_ranges

# divide into chunks and process TTD
print("Sorting alignment locations and processing TagTaxaDist")
all_markers = sorted(aligndict.keys())
# set up markers for each chunk
chunk_ranges = makeChunks(all_markers, nchunks)
for chnk in range(nchunks):
  these_markers = all_markers[chunk_ranges[chnk][0]:chunk_ranges[chnk][1]]
  print("Chunk {}".format(chnk + 1))
  print("{} to {}".format(these_markers[0], these_markers[-1]))
  nm = len(these_markers)
  tag_table = [] # to hold tag info in spreadsheet-like format
  # to index rows that correspond to given markers, same order as these_markers
  table_row_per_marker = [[-1] for i in range(len(these_markers))]
  curr_row = 0
  # loop to build table of tag info
  for mi in range(nm):
    m = these_markers[mi]
    # pad out marker names to maximum
    npad = maxalign - len(m)
    m_exp = list(m) + ['' for i in range(npad)]
    # get tags, NM, and CIGAR
    m_tags = [tup[0] for tup in aligndict[m]]
    m_NM = [list(tup[1]) + ['' for i in range(npad)] for tup in aligndict[m]]
    m_CIGAR = [list(tup[2]) + ['' for i in range(npad)] for tup in aligndict[m]]
    nt = len(m_tags) # number of tags for this marker
    assert len(m_NM) == nt
    # add to table
    tag_table.extend([m_exp + [m_tags[i]] + m_NM[i] + m_CIGAR[i] for i in range(nt)])
    table_row_per_marker[mi] = range(curr_row, curr_row + nt)
    curr_row += nt
  # extract all tag sequences
  these_tags = [tt[maxalign] for tt in tag_table]
  # extract read depth for these tags
  ttdheader, this_ttd = extractTTD(these_tags, myttd, sample_index)
  # determine which markers (rows) to keep
  markers_kept = [keepMarker([this_ttd[ti] for ti in table_row_per_marker[mi]],
                             min_ind_with_reads) for mi in range(nm)]
  rows_kept = [ti for mi in range(nm) for ti in table_row_per_marker[mi] \
               if markers_kept[mi]]
  # subset tag table and depth table
  tag_table = [tag_table[ti] for ti in rows_kept]
  this_ttd = [this_ttd[ti] for ti in rows_kept]
  # export to files
  alignout = "{}_{}_align.csv".format(outbase, chnk + 1)
  ttdout = "{}_{}_depth.csv".format(outbase, chnk + 1)
  with open(alignout, mode = 'w', newline = '') as outcon:
    mywriter = csv.writer(outcon)
    alignheader = ["Alignment {}".format(i + 1) for i in range(maxalign)] + \
    ["Tag sequence"] + ["NM {}".format(i + 1) for i in range(maxalign)] + \
    ["CIGAR {}".format(i + 1) for i in range(maxalign)]
    mywriter.writerow(alignheader)
    [mywriter.writerow(tt) for tt in tag_table]
  with open(ttdout, mode = 'w', newline = '') as outcon:
    mywriter = csv.writer(outcon)
    mywriter.writerow(["Tag sequence"] + ttdheader)
    [mywriter.writerow(tt) for tt in this_ttd]
