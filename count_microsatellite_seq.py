#!/usr/bin/env python

import sys
import re

# fasta from stdin
# fraction msat to stdout

def count_microsatellite_fraction(s):
	hpc_sequence = re.sub("([ATCG])\\1+", "\\1", s)
	hpc_poses = [0]
	for i in range(1, len(s)):
		if s[i-1] != s[i]:
			hpc_poses.append(i)
	assert len(hpc_poses) == len(hpc_sequence)
	hpc_poses.append(len(s))
	is_microsatellite = []
	for i in range(0, len(hpc_sequence)):
		is_microsatellite.append(False)
	for motif_length in range(2, 7):
		for i in range(2*motif_length, len(hpc_sequence)):
			if hpc_sequence[i-2*motif_length:i-motif_length] == hpc_sequence[i-motif_length:i]:
				for j in range(i-2*motif_length, i+1):
					is_microsatellite[j] = True
	msat_count = 0
	for i in range(0, len(hpc_sequence)):
		if is_microsatellite[i]:
			msat_count += hpc_poses[i+1] - hpc_poses[i]
	return float(msat_count)/float(len(s))

current_name = ""
current_seq = ""
for l in sys.stdin:
	if l[0] == '>':
		if len(current_seq) > 0: print(count_microsatellite_fraction(current_seq))
		current_name = l[1:].strip()
		current_seq = ""
	else:
		current_seq += l.strip()
if len(current_seq) > 0: print(count_microsatellite_fraction(current_seq))
