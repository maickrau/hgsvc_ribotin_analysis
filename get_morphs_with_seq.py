#!/usr/bin/env python

import sys

sequence = sys.argv[1]
# fasta from stdin
# list of names to stdout

current_name = ""
current_seq = ""
for l in sys.stdin:
	if l[0] == '>':
		if len(current_seq) > 0:
			if sequence in current_seq: print(current_name)
		current_name = l[1:].strip()
		current_seq = ""
	else:
		current_seq += l.strip()
if len(current_seq) > 0:
	if sequence in current_seq: print(current_name)


