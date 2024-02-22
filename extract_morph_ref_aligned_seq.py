#!/usr/bin/env python

import sys

morph_fasta_file = sys.argv[1] # uncompressed .fa
morph_paf_file = sys.argv[2] # with -c --eqx
ref_start_position = int(sys.argv[3])
ref_end_position = int(sys.argv[4]) # exclusive
# fasta to stdout

def find_cigar_query_poses(cigarstr, relative_ref_start, relative_ref_end):
	cigarparts = cigarstr.replace("M", "M\t").replace("I", "I\t").replace("D", "D\t").replace("=", "=\t").replace("X", "X\t").strip().split("\t")
	cigarparts = [(int(part[:-1]), part[-1]) for part in cigarparts]
	result_start = None
	result_end = None
	current_ref_pos = 0
	current_query_pos = 0
	for part in cigarparts:
		if part[1] == 'M' or part[1] == '=' or part[1] == 'X':
			if result_start is None and current_ref_pos + part[0] > relative_ref_start:
				result_start = current_query_pos + (relative_ref_start - current_ref_pos)
			if result_end is None and current_ref_pos + part[0] > relative_ref_end:
				result_end = current_query_pos + (relative_ref_end - current_ref_pos)
			current_query_pos += part[0]
			current_ref_pos += part[0]
		elif part[1] == 'I':
			current_query_pos += part[0]
		elif part[1] == 'D':
			if result_start is None and current_ref_pos + part[0] > relative_ref_start:
				result_start = current_query_pos
			if result_end is None and current_ref_pos + part[0] > relative_ref_end:
				result_end = current_query_pos
			current_ref_pos += part[0]
	return (result_start, result_end)

morph_seqs = {}
current_name = ""
current_seq = ""
with open(morph_fasta_file) as f:
	for l in f:
		if l[0] == '>':
			if len(current_seq) > 0: morph_seqs[current_name] = current_seq
			current_name = l[1:].strip()
			current_seq = ""
		else:
			current_seq += l.strip()
if len(current_seq) > 0: morph_seqs[current_name] = current_seq

with open(morph_paf_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		assert parts[4] == "+" # backwards alignments need extra implementation, just assume it doesn't happen
		if int(parts[7]) > ref_start_position: continue
		if int(parts[8]) < ref_end_position: continue
		tags = parts[13:]
		query_start_pos = None
		query_end_pos = None
		for tag in tags:
			if tag[0:5] == "cg:Z:":
				(query_start_pos, query_end_pos) = find_cigar_query_poses(tag[5:], ref_start_position - int(parts[7]), ref_end_position - int(parts[7]))
		assert parts[0] in morph_seqs
		assert query_start_pos is not None
		assert query_end_pos is not None
		query_start_pos += int(parts[2])
		query_end_pos += int(parts[2])
		if query_end_pos > query_start_pos:
			print(">" + parts[0] + "_" + str(query_start_pos) + "_" + str(query_end_pos))
			print(morph_seqs[parts[0]][query_start_pos:query_end_pos])
