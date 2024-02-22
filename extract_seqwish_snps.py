#!/usr/bin/env python

import sys

# seqwish gfa from stdin
# phylo file to stdout
# name mapping to namemap.txt

def revcomp(s):
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return "".join(comp[c] for c in s[::-1])

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(n1, n2):
	fwstr = n1 + n2
	bwstr = revnode(n2) + revnode(n1)
	if bwstr < fwstr: return (revnode(n2), revnode(n1))
	return (n1, n2)

def getone(s):
	for c in s: return c

def is_valid_snp_bubble(bubble, node_coverage_in_paths, node_repeats_in_path, nodeseqs, edges, paths):
	if node_coverage_in_paths[bubble[0][1:]] != len(paths): return False
	if node_coverage_in_paths[bubble[1][1:]] != len(paths): return False
	if bubble[0][1:] in node_repeats_in_path: return False
	if bubble[1][1:] in node_repeats_in_path: return False
	allele_coverage_sum = 0
	for edge in edges[bubble[0]]:
		if edge[1:] in node_repeats_in_path: return False
		allele_coverage_sum += node_coverage_in_paths[edge[1:]]
		if len(nodeseqs[edge[1:]]) != 1: return False
	if allele_coverage_sum != len(paths): return False
	return True

# only simple bubbles, possibly multiallelic
def find_bubble(start, edges):
	if start not in edges: return None
	if len(edges[start]) < 2: return None
	otherside = None
	for node in edges[start]:
		if len(edges[revnode(node)]) != 1: return None
		if node not in edges: return None
		if len(edges[node]) != 1: return None
		if otherside is None:
			otherside = getone(edges[node])
		else:
			if getone(edges[node]) != otherside: return None
	if len(edges[revnode(otherside)]) != len(edges[start]): return None
	return (start, otherside)

def format_name(name):
	parts = name.split("_")
	sample = parts[0]
	tangle = parts[1].replace("tangle", "")
	index = parts[2].replace("morphconsensus", "")
	if len(index) == 1: index = "0" + index
	result = sample + tangle + index
	assert len(result) <= 10
	while len(result) < 10: result += " "
	return result

edges = {}
nodeseqs = {}
paths = []

for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "S":
		nodeseqs[parts[1]] = parts[2]
	elif parts[0] == "L":
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		if fromnode not in edges: edges[fromnode] = set()
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[fromnode].add(tonode)
		edges[revnode(tonode)].add(revnode(fromnode))
	elif parts[0] == "P":
		pathname = parts[1]
		path = [(">" if c[-1] == "+" else "<") + c[:-1] for c in parts[2].split(",")]
		paths.append((pathname, path))

node_coverage_in_paths = {}
node_repeats_in_path = set()

for pathpair in paths:
	coverage = {}
	for node in pathpair[1]:
		if node[1:] not in coverage: coverage[node[1:]] = 0
		coverage[node[1:]] += 1
	for node in coverage:
		if node not in node_coverage_in_paths: node_coverage_in_paths[node] = 0
		node_coverage_in_paths[node] += 1
		if coverage[node] != 1: node_repeats_in_path.add(node)

bubble_indices = {}
bubble_count = 0
for edge in edges:
	bubble = find_bubble(edge, edges)
	if not bubble: continue
	if bubble != canon(bubble[0], bubble[1]): continue
	if not is_valid_snp_bubble(bubble, node_coverage_in_paths, node_repeats_in_path, nodeseqs, edges, paths): continue
	for node in edges[edge]:
		bubble_indices[node[1:]] = bubble_count
	bubble_count += 1

num_rows = len(paths)
num_columns = bubble_count

print(str(num_rows) + " " + str(num_columns))

namemap = open("namemap.txt", "w")

for path in paths:
	alleles = []
	for i in range(0, bubble_count): alleles.append(None)
	for node in path[1]:
		if node[1:] not in bubble_indices: continue
		index = bubble_indices[node[1:]]
		assert alleles[index] is None
		seq = nodeseqs[node[1:]]
		if node[0] == "<": seq = revcomp(seq)
		assert len(seq) == 1
		alleles[index] = seq
	for i in range(0, bubble_count):
		assert alleles[i] is not None
		assert len(alleles[i]) == 1
	name = format_name(path[0])
	namemap.write(path[0] + "\t" + name + "\n")
	assert len("".join(alleles)) == bubble_count
	print(name + " " + "".join(alleles))
