#!/usr/bin/env python

import sys

nodesfile = sys.argv[1]
# gfa from stdin
# gfa to stdout

nodes = set()
with open(nodesfile) as f:
	for l in f:
		nodes.add(l.strip())

for l in sys.stdin:
	parts = l.strip().split("\t")
	if parts[0] == "S" and parts[1] in nodes: print(l.strip())
	if parts[0] == "L" and parts[1] in nodes and parts[3] in nodes: print(l.strip())
	if parts[0] == "P" and parts[2].split(",")[0][:-1] in nodes: print(l.strip()) # kinda wrong, only checks first node.
