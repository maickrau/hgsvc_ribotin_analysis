
exit # not meant to run automatically, run commands manually one at a time instead

# merge to single files, copy consensus
cat hgsvc_ribotin/*/consensus.fa | sed 's/>/\n>/g' | grep -vP '^$' > merged_consensuses.fa
cat hgsvc_ribotin/*/morphs.fa > merged_morphs.fa
cp ribotin/template_seqs/rDNA_one_unit.fasta .

# align all-vs-all and to ref and to consensuses
minimap2 -x asm20 -X -c --eqx -t 24 merged_morphs.fa merged_morphs.fa > alns_ava.paf
minimap2 -x asm20 -c --eqx -t 24 rDNA_one_unit.fasta merged_morphs.fa > alns_to_ref.paf
minimap2 -x asm20 -c --eqx -t 24 merged_consensuses.fa merged_morphs.fa > alns_to_consensus.paf

# morph most similar to reference
awk '$4-$3>$2*0.9' < alns_to_ref.paf | cut -f 6,7,13 | sed 's/NM:i://g' | awk '{error_rate = $3/$2; print error_rate "\t" $1}' | sort -nr | awk '{rate_per_ref[$2] = $1;}END{for (ref in rate_per_ref){print ref "\t" rate_per_ref[ref];}}' | less
# morph most similar to consensuses in the same sample & tangle
awk 'substr($1,1,16)==substr($6,1,16)' | awk '$4-$3>$2*0.9' < alns_to_consensus.paf | cut -f 6,7,13 | sed 's/NM:i://g' | awk '{error_rate = $3/$2; print error_rate "\t" $1}' | sort -nr | awk '{rate_per_ref[$2] = $1;}END{for (ref in rate_per_ref){print ref "\t" rate_per_ref[ref];}}' | less
# some of them are reasonableor even identical, probably because of the verkko mode per-chromosome tangle separation

# extract 5.8s, 18s, 28s sequences per morph
./extract_morph_ref_aligned_seq.py merged_morphs.fa alns_to_ref.paf 6597 6753 > 5.8s.fa
./extract_morph_ref_aligned_seq.py merged_morphs.fa alns_to_ref.paf 3658 5526 > 18s.fa
./extract_morph_ref_aligned_seq.py merged_morphs.fa alns_to_ref.paf 7921 12971 > 28s.fa

# major 5.8s allele, count 2031/2079. 156bp long
grep -v '>' < 5.8s.fa | sort | uniq -c | sort -nr | head -n 1 | sed 's/^ \+[0-9]\+ //g' | less
# major 18s allele, count 1862/2079. 1868bp long
grep -v '>' < 18s.fa | sort | uniq -c | sort -nr | head -n 1 | sed 's/^ \+[0-9]\+ //g' | less
# major 28s allele does not exist. most abundant allele has count 12/2079, 5065bp long
grep -v '>' < 28s.fa | sort | uniq -c | sort -nr | head -n 1 | sed 's/^ \+[0-9]\+ //g' | less

# could the lack of a 28s major allele be simply due to error rate?
# assume the ground truth is they are all identical, what kind of error rate would cause this result?
# sequence     # in major allele     length     correct/sequence     correct/bp       error % /bp
# 18s          1862                  1868bp     0.895622896          0.999940989      0.0059011 %
# 5.8s         2031                  156bp      0.976911977          0.999850276      0.0149724 %
# 28s          12                    5070bp     0.005772006          0.998983804      0.1016196 %
# unclear but hints towards more genuine variation in 28s
./count_microsatellite_seq.py < 5.8s.fa | awk '{sum += $1; count += 1;}END{print sum/count;}'
# 0.576937
./count_microsatellite_seq.py < 18s.fa | awk '{sum += $1; count += 1;}END{print sum/count;}'
# 0.568487
./count_microsatellite_seq.py < 28s.fa | awk '{sum += $1; count += 1;}END{print sum/count;}'
# 0.720986
# 28s has a higher microsatellite fraction which makes it plausible that the error rate really is higher
# but! probably not that much higher. ~1.25x as much microsatellites, so naively error rate should be 1.25x ??
# and 18s has no more microsatellites than 5.8s, so it probably has more real variation than 5.8s
# assuming "truth" error rate is proportional to msat fraction and 5.8s has zero real variation
#  then "real" error rate in 28s would be ~0.007% and ~70% of alleles would survive error-free
# so 28s probably just has more variation??

grep -v '>' < 28s.fa | awk '{print length($0);}' | sort -n | uniq -c | awk '{sum += $1; print $1 "\t" $2 "\t" sum;}' | less
# 28s median length is 5070bp, but only 76/2079 have that length

./get_morphs_with_seq.py GACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATCGACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGCTT < merged_morphs.fa > morphs_with_major_5.8s.txt
./get_morphs_with_seq.py ACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCTTCGCGGGGGGGATGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCTCTCCGGCCCCGGCCGGGGGGCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCCCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGCCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACCGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA < merged_morphs.fa > morphs_with_major_18s.txt

grep '>' < merged_morphs.fa | awk '{print substr($1,2,7);}' | sort | uniq | wc -l
awk '{print substr($1,1,7);}' < morphs_with_major_5.8s.txt | sort | uniq > samples_with_major_5.8s.txt
awk '{print substr($1,1,7);}' < morphs_with_major_18s.txt | sort | uniq > samples_with_major_18s.txt

cat morphs_with_major_5.8s.txt morphs_with_major_18s.txt | sort | uniq -c | grep -P '^\s*2\s' | sed 's/^ \+[0-9]\+ //g' > morphs_with_major_5.8s_and_18s.txt
awk '{print substr($1,1,7);}' < morphs_with_major_5.8s_and_18s.txt | sort | uniq > samples_with_morph_major_5.8s_and_18s.txt
cut -d '_' -f 1,4 < morphs_with_major_5.8s_and_18s.txt | sed 's/coverage//g' | tr '_' '\t' | awk '{count[$1] += $2;}END{for (sample in count) { print sample "\t" count[sample];}}' > total_coverage_major_5.8s_and_18s_per_sample.txt

wc -l samples_with_m*
#  62 samples_with_major_18s.txt
#  62 samples_with_major_5.8s.txt
#  62 samples_with_morph_major_5.8s_and_18s.txt
# all samples have a morph with the major 5.8s and 18s allele !

# most covered individual morph with major 5.8s and 18s alleles per sample
cut -f 1,4 -d '_' < morphs_with_major_5.8s_and_18s.txt | sed 's/coverage//g' | tr '_' '\t' | awk '{print $2 "\t" $1}' | sort -n | awk '{coverage[$2] = $1;}END{for (sample in coverage) { print sample "\t" coverage[sample];}}' | less

# sample with most coverage in major 5.8s & 18s morphs
awk '{print $2 "\t" $1;}' < total_coverage_major_5.8s_and_18s_per_sample.txt | sort -nr | head -n 1
# 10042 NA24385
# sample with least coverage in major 5.8s & 18s morphs
awk '{print $2 "\t" $1;}' < total_coverage_major_5.8s_and_18s_per_sample.txt | sort -nr | tail -n 1
# 1181 HG02282
