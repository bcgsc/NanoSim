#! /usr/bin/env python

import sys
import argparse
import collections


def parse_paf(line):
    out = dict()
    fields = line.rstrip().split()
    out["query_name"] = fields[0]
    out["query_length"] = int(fields[1])
    out["target_name"] = fields[5]
    out["target_start"] = int(fields[7])
    out["num_matches"] = int(fields[9])
    out["alignment_block_length"] = int(fields[10])
    return out

# Process the alignment records for a single read and fractionally assign
# that read to the set of transcripts the read is "compatible" with.
# Compatibility is defined by the length of an alignment relative
# to the longest alignment for the read.
def get_compatibility(records):

    full_length_min_distance = 20
    min_read_length = 0

    is_full_length = lambda p : p < full_length_min_distance

    # Determine best match
    read_length = records[0]["query_length"]
    best_match_align_len = 0
    best_num_matches = 0
    best_is_full_length = False

    for r in records:
        fl = is_full_length(r["target_start"])
        #sys.stderr.write("%s match: %d abl: %d ts: %d full: %d\n" % (r["query_name"], r["num_matches"], r["alignment_block_length"], r["target_start"], b))
        if r["num_matches"] > best_num_matches or (r["num_matches"] == best_num_matches and fl):
            best_match_align_len = r["alignment_block_length"]
            best_num_matches = r["num_matches"]
            best_is_full_length = fl

    fraction_aligned = best_match_align_len / float(read_length)
    if fraction_aligned < 0.5 or read_length < min_read_length:
        #sys.stderr.write("Skip %s %d %.2f\n" % (r["query_name"], read_length, fraction_aligned))
        return

    # All records within threshold of the best score are considered to be compatible with the read
    threshold = 0.95

    def is_equivalent_hit(x):
        f = float(x["num_matches"]) / best_num_matches
        l = is_full_length(x["target_start"])
        return f > threshold and l == best_is_full_length

    # Count equivalent hits
    num_hits = 0
    for r in records:
        if is_equivalent_hit(r):
            num_hits += 1

    for r in records:
        if is_equivalent_hit(r):
            transcript_compatibility[r["query_name"]].append((r["target_name"], 1.0 / num_hits))

# Calculate the abundance of the transcript set based on read-transcript compatibilities
def calculate_abundance(compatibility):
    abundance = collections.defaultdict(float)
    total = 0
    for read in compatibility:

        if args.verbose > 1:
            sys.stderr.write("[compatibility] %s: %s\n" % (read, compatibility[read]))

        for t in compatibility[read]:
            abundance[t[0]] += t[1]
            total += t[1]

    for transcript in abundance:
        abundance[transcript] = abundance[transcript] / total
        if args.verbose > 0:
            sys.stderr.write("[abundance] %s: %.4f\n" % (transcript, abundance[transcript]))
    return abundance

# Update read-transcript compatibility based on transcript abundances
def update_compatibility(compatibility, abundance):
    for read in compatibility:

        ids = list()
        total = 0
        for t in compatibility[read]:
            total += abundance[t[0]]
            ids.append(t[0])

        compatibility[read] = list()
        for i in ids:
            compatibility[read].append((i, abundance[i] / total))
#
# Read arguments
#
parser = argparse.ArgumentParser( description='Calculate transcript abundance from minimap2 alignment to transcripts')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-c', '--compatibility', type=str, required=False, default="")
parser.add_argument('-n', '--em-iterations', type=int, default=10)
parser.add_argument('-v', '--verbose', type=int, default=0)
args = parser.parse_args()

#
# Read the input PAF file and calculate the initial read-transcript compatibility
#
transcript_compatibility = collections.defaultdict(list)
prev_read_name = ""
curr_records = list()

fh = open(args.input)
for line in fh:
    alignment = parse_paf(line)
    if alignment["query_name"] != prev_read_name and len(curr_records) > 0:
        get_compatibility(curr_records)
        curr_records = list()
    prev_read_name = alignment["query_name"]
    curr_records.append(alignment)

# Process last batch
get_compatibility(curr_records)

#
# Run EM to calculate abundance and update read-transcript compatibility
#
for i in range(0, args.em_iterations):
    sys.stderr.write("EM iteration %d\n" % (i))

    # Calculate abundance from compatibility assignments
    abundance = calculate_abundance(transcript_compatibility)

    # Update compatibility assignments
    update_compatibility(transcript_compatibility, abundance)

# Write results as a TSV file
total_reads = len(transcript_compatibility)
sys.stderr.write("Parsed alignments for %s reads\n" % (total_reads))
print("target_id\test_counts\ttpm")
for transcript_id in abundance:
    a = abundance[transcript_id]
    tc = a * total_reads
    tpm = a * 1000000
    print("%s\t%.4lf\t%.4lf" % (transcript_id, tc, tpm))

# Write read-transcript assignments
if args.compatibility != "":
    compatibility_writer = open(args.compatibility, "w")
    for read in transcript_compatibility:

        num_compat = len(transcript_compatibility[read])
        ids = list()
        probs = list()
        tpms = list()
        for t in transcript_compatibility[read]:

            if t[1] > 0.1:
                ids.append(t[0])
                probs.append(t[1])
                tpms.append(abundance[t[0]] * 1000000)

        compatibility_writer.write("%s\t%d\t%s\t%s\t%s\n" % (read, num_compat, ",".join(ids), ",".join([str(x) for x in probs]), ",".join([str(x) for x in tpms])))
