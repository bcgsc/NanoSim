import sys
from time import strftime
import re
import pysam
from scipy.stats import rv_discrete, lognorm
import numpy as np


class trunc_lognorm_gen(rv_discrete):
    def _cdf(self, x, s, scale):
        lncdf_x = lognorm.cdf(x, s, scale=scale)
        lncdf_a = lognorm.cdf(self.a, s, scale=scale)
        lncdf_b = lognorm.cdf(self.b, s, scale=scale)
        return (lncdf_x - lncdf_a) / (lncdf_b - lncdf_a)

    def _ppf(self, q, s, scale):
        lncdf_a = lognorm.cdf(self.a, s, scale=scale)
        lncdf_b = lognorm.cdf(self.b, s, scale=scale)
        ln_q = q * (lncdf_b - lncdf_a) + lncdf_a
        return lognorm.ppf(ln_q, s, scale=scale)


def convert_cs(cs_string):
    cs_arr = []
    for item in re.findall('(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)', cs_string):
        op = item[0]
        if op == ":":
            size = int(item[1:])
        elif op == "+":
            size = len(item) - 1
        elif op == "*":
            size = 1
        else:
            continue
        cs_arr += op * size
    return cs_arr


def analyze_ht_quals(alnm):
    ht_quals = []

    head = alnm.query_sequence[0: alnm.query_alignment_start]
    tail = alnm.query_sequence[alnm.query_alignment_end: len(alnm.query_sequence)]
    if head:
        head_quals = alnm.query_qualities[0: alnm.query_alignment_start]
        ht_quals += head_quals.tolist()

    if tail:
        tail_quals = alnm.query_qualities[alnm.query_alignment_end: len(alnm.query_sequence)]
        ht_quals += tail_quals.tolist()

    return ht_quals


def analyze_aligned_base_qualities(primary_alnm_file):
    """
    Based on cs tag for each alignment, extract qualities for each error type (mis/ins), match, and head/tail base.
    :param primary_alnm_file: Primary alignment file
    :return Dictionary of qualities per type {"type": [base_qualities, ...], ...}
    """
    alignments = pysam.AlignmentFile(primary_alnm_file, "rb")
    quals_per_type = {"mis": [], "ins": [], "match": [], "ht": [], "unmapped": []}

    cs_dict = {":": "match", "+": "ins", "*": "mis"}
    for alnm in alignments:
        if alnm.is_secondary:
            print(f"Skipping secondary alignment: {alnm.query_name}", file=sys.stderr)
            continue
        quals = alnm.query_alignment_qualities.tolist()

        read = alnm.query_alignment_sequence
        cs_arr = convert_cs(alnm.get_tag("cs"))  # expand cs tag into individual characters
        for i, base in enumerate(read):
            err = cs_dict[cs_arr[i]]
            quals_per_type[err].append(quals[i])

        quals_per_type["ht"] += analyze_ht_quals(alnm)

    return quals_per_type


def fit_lognorm(quals_per_type, prefix):
    out_model_params = open(prefix + "_base_qualities_model_parameters.tsv", "w")
    out_model_params.write("type\tsd\tloc\tmu\n")
    for type in quals_per_type.keys():
        sample_size = 500000
        if len(quals_per_type[type]) > sample_size:
            quals_sample = np.random.choice(quals_per_type[type],
                                            sample_size)  # sampled because fitting takes a long time
        else:
            quals_sample = quals_per_type[type]

        sd, loc, scale = lognorm.fit(quals_sample, floc=0)  # fix loc to 0 for now, might want to change later
        mu = np.log(scale)
        out_model_params.write(type + "\t" + str(sd) + "\t" + str(loc) + "\t" + str(mu) + "\n")
    out_model_params.close()


def model_base_qualities(alnm_file, prefix, unmapped_base_quals):
    """
    Estimates model parameters based on base qualities identified in alignment file. Model parameters are calculated
    separately for matches and per error type (mis, ins).
    :param alnm_file: Alignment file
    :param prefix: Prefix of output files
    :param unmapped_base_quals: List of base qualities from unmapped bases
    """
    # Analyze base qualities
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parsing alignment file for base qualities "
                     "relative to matches and each error type\n")
    sys.stdout.flush()
    quals_per_type = analyze_aligned_base_qualities(alnm_file)
    quals_per_type["unmapped"] = unmapped_base_quals

    # Estimate model parameters and write to file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Estimating model parameters\n")
    sys.stdout.flush()
    fit_lognorm(quals_per_type, prefix)


def predict_base_qualities(sd, loc, scale, n):
    """
    Samples base qualities log normal distribution based on previously estimated model parameters.
    :param sd: Standad deviation
    :param loc: Location
    :param scale: Scale (exp(mu))
    :param n: Sample size
    :return Array of sampled base qualities
    """
    truncln_gen = trunc_lognorm_gen(name="truncln", a=1, b=93, shapes="s, scale")
    return truncln_gen.rvs(s=sd, loc=loc, scale=scale, size=n).tolist()

