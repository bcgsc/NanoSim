import sys
from time import strftime
import regex
import numpy as np
from sklearn import linear_model
import piecewise_regression


def calc_homopolymer_mis_rate(alnms):
    """
    Calculate mismatch rates in previously identified homopolymers
    :param alnms: Alignment file in MAF format
    :return: Mismatch rate in homopolymers
    """
    err_dict = {"ins": 0, "del": 0, "mis": 0, "match": 0}

    for alnm in alnms:
        ref_seq = alnm[0]
        read_seq = alnm[1]

        err_dict["ins"] += ref_seq.count("-")  # number of insertions
        err_dict["del"] += read_seq.count("-")  # number of deletions

        i = 0
        while i < len(ref_seq):
            if ref_seq[i] != "-" and read_seq[i] != "-":
                if ref_seq[i] != read_seq[i]:  # mismatch
                    err_dict["mis"] += 1
                else:  # match
                    err_dict["match"] += 1
            i += 1

    return err_dict["mis"] / (err_dict["del"] + err_dict["mis"] + err_dict["match"])


def analyze_homopolymers(alnm_file, min_hp_len, prefix):
    """
    Identify homopolymers with length >= min_homopolymer_len in reference genome and calculate corresponding length
    of homopolymer in reads from alignment file (MAF).
    :param alnm_file: Alignment file in MAF format
    :param min_hp_length: Minimum homopolymer length in reference genome to analyze
    :param prefix: Prefix of output files
    :return:
        Dictionary of homopolymer lengths per base {"base": {ref_len: [read_len_1, read_len_2, ...], ...}, ...}
        Array of ref and read homopolymer sequences in MAF format [[ref_hp_seq, read_hp_seq, base], ...]
    """
    hp_lengths = []  # [[ref_pos, type, ref_homopolymer_len, read_homopolymer_len], ... ]
    hp_lengths_per_base = {"AT": {}, "CG": {}}  # {"AT": {"ref_hp_len": [read_hp_len_1, read_hp_len_2, ...], ...}, ...}
    hp_alnms = []  # [[ref_seq, read_seq, base], ...]

    i = 1
    with open(alnm_file, 'r') as f:
        for line in f:
            ref = line
            ref_info = ref.split()
            ref_name = ref_info[1]
            ref_start = int(ref_info[2])
            ref_seq = ref_info[6]

            read = next(f)
            read_info = read.split()
            read_seq = read_info[6]

            ref_coords = []  # [[base, start, end], ...]: Base, start pos, and end pos of a homopolymer in the reference
            ref_seq_without_dash = ref_seq.replace("-", "")
            pattern = "A{" + min_hp_len + ",}|C{" + min_hp_len + ",}|G{" + min_hp_len + ",}|T{" + min_hp_len + ",}"
            for match in regex.finditer(pattern, ref_seq_without_dash):
                ref_coords.append([match.group()[0], ref_start + match.start(), ref_start + match.end()])

            aligned_coords = []  # [[start, end], ...]
            pattern = "(-*A-*){" + min_hp_len + ",}|(-*C-*){" + min_hp_len + ",}|(-*G-*){" + min_hp_len + ",}|(-*T-*){" + min_hp_len + ",}"
            for match in regex.finditer(pattern, ref_seq):
                aligned_coords.append([match.start(), match.end()])

            non_hp_ref_seq = ""
            non_hp_read_seq = ""
            last_pos = 0
            for item1, item2 in zip(ref_coords, aligned_coords):
                base = item1[0]
                ref_start = item1[1]
                ref_end = item1[2]
                start = item2[0]
                end = item2[1]

                ref_homopolymer_len = len(ref_seq[start: end].replace("-", ""))
                read_homopolymer = read_seq[start: end].replace("-", "")

                # Find max length matching homopolymer in read, allowing for one mismatch
                read_homopolymer_len = 0
                for match in regex.finditer("(" + base * 2 + "+){s<=1}", read_homopolymer):
                # for match in regex.finditer("((" + base * 2 + "+){s<=1})+", read_homopolymer): # might be a better regex pattern to use
                    if match:
                        match_seq = match.group()
                        match_len = len(match.group())

                        # Accounting for mismatches at beginning and end of matches (i.e. should not be part of match)
                        if match_seq[0] != base:
                            match_seq = match_seq[1:]
                            match_len -= 1

                        if match_seq[-1] != base:
                            match_seq = match_seq[:-1]
                            match_len -= 1

                        if read_homopolymer_len < match_len:
                            read_homopolymer_len = match_len

                non_hp_ref_seq += ref_seq[last_pos: start]
                non_hp_read_seq += read_seq[last_pos: start]
                last_pos = end

                hp_alnms.append([ref_seq[start: end], read_seq[start: end], base])
                hp_lengths.append([ref_name + ":" + str(ref_start + 1) + "-" + str(ref_end + 1), str(base),
                                   str(ref_homopolymer_len), str(read_homopolymer_len)])

                key = "AT" if base in ["A", "T"] else "CG"
                if ref_homopolymer_len not in hp_lengths_per_base[key].keys():
                    hp_lengths_per_base[key][ref_homopolymer_len] = []
                hp_lengths_per_base[key][ref_homopolymer_len].append(read_homopolymer_len)

            # non_hp_alnms.append([non_hp_ref_seq, non_hp_read_seq])
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished alignment >> " + str(i) + "\r")
            sys.stdout.flush()
            i += 1

    # Counting number of homopolymers at each position
    hp_lengths = np.unique(hp_lengths, axis=0, return_counts=True)
    hp_lengths_by_pos = []
    for length, count in zip(hp_lengths[0], hp_lengths[1]):
        hp_lengths_by_pos.append(np.append(length, count))

    # Writing to file
    out_lengths = open(prefix + "_hp_lengths.tsv", "w")
    out_lengths.write("Chrom:Ref pos\tType\tRef length\tRead length\tCount\n")
    for info in hp_lengths_by_pos:
        out_lengths.write(str(info[0]) + "\t" + str(info[1]) + "\t" + str(info[2]) + "\t" + str(info[3]) + "\t" +
                          str(info[4]) + "\n")
    out_lengths.close()
    return hp_lengths_per_base, hp_alnms


def fit_piecewise(hp_lengths_per_base):
    header = ""
    estimates_per_base = {"AT": "", "CG": ""}
    for base in hp_lengths_per_base.keys():
        xx = list(hp_lengths_per_base[base].keys())  # ref hp lengths
        yy = list(np.mean(hp_read_lengths, dtype=np.float64) for hp_read_lengths in hp_lengths_per_base[base].values())
        pw_fit = piecewise_regression.Fit(xx, yy, n_breakpoints=1)  # n_breakpoints set based on previous analysis; may want to optimize this later
        pw_estimates = pw_fit.get_results()["estimates"]

        if header == "":
            header = "\t".join(list(pw_estimates.keys()))

        estimates = []
        for param in pw_estimates.keys():
            estimates.append(str(pw_estimates[param]["estimate"]))
        estimates_per_base[base] = "\t".join(estimates)

    return header, estimates_per_base


# Adapted from https://github.com/chasmani/piecewise-regression/pull/10
def predict_piecewise(hp_length, piecewise_estimates):
    piecewise_params = piecewise_estimates.keys()
    intercept_hat = float(piecewise_estimates["const"])
    alpha_hat = float(piecewise_estimates["alpha1"])
    xx = hp_length
    beta_hats = []
    breakpoints = []
    for param in piecewise_params:
        if "breakpoint" in param:
            breakpoints.append(piecewise_estimates[param])
        elif "beta" in param:
            beta_hats.append(piecewise_estimates[param])

    # Build the fit plot segment by segment. Betas are defined as
    # difference in gradient from previous section
    yy_pred = intercept_hat + alpha_hat * xx
    for bp_count in range(len(breakpoints)):
        yy_pred += beta_hats[bp_count] * \
                   np.maximum(xx - breakpoints[bp_count], 0)
    return yy_pred


def fit_lr(hp_lengths_per_base):
    header = "intercept\tslope"
    estimates_per_base = {"AT": "", "CG": ""}
    for base in hp_lengths_per_base.keys():
        x = np.fromiter(hp_lengths_per_base[base].keys(), dtype=int).reshape(-1, 1)  # ref hp lengths
        y = np.fromiter([np.std(hp_read_lengths, dtype=np.float64) for hp_read_lengths in hp_lengths_per_base[base].values()],
                        dtype=np.float64)
        lr = linear_model.LinearRegression(fit_intercept=False)
        lr.fit(x, y)

        estimates_per_base[base] = str(lr.intercept_) + "\t" + str(lr.coef_[0])

    return header, estimates_per_base


def predict_lr(hp_length, lr_estimates):
    intercept = lr_estimates["intercept"]
    slope = lr_estimates["slope"]
    x = hp_length
    y_pred = intercept + slope * x
    return y_pred


def model_homopolymer_lengths(alnm_file, min_hp_len, prefix):
    """
    Estimates model parameters based on homopolymer lengths identified in alignment file (MAF)
    :param alnm_file: Alignment file in MAF format
    :param min_hp_length: Minimum homopolymer length in reference genome to analyze
    :param prefix: Prefix of output files
    """

    # Identify homopolymers in alignment file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parsing alignment file for homopolymers\n")
    sys.stdout.flush()
    hp_lengths_per_base, hp_alnms = analyze_homopolymers(alnm_file, min_hp_len, prefix)

    # Calculate mismatch rates in homopolymers
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Calculating mismatch rates in homopolymers\n")
    sys.stdout.flush()
    mismatch_rate = calc_homopolymer_mis_rate(hp_alnms)

    # Estimate model parameters and write to file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Estimating model parameters\n")
    sys.stdout.flush()
    pw_header, pw_estimates = fit_piecewise(hp_lengths_per_base)
    lr_header, lr_estimates = fit_lr(hp_lengths_per_base)

    out_model_params = open(prefix + "_hp_lengths_model_parameters.tsv", "w")
    header = "base\t" + pw_header + "\t" + lr_header + "\n"
    all_estimates = {base: pw_estimates[base] + "\t" + lr_estimates[base] for base in pw_estimates.keys()}
    out_model_params.write("#Homopolymer mismatch rate: " + str(mismatch_rate) + "\n")
    out_model_params.write(header)
    for base in all_estimates:
        out_model_params.write(base + "\t" + all_estimates[base] + "\n")
    out_model_params.close()


def get_nd_par(ref_hp_length, pw_estimates, lr_estimates):
    """
    Samples homopolymer lengths from normal distribution based on reference homopolymer length and previously estimated model parameters.
    :param alnm_file: Alignment file in MAF format
    :param min_hp_length: Minimum homopolymer length in reference genome to analyze
    :param prefix: Prefix of output files
    :return: Normal distribution parameters for each base
    """
    parameters = {"A": {}, "T": {}, "C": {}, "G": {}}
    for base in pw_estimates.keys():
        parameters[base]["mu"] = predict_piecewise(ref_hp_length, pw_estimates[base])
        parameters[base]["sigma"] = predict_lr(ref_hp_length, lr_estimates[base])

    return parameters["AT"]["mu"], parameters["AT"]["sigma"], parameters["AT"]["mu"], parameters["AT"]["sigma"], \
           parameters["CG"]["mu"], parameters["CG"]["sigma"], parameters["CG"]["mu"], parameters["CG"]["sigma"]

