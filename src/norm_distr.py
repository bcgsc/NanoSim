def get_nd_par(ref_len, basecaller):
    if basecaller == "albacore":
        return albacore_nd_par(ref_len)
    elif basecaller == "guppy":
        return guppy_nd_par(ref_len)
    elif basecaller == "guppy-flipflop":
        return guppy_ff_nd_par(ref_len)


def seg_par(ref_len, changepoint):
    if ref_len > changepoint:
        xe2 = 0
        xe3 = ref_len - changepoint
    else:
        xe2 = ref_len - changepoint
        xe3 = 0

    xe1 = ref_len - changepoint
    return xe1, xe2, xe3


def albacore_nd_par(ref_len):
    at_xe1, at_xe2, at_xe3 = seg_par(ref_len, 12)
    at_mu = 9.32968110 + 0.75215056 * at_xe1 + 0.01234263 * at_xe2**2 - 0.02699184 * at_xe3**2
    at_sigma = 0.2507 * ref_len - 0.1510

    cg_xe1, cg_xe2, cg_xe3 = seg_par(ref_len, 7)
    cg_mu = 4.76783156 + 0.21751512 * cg_xe1 - 0.10414613 * cg_xe2**2 - 0.01647626 * cg_xe3**2
    cg_sigma = 0.1109 * ref_len + 0.8959

    return at_mu, at_sigma, cg_mu, cg_sigma


def guppy_nd_par(ref_len):
    at_xe1, at_xe2, at_xe3 = seg_par(ref_len, 16)
    at_mu = 12.13283713 + 0.71843395 * at_xe1 + 0.00127124 * at_xe2 ** 2 - 0.01429113 * at_xe3 ** 2
    at_sigma = 0.2138 * ref_len + 0.4799

    cg_xe1, cg_xe2, cg_xe3 = seg_par(ref_len, 12)
    cg_mu = 6.13526513 + 0.21219760 * cg_xe1 - 0.01249273 * cg_xe2 ** 2 - 0.04870821 * cg_xe3 ** 2
    cg_sigma = 0.3021 * ref_len - 0.06803

    return at_mu, at_sigma, cg_mu, cg_sigma


def guppy_ff_nd_par(ref_len):
    at_xe1, at_xe2, at_xe3 = seg_par(ref_len, 12)
    at_mu = 9.65577227 + 0.92524095 * at_xe1 + 0.02814258 * at_xe2 ** 2 - 0.01666699 * at_xe3 ** 2
    at_sigma = 0.2013 * ref_len + 0.2159

    cg_xe1, cg_xe2, cg_xe3 = seg_par(ref_len, 14)
    cg_mu = 5.5402163422 - 0.0163232962 * cg_xe1 - 0.0205230566 * cg_xe2 ** 2 + 0.0009633945 * cg_xe3 ** 2
    cg_sigma = 0.1514 * ref_len + 0.5611

    return at_mu, at_sigma, cg_mu, cg_sigma


def get_hpmis_rate(basecaller):
    if basecaller == "albacore":
        return 0.02204
    elif basecaller == "guppy":
        return 0.02166
    elif basecaller == "guppy-flipflop":
        return 0.02215
