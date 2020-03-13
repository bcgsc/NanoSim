def get_nd_par(ref_len, read_type, basecaller):
    if read_type:
        if read_type == "dRNA":
            return drna_nd_par(ref_len, basecaller)
        else:
            return cdna_nd_par(ref_len, read_type)
    else:
        return dna_nd_par(ref_len, basecaller)


def seg_par(ref_len, changepoint):
    if ref_len > changepoint:
        xe2 = 0
        xe3 = ref_len - changepoint
    else:
        xe2 = ref_len - changepoint
        xe3 = 0

    xe1 = ref_len - changepoint
    return xe1, xe2, xe3


def dna_nd_par(ref_len, basecaller):
    if basecaller == "albacore":
        at_xe1, at_xe2, at_xe3 = seg_par(ref_len, 12)
        at_mu = 9.32968110 + 0.75215056 * at_xe1 + 0.01234263 * at_xe2 ** 2 - 0.02699184 * at_xe3 ** 2
        at_sigma = 0.2507 * ref_len - 0.1510

        cg_xe1, cg_xe2, cg_xe3 = seg_par(ref_len, 7)
        cg_mu = 4.76783156 + 0.21751512 * cg_xe1 - 0.10414613 * cg_xe2 ** 2 - 0.01647626 * cg_xe3 ** 2
        cg_sigma = 0.1109 * ref_len + 0.8959
    elif basecaller == "guppy":
        at_xe1, at_xe2, at_xe3 = seg_par(ref_len, 16)
        at_mu = 12.13283713 + 0.71843395 * at_xe1 + 0.00127124 * at_xe2 ** 2 - 0.01429113 * at_xe3 ** 2
        at_sigma = 0.2138 * ref_len + 0.4799

        cg_xe1, cg_xe2, cg_xe3 = seg_par(ref_len, 12)
        cg_mu = 6.13526513 + 0.21219760 * cg_xe1 - 0.01249273 * cg_xe2 ** 2 - 0.04870821 * cg_xe3 ** 2
        cg_sigma = 0.3021 * ref_len - 0.06803
    else:  # guppy-flipflop
        at_xe1, at_xe2, at_xe3 = seg_par(ref_len, 12)
        at_mu = 9.65577227 + 0.92524095 * at_xe1 + 0.02814258 * at_xe2 ** 2 - 0.01666699 * at_xe3 ** 2
        at_sigma = 0.2013 * ref_len + 0.2159

        cg_xe1, cg_xe2, cg_xe3 = seg_par(ref_len, 14)
        cg_mu = 5.5402163422 - 0.0163232962 * cg_xe1 - 0.0205230566 * cg_xe2 ** 2 + 0.0009633945 * cg_xe3 ** 2
        cg_sigma = 0.1514 * ref_len + 0.5611

    return at_mu, at_sigma, at_mu, at_sigma, cg_mu, cg_sigma, cg_mu, cg_sigma


def drna_nd_par(ref_len, basecaller):
    if basecaller == "albacore":
        a_xe1, a_xe2, a_xe3 = seg_par(ref_len, 15)
        a_mu = 7.920954923 + 0.060905864 * a_xe1 - 0.031587949 * a_xe2 ** 2 + 0.008346099 * a_xe3 ** 2
        a_sigma = 0.21924 * ref_len + 0.08617

        t_xe1, t_xe2, t_xe3 = seg_par(ref_len, 10)
        t_mu = 5.57952387 + 0.11981070 * t_xe1 - 0.05381813 * t_xe2 ** 2 - 0.01851555 * t_xe3 ** 2
        t_sigma = 0.1438 * ref_len + 0.9004

        c_xe1, c_xe2, c_xe3 = seg_par(ref_len, 6)
        c_mu = 4.34230221 + 0.23685824 * c_xe1 - 0.08072337 * c_xe2 ** 2 - 0.01720295 * c_xe3 ** 2
        c_sigma = 0.06822 * ref_len + 0.87553

        g_xe1, g_xe2, g_xe3 = seg_par(ref_len, 10)
        g_mu = 4.5110033642 + 0.1757467623 * g_xe1 - 0.0009943928 * g_xe2 ** 2 - 0.0915431864 * g_xe3 ** 2
        g_sigma = 0.1066 * ref_len + 0.8223
    else:  # guppy
        a_xe1, a_xe2, a_xe3 = seg_par(ref_len, 16)
        a_mu = 10.619825093 + 0.367474710 * a_xe1 - 0.018351943 * a_xe2 ** 2 - 0.001138652 * a_xe3 ** 2
        a_sigma = 0.3128 * ref_len - 0.6731

        t_xe1, t_xe2, t_xe3 = seg_par(ref_len, 18)
        t_mu = 9.819245514 + 0.212176633 * t_xe1 - 0.016791683 * t_xe2 ** 2 - 0.001310379 * t_xe3 ** 2
        t_sigma = 0.23963 * ref_len + 0.02851

        c_xe1, c_xe2, c_xe3 = seg_par(ref_len, 9)
        c_mu = 4.873485405 + 0.006877753 * c_xe1 - 0.043582044 * c_xe2 ** 2 + 0.018552245 * c_xe3 ** 2
        c_sigma = 0.07976 * ref_len + 0.67479

        g_xe1, g_xe2, g_xe3 = seg_par(ref_len, 7)
        g_mu = 4.431945507 + 0.162662708 * g_xe1 - 0.070326449 * g_xe2 ** 2 + 0.004705819 * g_xe3 ** 2
        g_sigma = 0.08815 * ref_len + 0.81112

    return a_mu, a_sigma, t_mu, t_sigma, c_mu, c_sigma, g_mu, g_sigma


def cdna_nd_par(ref_len, read_type):
    if read_type == "cDNA_1D":
        a_xe1, a_xe2, a_xe3 = seg_par(ref_len, 17)
        a_mu = 12.242619571 + 0.668226076 * a_xe1 + 0.001965846 * a_xe2 ** 2 - 0.026708831 * a_xe3 ** 2
        a_sigma = 0.3703 * ref_len - 0.8779

        t_xe1, t_xe2, t_xe3 = seg_par(ref_len, 14)
        t_mu = 11.979887272 + 1.005918453 * t_xe1 + 0.018442004 * t_xe2 ** 2 - 0.001733806 * t_xe3 ** 2
        t_sigma = 0.2374 * ref_len + 0.0647

        c_xe1, c_xe2, c_xe3 = seg_par(ref_len, 7)
        c_mu = 4.63015250 + 0.02288890 * c_xe1 - 0.13258183 * c_xe2 ** 2 + 0.01859354 * c_xe3 ** 2
        c_sigma =  0.1805 * ref_len + 0.3770

        g_xe1, g_xe2, g_xe3 = seg_par(ref_len, 6)
        g_mu = 4.26732085 + 0.18475982 * g_xe1 - 0.14151564 * g_xe2 ** 2 - 0.03719026 * g_xe3 ** 2
        g_sigma = 0.1065 * ref_len + 0.8318
    else:  # cDNA 1D2
        a_xe1, a_xe2, a_xe3 = seg_par(ref_len, 16)
        a_mu = 14.807216375 + 0.957883839 * a_xe1 + 0.003869761 * a_xe2 ** 2 - 0.027664685 * a_xe3 ** 2
        a_sigma = 0.1842 * ref_len + 0.5642

        t_xe1, t_xe2, at_xe3 = seg_par(ref_len, 12)
        t_mu = 11.281353708 + 1.039742281 * t_xe1 + 0.015074972 * t_xe2 ** 2 - 0.004161978 * at_xe3 ** 2
        t_sigma = 0.1842 * ref_len + 0.5642

        c_xe1, c_xe2, c_xe3 = seg_par(ref_len, 11)
        c_mu = 6.925880017 + 0.383403249 * c_xe1 - 0.003611025 * c_xe2 ** 2 - 0.038495472 * c_xe3 ** 2
        c_sigma = 0.32113 * ref_len - 0.04017

        g_xe1, g_xe2, g_xe3 = seg_par(ref_len, 14)
        g_mu = 9.10846385 + 0.68454921 * g_xe1 + 0.01769293 * g_xe2 ** 2 - 0.07252329 * g_xe3 ** 2
        g_sigma = 0.32113 * ref_len - 0.04017

    return a_mu, a_sigma, t_mu, t_sigma, c_mu, c_sigma, g_mu, g_sigma


def get_hpmis_rate(read_type, basecaller):
    if read_type:
        if read_type == "dRNA":
            if basecaller == "albacore":
                return 0.041483
            elif basecaller == "guppy":
                return 0.027234
        elif read_type == "cDNA_1D":
            return 0.036122
        else:  # cDNA 1D2
            return 0.040993
    else:  #DNA
        if basecaller == "albacore":
            return 0.02204
        elif basecaller == "guppy":
            return 0.02166
        elif basecaller == "guppy-flipflop":
            return 0.02215

