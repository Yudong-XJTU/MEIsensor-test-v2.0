#!/usr/bin/env python3
#
# postprocessing

import math

import utils
from config import Config
from SV import SVCall


def annotate_sv(svcall: SVCall, config):
    genotype_sv(svcall, config)

    if svcall.svtype == "INS" and not config.symbolic:
        merged_signals = [signal for signal in svcall.postprocess.cluster.signals if signal.seq is not None]

        if len(merged_signals):
            best_signal = merged_signals[0]
            best_index = 0
            best_diff = abs(len(best_signal.seq) - svcall.svlen) + abs(best_signal.ref_start - svcall.pos) * 1.5
            for i, ld in enumerate(merged_signals):
                if i == 0:
                    continue
                curr_diff = abs(len(ld.seq) - svcall.svlen) + abs(ld.ref_start - svcall.pos) * 1.5
                if curr_diff < best_diff:
                    best_signal = ld
                    best_index = i
                    best_diff = curr_diff

            merged_signals.pop(best_index)

            if len(merged_signals) >= config.consensus_min_reads and not config.no_consensus:
                kmer_len = config.consensus_kmer_len
                skip = config.consensus_kmer_skip_base + int(len(best_signal.seq)*config.consensus_kmer_skip_seqlen_mult)
                skip_repetitive = skip

                svcall.alt = novel_from_reads(best_signal, merged_signals, klen=kmer_len, skip=skip,
                                              skip_repetitive=skip_repetitive)
            else:
                svcall.alt = best_signal.seq


def add_request(svcall, field, pos, requests_for_coverage, config):
    bin_index = int(pos / config.coverage_binsize) * config.coverage_binsize
    if bin_index not in requests_for_coverage:
        requests_for_coverage[bin_index] = []
    requests_for_coverage[bin_index].append((svcall, field))


def coverage(calls, signal_provider, config):
    requests_for_coverage = coverage_build_requests(calls, config)
    return coverage_fulfill(requests_for_coverage, signal_provider, config)


def coverage_build_requests(calls, config: Config):
    """
    Updated: changed coverage calculation to be inside the SV for INV/DEL/DUP, outside (like it is now) for INS/BND?
    BIN: ---
    INS/BND    |-------------V---------------|
              U---       S--- E---            D---
                        C---
    DUP/DEL    |------------===============------------|
    INV       U---          S---       E---            D---
                                  C---
    """
    requests_for_coverage = {}
    for svcall in calls:
        start = svcall.pos
        if svcall.svtype == "INS":
            end = start + 1
        else:
            end = svcall.pos + abs(svcall.svlen)
        if svcall.svtype in ["INS", "BND"]:
            add_request(svcall, "coverage_start", start - config.coverage_binsize, requests_for_coverage, config)
            add_request(svcall, "coverage_center", int((start + end - config.coverage_binsize) / 2),
                        requests_for_coverage, config)
            add_request(svcall, "coverage_end", end + config.coverage_binsize, requests_for_coverage, config)
            add_request(svcall, "coverage_upstream", start - config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage, config)
            add_request(svcall, "coverage_downstream", end + config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage, config)
        else:
            add_request(svcall, "coverage_start", start, requests_for_coverage, config)
            add_request(svcall, "coverage_center", int((start + end - config.coverage_binsize) / 2),
                        requests_for_coverage, config)
            add_request(svcall, "coverage_end", end - config.coverage_binsize, requests_for_coverage, config)
            add_request(svcall, "coverage_upstream", start - config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage, config)
            add_request(svcall, "coverage_downstream", end + config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage, config)

    return requests_for_coverage


def coverage_fulfill(requests_for_coverage, signal_provider, config: Config):
    if len(requests_for_coverage) == 0:
        return -1, -1

    start_bin = signal_provider.coverage_bins_min_bin
    end_bin = int(signal_provider.end / config.coverage_binsize) * config.coverage_binsize
    coverage_fwd = 0
    coverage_rev = 0
    coverage_fwd_total = 0
    coverage_rev_total = 0
    n = 0

    for bin_index in range(start_bin, end_bin + config.coverage_binsize, config.coverage_binsize):
        n += 1

        if bin_index in signal_provider.coverage_bins_fwd:
            coverage_fwd += signal_provider.coverage_bins_fwd[bin_index]

        if bin_index in signal_provider.coverage_bins_rev:
            coverage_rev += signal_provider.coverage_bins_rev[bin_index]

        if bin_index in requests_for_coverage:
            coverage_total_curr = coverage_fwd + coverage_rev
            for svcall, field in requests_for_coverage[bin_index]:
                setattr(svcall, field, coverage_total_curr)

        coverage_fwd_total += coverage_fwd
        coverage_rev_total += coverage_rev

    average_coverage_fwd = coverage_fwd_total / float(n) if n > 0 else 0
    average_coverage_rev = coverage_rev_total / float(n) if n > 0 else 0
    return average_coverage_fwd, average_coverage_rev


def qc_sv_support(svcall, coverage_global, config):
    if config.mosaic:
        return True
    if config.minsupport == "auto":
        if not qc_support_auto(svcall, coverage_global, config):
            svcall.filter = "SUPPORT_MIN"
            return False
    else:
        if not qc_support_const(svcall, config):
            svcall.filter = "SUPPORT_MIN"
            return False
    return True


def rescale_support(svcall, config: Config):
    if svcall.svtype != "INS" or svcall.svlen < config.long_ins_length:
        return svcall.support
    else:
        base = svcall.support
        scale_factor = config.long_ins_rescale_mult * (float(svcall.svlen) / config.long_ins_length)
        return round(base * (config.long_ins_rescale_base + scale_factor))


def qc_support_auto(svcall, coverage_global, config):
    support = rescale_support(svcall, config)
    # if svcall.svtype=="INS":
    #    coverage_list=[svcall.coverage_center]
    # else:
    coverage_list = [each_coverage for each_coverage in [svcall.coverage_upstream, svcall.coverage_downstream] if
                     each_coverage != 0]
    if len(coverage_list) == 0:
        coverage_list = [each_coverage for each_coverage in [svcall.coverage_start, svcall.coverage_center,
                                                             svcall.coverage_end] if each_coverage != 0]
    if len(coverage_list) == 0:
        coverage_regional = coverage_global
    else:
        coverage_regional = round(sum(coverage_list) / len(coverage_list))
        if coverage_regional == 0:
            coverage_regional = coverage_global
    coverage_global_weight = (1.0 - config.minsupport_auto_regional_coverage_weight)
    coverage = (coverage_regional * config.minsupport_auto_regional_coverage_weight +
                coverage_global * coverage_global_weight)
    min_support = round(config.minsupport_auto_base + config.minsupport_auto_mult * coverage)
    return support >= min_support
    # return True


def qc_support_const(svcall, config):
    # svcall.set_info("MINSUPPORT",config.minsupport)
    return svcall.support >= config.minsupport


def qc_sv(svcall: SVCall, config: Config):
    af = svcall.get_info("VAF")
    af = af if af is not None else 0
    sv_is_mosaic = af <= config.mosaic_af_max

    if config.qc_stdev:
        stdev_pos = svcall.get_info("STDEV_POS")
        if stdev_pos > config.qc_stdev_abs_max:
            svcall.filter = "STDEV_POS"
            return False
        if svcall.svtype != "BND" and stdev_pos / abs(svcall.svlen) > 2.0:
            svcall.filter = "STDEV_POS"
            return False

        stdev_len = svcall.get_info("STDEV_LEN")
        if stdev_len is not None:
            if svcall.svtype != "BND" and stdev_len / abs(svcall.svlen) > 1.0:
                svcall.filter = "STDEV_LEN"
                return False
            if stdev_len > config.qc_stdev_abs_max:
                svcall.filter = "STDEV_LEN"
                return False

    if config.mosaic and sv_is_mosaic:
        if svcall.support < config.mosaic_min_reads:
            svcall.filter = "SUPPORT_MIN"
            return False

    if abs(svcall.svlen) < config.minsvlen:
        svcall.filter = "SVLEN_MIN"
        return False

    if not (config.mosaic and sv_is_mosaic) and config.qc_strand:
        is_long_ins = (svcall.svtype == "INS" and svcall.svlen >= config.long_ins_length)
        if not is_long_ins and len(set(signals.strand for signals in svcall.postprocess.cluster.signals)) < 2:
            svcall.filter = "STRAND"
            return False
    elif (config.mosaic and sv_is_mosaic) and config.mosaic_qc_strand:
        is_long_ins = (svcall.svtype == "INS" and svcall.svlen >= config.long_ins_length)
        if (not is_long_ins and len(set(signals.strand for signals in svcall.postprocess.cluster.signals)) < 2
                and svcall.support >= config.mosaic_use_strand_thresholds):
            svcall.filter = "STRAND_MOSAIC"
            return False

    upstream_downstream_max_coverage_diff = 0.7  # 70%
    upstream_downstream_diff = 0.5  # add to config
    if svcall.svtype == "INS" and (svcall.coverage_upstream < config.qc_coverage or
                                     svcall.coverage_downstream < config.qc_coverage):
        svcall.filter = "COV_CHANGE_INS"
        return False

    qc_coverage_max_change_frac = config.qc_coverage_max_change_frac
    if config.mosaic and sv_is_mosaic:
        qc_coverage_max_change_frac = config.mosaic_qc_coverage_max_change_frac
    if qc_coverage_max_change_frac != -1.0:
        if svcall.coverage_upstream != 0:
            u = float(svcall.coverage_upstream)
        else:
            u = 1.0

        if svcall.coverage_start != 0:
            s = float(svcall.coverage_start)
        else:
            s = 1.0

        if svcall.coverage_center != 0:
            c = float(svcall.coverage_center)
        else:
            c = 1.0

        if svcall.coverage_end != 0:
            e = float(svcall.coverage_end)
        else:
            e = 1.0

        if svcall.coverage_downstream != 0:
            d = float(svcall.coverage_downstream)
        else:
            d = 1.0

        if abs(u - s) / max(u, s) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_US"
            return False

        if abs(s - c) / max(s, c) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_SC"
            return False

        if abs(c - e) / max(c, e) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_CE"
            return False

        if abs(e - d) / max(e, d) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_ED"
            return False

    return True


def qc_sv_post_annotate(svcall: SVCall, config: Config):
    af = svcall.get_info("VAF")
    af = af if af is not None else 0
    sv_is_mosaic = af <= config.mosaic_af_max

    if (svcall.coverage_center < config.qc_coverage and
            (len(svcall.genotypes) == 0 or (svcall.genotypes[0][0] != "." and
                                            svcall.genotypes[0][0] + svcall.genotypes[0][1] < 2))):
        svcall.filter = "COV_MIN_GT"
        return False

    qc_nm = config.qc_nm
    qc_nm_threshold = config.qc_nm_threshold * config.qc_nm_mult
    if config.mosaic and sv_is_mosaic:
        qc_nm = config.mosaic_qc_nm
        qc_nm_threshold = config.qc_nm_threshold * config.qc_nm_mult
    if qc_nm and svcall.nm > qc_nm_threshold and (len(svcall.genotypes) == 0 or svcall.genotypes[0][1] == 0):
        svcall.filter = "ALN_NM"
        return False

    if config.mosaic:
        if sv_is_mosaic and (af < config.mosaic_af_min or af > config.mosaic_af_max):
            svcall.filter = "MOSAIC_VAF"
            return False
        elif not sv_is_mosaic and not config.mosaic_include_germline:
            svcall.filter = "NOT_MOSAIC_VAF"
            return False

    return True


def binomial_coef(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


#
# Consensus
#

def iter_kmers(seq, klen, skip):
    for i in range(0, len(seq) - klen, skip):
        yield (i, seq[i:i + klen])


def novel_from_reads(best_signal, other_signals, klen, skip, skip_repetitive, debug=False):
    consensus_min = 2
    maxshift = klen
    minspan = 0.2
    minalns = 0.25
    minident = 0.5
    minident_abs = 5
    minbestdiff = 3

    alignments = []
    anchors = {}
    taboo = set()
    for i, kmer in iter_kmers(best_signal.seq, klen=klen, skip=skip_repetitive):
        if kmer in taboo:
            continue
        if kmer in anchors:
            del anchors[kmer]
            taboo.add(kmer)
            continue
        anchors[kmer] = i

    for signal_i, signal in enumerate(other_signals):
        last_i = None
        last_j = None
        conseq = ""
        span = 0
        for j, kmer in iter_kmers(signal.seq, klen=klen, skip=skip):
            if kmer not in anchors:
                continue
            i = anchors[kmer]
            if abs(i - j) > maxshift:
                continue
            if last_i is not None and i <= last_i:
                continue

            if last_i is None:
                if j > 0:
                    conseq = "-" * i
            else:
                fwd_i = i - last_i
                fwd_j = j - last_j
                if len(conseq) + fwd_j > len(best_signal.seq):
                    fwd_j = len(best_signal.seq) - len(conseq)

                if fwd_i == fwd_j and fwd_j > 0:
                    span += (j - last_j)
                    m = 0
                    for l in range(1, (j - last_j) + 1):
                        if signal.seq[last_j + l] == best_signal.seq[last_i + l]:
                            m += 1
                    ident = m / float((j - last_j))
                    if ident >= minident:
                        conseq += signal.seq[last_j:j][:fwd_j]
                    else:
                        conseq += "-" * (fwd_j)
                else:
                    conseq += "-" * (fwd_j)
            last_i = i
            last_j = j

        if len(conseq) < len(best_signal.seq):
            conseq += "-" * (len(best_signal.seq) - len(conseq))

        conseq_new = []
        h = 0
        while h < len(best_signal.seq):
            if conseq[h] == "-":
                conseq_new.append("-")
                h += 1
            else:
                buffer = []
                ident = 0
                while h < len(best_signal.seq) and conseq[h] != "-":
                    ident += (best_signal.seq[h] == conseq[h])
                    buffer.append(conseq[h])
                    h += 1
                if ident / float(len(buffer)) > minident and ident > minident_abs:
                    conseq_new.append("".join(buffer))
                else:
                    conseq_new.append("-" * len(buffer))
        conseq = "".join(conseq_new)

        if span / float(len(best_signal.seq)) > minspan:
            alignments.append(conseq)

    maxal = 1
    for i in range(len(best_signal.seq)):
        maxal = max(maxal, len([best_signal.seq[i]] + [a[i] for a in alignments if not a[i] in "^_"]))
    maxal = float(maxal)

    flattened = ""
    for i in range(len(best_signal.seq)):
        al = [a[i] for a in alignments if not a[i] == "-"]
        if len(al) < consensus_min or len(al) / maxal < minalns:
            flattened += best_signal.seq[i]
        else:
            top = utils.most_common([best_signal.seq[i]] + al)
            if len(top) > 1 and top[0][0] - top[1][0] >= minbestdiff:
                flattened += top[0][1]
            else:
                flattened += best_signal.seq[i]

            # if top[0][0]/float(len(al)+1) < 0.75:
            #    flattened+=best_signal.seq[i]
            # else:
            #    flattened+=top[0][1]

    # print("FLT",flattened)
    # if debug:
    #    print("B",best_signal.seq)
    #    for o in other_signals:
    #        print("O",o.seq)
    #    print("F",flattened)
    #    print("=====")
    return flattened


#
# Genotyping
#

class UnknownGenotype(Exception):
    """
    Unable to determine genotype
    """


def binomial_probability(k, n, p):
    try:
        # Binomial coef cancels out for likelihood ratios
        return (p ** k) * ((1.0 - p) ** (n - k))
    except OverflowError:
        return 1.0


def likelihood_ratio(q1, q2):
    if q1 / q2 > 0:
        try:
            return math.log(q1 / q2, 10)
        except ValueError:
            return 0
    else:
        return 0


class Genotyper:
    """
    Generic genotyping
    """
    _support: int
    _coverage: float

    def __init__(self, svcall: SVCall, config):
        self.svcall = svcall
        self.config = config



    def _calculate_support(self) -> int:
        return self.svcall.support

    def _calculate_coverage(self, support: int) -> int:
        return self._get_coverage_from_list()

    def _calculate_af(self, support: int, coverage: int) -> float:
        """
        Calculate allele fraction
        """
        return support / float(coverage)

    def _get_coverage_from_list(self, coverage_list: list = None) -> int:
        """
        Coverage here is NOT the same as in IGV, but the number of reads spanning the SV
        """
        svcall = self.svcall
        if coverage_list is None:
            coverage_list = [svcall.coverage_start, svcall.coverage_center, svcall.coverage_end]

        coverage_list = [each_coverage for each_coverage in coverage_list if each_coverage != 0]

        if len(coverage_list) > 0:
            return round(sum(coverage_list) / len(coverage_list))
        else:
            raise UnknownGenotype()

    def _filter_by_z_score(self, z_score: float) -> bool:
        """
        Should this call be filtered due to z score?
        """
        return z_score < self.config.genotype_min_z_score and not self.config.mosaic

    def calculate(self):
        config = self.config
        normalization_target = 250
        hom_ref_p = config.genotype_error
        het_p = (1.0 / config.genotype_ploidy)  # - config.genotype_error
        hom_var_p = 1.0 - config.genotype_error
        svcall = self.svcall

        support = self._calculate_support()
        try:
            coverage = self._calculate_coverage(support)
        except UnknownGenotype:
            return

        if support > coverage:
            coverage = support

        af = self._calculate_af(support, coverage)

        genotype_p = [((0, 0), hom_ref_p),
                      ((0, 1), het_p),
                      ((1, 1), hom_var_p)]

        max_signal = max(support, coverage)
        if max_signal > normalization_target:
            norm = normalization_target / float(max_signal)
            normalized_support = round(support * norm)
            normalized_coverage = round(coverage * norm)
        else:
            normalized_support = support
            normalized_coverage = coverage

        genotype_likelihoods = []
        for gt, p in genotype_p:
            q = binomial_probability(normalized_support, normalized_coverage, p)
            genotype_likelihoods.append((gt, q))
        genotype_likelihoods.sort(key=lambda k: k[1], reverse=True)

        sum_likelihoods = sum(q for gt, q in genotype_likelihoods)
        normalized_likelihoods = [(gt, (q / sum_likelihoods)) for gt, q in genotype_likelihoods]

        gt1, q1 = normalized_likelihoods[0]
        gt2, q2 = normalized_likelihoods[1]
        qz = [q for gt, q in normalized_likelihoods if gt == (0, 0)][0]
        genotype_z_score = min(60, int((-10) * likelihood_ratio(qz, q1)))
        genotype_quality = min(60, int((-10) * likelihood_ratio(q2, q1)))

        if svcall.filter == "PASS" and self._filter_by_z_score(genotype_z_score):
            svcall.filter = "GT"

        a, b = gt1
        svcall.genotypes[0] = (a, b, genotype_quality, coverage - support, support)
        svcall.set_info("VAF", af)


class InsertionGenotyper(Genotyper):
    def _calculate_support(self):
        """
        Rescale support for long insertions. This skews our support way higher than the number of reads we have!
        """
        return rescale_support(self.svcall, self.config)

    def _calculate_coverage(self, coverage_list: list = None) -> float:
        return self._get_coverage_from_list([self.svcall.coverage_center])

    def _filter_by_z_score(self, z_score: float) -> bool:
        """
        When detecting large insertions, we allow for below threshold z scores
        """
        flt = super()._filter_by_z_score(z_score)
        if flt and self.svcall.svlen >= self.config.long_ins_length and self.config.detect_large_ins:
            return False
        return flt


class DeletionGenotyper(Genotyper):
    ...


GENOTYPER_BY_TYPE = {
    'INS': InsertionGenotyper,
    'DEL': DeletionGenotyper,
}


def genotype_sv(svcall: SVCall, config):
    GENOTYPER_BY_TYPE.get(svcall.svtype, Genotyper)(svcall, config).calculate()
