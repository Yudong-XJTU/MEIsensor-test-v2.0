#!/usr/bin/env python3
#
# sv
from dataclasses import dataclass
import statistics
from typing import Optional

import utils

TYPES = ["INS", "DEL"]



@dataclass
class SVCallPostprocessingInfo:
    cluster: list


@dataclass
class SVCall:
    contig: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: int
    filter: str
    info: dict

    svtype: str
    svlen: int
    end: int
    genotypes: dict[int, tuple]

    precise: bool
    support: int
    rnames: list[str]

    qc: bool
    nm: float
    postprocess: Optional[SVCallPostprocessingInfo]

    svlens: list[int] = None

    fwd: int = None
    rev: int = None

    coverage_upstream: int = 0
    coverage_downstream: int = 0
    coverage_start: int = 0
    coverage_center: int = 0
    coverage_end: int = 0

    sample_internal_id: int = None

    raw_vcf_line: Optional[str] = None
    raw_vcf_line_index: Optional[int] = None

    _counter = 0

    def set_info(self, k, v):
        self.info[k] = v

    def get_info(self, k):
        return self.info[k] if k in self.info else None

    def has_info(self, k):
        return k in self.info

    def finalize(self):
        self.postprocess = None


@dataclass
class Cluster:
    id: str
    svtype: str
    contig: str
    start: int
    end: int
    seed: int
    signals: list
    signals_long: Optional[list]

    @property
    def span(self) -> Optional[int]:
        if self.end is None or self.start is None:
            return None

        return self.end - self.start

    def compute_metrics(self, max_n=100):
        n = min(len(self.signals), max_n)
        if n == 0:
            self.mean_svlen = 0
            self.stdev_start = 0
            return

        step = int(len(self.signals) / n)
        if n > 1:
            self.mean_svlen = sum(self.signals[i].svlen for i in range(0, len(self.signals), step)) / float(n)
            self.stdev_start = statistics.stdev(self.signals[i].ref_start for i in range(0, len(self.signals), step))
        else:
            self.mean_svlen = self.signals[0].svlen
            self.stdev_start = 0


def merge_inner(cluster, threshold):
    read_seq = {}
    for ld in cluster.signals:
        if ld.read_qname not in read_seq:
            read_seq[ld.read_qname] = []
        read_seq[ld.read_qname].append(ld)

    cluster.signals = []
    for qname in read_seq:
        read_seq[qname].sort(key=lambda k: k.ref_start)
        to_merge = read_seq[qname][0]

        curr_signal = to_merge

        last_ref_end = to_merge.ref_end
        last_qry_end = to_merge.qry_end
        last_ref_start = to_merge.ref_start
        last_qry_start = to_merge.qry_start

        for to_merge in read_seq[qname][1:]:
            merge = (threshold == -1) or ((abs(to_merge.ref_start - last_ref_end) < threshold or abs(to_merge.ref_start - last_ref_start) < threshold) and (
                        abs(to_merge.qry_start - last_qry_end) < threshold or abs(to_merge.qry_start - last_qry_start) < threshold))
            if merge:
                curr_signal.svlen += to_merge.svlen
                if to_merge.seq is None or curr_signal.seq is None:
                    curr_signal.seq = None
                else:
                    curr_signal.seq += to_merge.seq
            else:
                cluster.signals.append(curr_signal)
                curr_signal = to_merge
            last_ref_end = to_merge.ref_end
            last_qry_end = to_merge.qry_end
            last_ref_start = to_merge.ref_start
            last_qry_start = to_merge.qry_start

        cluster.signals.append(curr_signal)
    return cluster


def resplit(cluster, prop, binsize, merge_threshold_min, merge_threshold_frac):
    bins_signals = {}
    for signal in cluster.signals:
        bin = int(abs(prop(signal)) / binsize) * binsize
        if not bin in bins_signals:
            bins_signals[bin] = [signal]
        else:
            bins_signals[bin].append(signal)

    new_clusters = list(sorted(bins_signals.keys()))
    i = 1
    while len(new_clusters) > 1 and i < len(new_clusters):
        last_cluster = new_clusters[i - 1]
        curr_cluster = new_clusters[i]
        merge_threshold = max(merge_threshold_min, min(curr_cluster, last_cluster) * merge_threshold_frac)
        merge = abs(curr_cluster - last_cluster) <= merge_threshold
        if merge:
            bins_signals[new_clusters[i]].extend(bins_signals[new_clusters[i - 1]])
            new_clusters.pop(i - 1)
            i = max(0, i - 2)
        else:
            i += 1

    for cluster_index in new_clusters:
        new_cluster = Cluster(id=cluster.id + f".{cluster_index}",
                              svtype=cluster.svtype,
                              contig=cluster.contig,
                              start=cluster.start,
                              end=cluster.end,
                              seed=cluster.seed,
                              signals=bins_signals[cluster_index],
                              signals_long=cluster.signals_long)
        yield new_cluster



def resolve(svtype, signal_provider, config):
    signaltab = signal_provider.signaltab[svtype]
    seeds = sorted(signal_provider.signaltab[svtype])

    if len(seeds) == 0:
        return []

    clusters = []
    for seed_index, seed in enumerate(seeds):
        if config.dev_call_region != None:
            if seed < config.dev_call_region["start"] or seed > config.dev_call_region["end"]:
                continue
        if svtype == "INS":
            signals = [signal for signal in signaltab[seed] if signal.svlen != None]
            signals_long = [signal for signal in signaltab[seed] if signal.svlen == None]

        else:
            signals = signaltab[seed]
            signals_long = None

        cluster = Cluster(id=f"CL.{svtype}.{signal_provider.contig}.{signal_provider.start}.{seed_index}",
                          svtype=svtype,
                          contig=signal_provider.contig,
                          start=seed,
                          end=seed + config.cluster_binsize,
                          seed=seed,
                          signals=signals,
                          signals_long=signals_long)

        cluster.compute_metrics()
        clusters.append(cluster)

    # Merge clusters
    cluster_count_initial = len(clusters)
    i = 0
    while i < len(clusters) - 1:
        curr_cluster = clusters[i]
        next_cluster = clusters[i + 1]

        inner_dist = (next_cluster.start - curr_cluster.end)
        outer_dist = (next_cluster.end - curr_cluster.start)
        merge = inner_dist <= min(curr_cluster.stdev_start, next_cluster.stdev_start) * config.cluster_r
        merge = merge

        if merge:
            clusters.pop(i + 1)
            curr_cluster.signals += next_cluster.signals
            if svtype == "INS":
                curr_cluster.signals_long += next_cluster.signals_long
            curr_cluster.end = next_cluster.end
            curr_cluster.compute_metrics()
            i = max(0, i - 2)
        i += 1

    if config.dev_trace_read:
        for c in clusters:
            for ld in c.signals:
                if ld.read_qname == config.dev_trace_read:
                    print(f"[DEV_TRACE_READ [2/4] [cluster.resolve] Read signal {ld} is in cluster {c.id}, containing a total of {len(c.signals)} signals")

    if config.dev_dump_clusters:
        filename = f"{config.input}.clusters.{svtype}.{signal_provider.contig}.{signal_provider.start}.{signal_provider.end}.bed"
        print(f"Dumping clusters to {filename}")
        with open(filename, "w") as h:
            for c in clusters:
                info = f"ID={c.id}, #SIGNALS={len(c.signals)}; "
                for ld in c.signals:
                    info += f"(ref_start={ld.ref_start},svlen={ld.svlen},source={ld.source}); "
                h.write(f"{c.contig}\t{c.start}\t{c.end}\t\"{info}\"\n")

    for cluster in clusters:
        if len(cluster.signals) == 0:
            continue
        if svtype == "INS":
            merge_inner_threshold = config.cluster_merge_pos
            merge_inner(cluster, merge_inner_threshold)

        if not config.dev_no_resplit:
            for new_cluster in resplit(cluster,
                                       prop=lambda signal: signal.svlen,
                                       binsize=config.cluster_resplit_binsize,
                                       merge_threshold_min=config.minsvlen,
                                       merge_threshold_frac=config.cluster_merge_len):
                yield new_cluster
        else:
            yield cluster


def calculate_bounds(svtype, ref_start_mode, svlen_mode):
    if svtype == "INS":
        svstart = ref_start_mode
        svend = ref_start_mode
    elif svtype == "DEL":
        svstart = ref_start_mode + svlen_mode
        svend = ref_start_mode
    else:
        svstart = ref_start_mode
        svend = svstart + abs(svlen_mode)
    return svstart, svend


def call_from(cluster, config, keep_qc_fails, task):
    signals = cluster.signals

    svtype = cluster.svtype
    stdev_pos, stdev_len = None, None
    qc = True

    svlen = utils.center(v.svlen for v in signals)

    if config.dev_emit_sv_lengths:
        svlens = [v.svlen for v in signals]
    else:
        svlens = None

    if abs(svlen) < config.minsvlen_screen:
        return

    # Count inline events only once per read, but split events as individual alignments, as in coverage calculation
    # inline_qnames=set(k.read_qname for k in signals if k.source=="INLINE")
    # support=len(inline_qnames)+sum(1 for k in signals if k.source!="INLINE")
    support_set = set(k.read_qname for k in signals)
    if svtype == "INS" and svlen >= config.long_ins_length:
        support_long_set = set(signal.read_qname for signal in cluster.signals_long)
        support_long = len(support_long_set)
        support_set |= support_long_set
        support = len(support_set)
    else:
        support = len(support_set)
        support_long = 0
    ref_start = utils.center(v.ref_start for v in signals)
    stdev_pos = utils.stdev(utils.trim((v.ref_start for v in signals)))
    precise = (stdev_pos < config.precise)
    svstart, svend = calculate_bounds(svtype, ref_start, svlen)
    qual = int(utils.mean(v.mapq for v in signals))

    support_fwd = sum(signal.strand == "+" for signal in signals)
    support_rev = len(signals) - support_fwd

    filter_sv = "PASS"

    if config.qc_nm_measure:
        nm_mean = utils.mean(v.nm for v in signals)
    else:
        nm_mean = -1

    if not keep_qc_fails and not qc:
        return

    svpi = SVCallPostprocessingInfo(cluster=cluster)

    if config.output_rnames:
        rnames = list(support_set)
    else:
        rnames = None

    svcall = SVCall(contig=cluster.contig,
                    pos=svstart,
                    id=f"{svtype}.{task.sv_id:X}S{task.id:X}",
                    ref="N",
                    alt=f"<{svtype}>",
                    qual=qual,
                    filter=filter_sv,
                    info=dict(),
                    svtype=svtype,
                    svlen=svlen,
                    svlens=svlens,
                    end=svend,
                    genotypes=dict(),
                    precise=precise,
                    support=support,
                    rnames=rnames,
                    postprocess=svpi,
                    qc=qc,
                    nm=nm_mean,
                    fwd=support_fwd,
                    rev=support_rev)

    if svtype == "INS":
        svcall.set_info("SUPPORT_LONG", support_long)

    if stdev_pos is not None:
        svcall.set_info("STDEV_POS", stdev_pos)
    if stdev_len is not None:
        svcall.set_info("STDEV_LEN", stdev_len)

    # svcall.set_info("CLUSTER_ID",cluster.id)

    task.sv_id += 1
    yield svcall


def classify_splits(read, signals, config, main_contig):
    """
    Determines the SV type of a split read (read with supplementary alignments)
    """
    minsvlen_screen = config.minsvlen_screen
    maxsvlen_other = minsvlen_screen * config.dev_split_max_query_distance_mult

    signals.sort(key=lambda ld: ld.qry_start)
    last = signals[0]
    last.svtypes_starts_lens = []

    if last.qry_start >= config.long_ins_length * 0.5:
        last.svtypes_starts_lens.append(("INS", last.ref_start, None))

    for i in range(1, len(signals)):
        curr = signals[i]
        curr.svtypes_starts_lens = []
        qry_dist_abs = abs(curr.qry_start - last.qry_end)

        if curr.contig == last.contig:
            rev = (curr.strand == "-")
            fwd = not rev
            if curr.strand == last.strand:
                #
                # INS, DEL
                #
                if (fwd and (curr.qry_start - last.qry_end) >= minsvlen_screen
                        and (curr.ref_start - last.ref_end) < maxsvlen_other
                        and (curr.qry_start - last.qry_end) - (curr.ref_start - last.ref_end) >= minsvlen_screen):
                    # INS, FWD
                    svstart = curr.ref_start
                    svlen = (curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq = read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq = None
                    curr.svtypes_starts_lens.append(("INS", svstart, svlen))

                elif (rev and (curr.qry_start - last.qry_end) >= minsvlen_screen
                      and (last.ref_start - curr.ref_end) < maxsvlen_other
                      and (curr.qry_start - last.qry_end) - (last.ref_start - curr.ref_end) >= minsvlen_screen):
                    # INS, REV
                    svstart = last.ref_start
                    svlen = (curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq = read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq = None
                    curr.svtypes_starts_lens.append(("INS", svstart, svlen))

                elif (fwd and (curr.ref_start - last.ref_end) >= minsvlen_screen
                        and qry_dist_abs < maxsvlen_other
                        and (curr.ref_start - last.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                    # DEL, FWD
                    svstart = curr.ref_start
                    svlen = (curr.ref_start - last.ref_end)
                    curr.svtypes_starts_lens.append(("DEL", svstart, -svlen))

                elif (rev and (last.ref_start - curr.ref_end) >= minsvlen_screen
                      and qry_dist_abs < maxsvlen_other
                      and (last.ref_start - curr.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                    # DEL, REV
                    svstart = last.ref_start
                    svlen = (last.ref_start - curr.ref_end)
                    curr.svtypes_starts_lens.append(("DEL", svstart, -svlen))
        last = curr
