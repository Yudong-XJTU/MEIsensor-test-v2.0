#!/usr/bin/env python3
#
# preprocessing
from dataclasses import dataclass
from typing import Optional

import pysam

import utils
import SV
from utils import Region


@dataclass
class Signal:
    read_id: int = None
    read_qname: str = None
    contig: str = None
    ref_start: int = None
    ref_end: int = None
    qry_start: int = None
    qry_end: int = None
    strand: str = None
    mapq: int = None
    nm: float = None
    source: str = None
    svtype: str = None
    svlen: Optional[int] = None
    seq: Optional[str] = None
    svtypes_starts_lens: list = None


def CIGAR_analyze(cigar):
    buf = ""
    readspan = 0
    refspan = 0
    clip_start = None
    clip = 0
    for c in cigar:
        if c.isnumeric():
            buf += c
        else:
            oplen = int(buf)
            h = False
            if c in "MIX=":
                readspan += oplen
                h = True
            if c in "MDX=N":
                refspan += oplen
                h = True
            if not h:
                if c in "SH":
                    if clip_start == None and readspan + refspan > 0:
                        clip_start = clip
                    clip += oplen
                else:
                    raise f"Unknown CIGAR operation: '{c}'"
            buf = ""
    if clip_start == None:
        clip_start = clip
    return clip_start, clip - clip_start, refspan, readspan


OPTAB = {pysam.CMATCH: (1, 1, 0),
         pysam.CEQUAL: (1, 1, 0),
         pysam.CDIFF: (1, 1, 0),
         pysam.CINS: (1, 0, 1),
         pysam.CDEL: (0, 1, 1),
         pysam.CREF_SKIP: (0, 1, 0),
         pysam.CSOFT_CLIP: (1, 0, 1),
         pysam.CHARD_CLIP: (0, 0, 0),
         pysam.CPAD: (0, 0, 0)}

OPLIST = [(0, 0, 0) for i in range(max(int(k) for k in OPTAB.keys()) + 1)]
for k, v in OPTAB.items():
    OPLIST[int(k)] = v


def read_iterindels(read_id, read, contig, config, use_clips, read_nm):
    minsvlen = config.minsvlen_screen
    longinslen = config.long_ins_length / 2.0
    seq_cache_maxlen = config.dev_seq_cache_maxlen
    qname = read.query_name
    mapq = read.mapping_quality
    strand = "-" if read.is_reverse else "+"
    CINS = pysam.CINS
    # CDEL = pysam.CDEL
    CSOFT_CLIP = pysam.CSOFT_CLIP

    pos_read = 0
    pos_ref = read.reference_start

    for op, oplength in read.cigartuples:
        add_read, add_ref, event = OPLIST[op]
        if event and oplength >= minsvlen:
            if op == CINS:
                yield Signal(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read + oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           oplength,
                           seq=read.query_sequence[pos_read:pos_read + oplength] if oplength <= seq_cache_maxlen else None)
            elif use_clips and op == CSOFT_CLIP and oplength >= longinslen:
                yield Signal(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read + oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           None,
                           seq=None)
        pos_read += add_read * oplength
        pos_ref += add_ref * oplength


def get_cigar_indels(read_id, read, contig, config, use_clips, read_nm):
    minsvlen = config.minsvlen_screen
    longinslen = config.long_ins_length / 2.0
    seq_cache_maxlen = config.dev_seq_cache_maxlen
    qname = read.query_name
    mapq = read.mapping_quality
    strand = "-" if read.is_reverse else "+"
    CINS = pysam.CINS
    CSOFT_CLIP = pysam.CSOFT_CLIP

    INS_SUM = 0
    # DEL_SUM = 0

    pos_read = 0
    pos_ref = read.reference_start
    for op, oplength in read.cigartuples:
        add_read, add_ref, event = OPLIST[op]
        if event:
            if op == CINS:
                INS_SUM += oplength
        pos_read += add_read * oplength
        pos_ref += add_ref * oplength

    return INS_SUM


def read_itersplits(read_id, read, contig, config, read_nm):
    # SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_signals = []
    supps = [part.split(",") for part in read.get_tag("SA").split(";") if len(part) > 0]
    trace_read = config.dev_trace_read != False and config.dev_trace_read == read.query_name

    if len(supps) > config.max_splits_base + config.max_splits_kb * (read.query_length / 1000.0):
        return

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [SignalProvider.read_itersplits] [{read.query_name}] passed max_splits check")

    if read.is_reverse:
        qry_start = read.query_length - read.query_alignment_end
    else:
        qry_start = read.query_alignment_start

    curr_signal = Signal(read_id,
                     read.query_name,
                     contig,
                     read.reference_start,
                     read.reference_start + read.reference_length,
                     qry_start,
                     qry_start + read.query_alignment_length,
                     "-" if read.is_reverse else "+",
                     read.mapping_quality,
                     read_nm,
                     "SPLIT_PRIM",
                     "?")
    all_signals.append(curr_signal)

    for refname, pos, strand, cigar, mapq, nm in supps:
        mapq = int(mapq)
        nm = int(nm)

        is_rev = (strand == "-")

        try:
            readstart_fwd, readstart_rev, refspan, readspan = CIGAR_analyze(cigar)
        except Exception as e:
            utils.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero = int(pos) - 1
        split_qry_start = readstart_rev if is_rev else readstart_fwd

        all_signals.append(Signal(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start + readspan,
                              strand,
                              mapq,
                              read_nm,
                              "SPLIT_SUP",
                              "?"))

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [SignalProvider.read_itersplits] [{read.query_name}] all_signals: {all_signals}")

    SV.classify_splits(read, all_signals, config, contig)

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [SignalProvider.read_itersplits] [{read.query_name}] classify_splits(all_signals): {all_signals}")

    for signal_i, signal in enumerate(all_signals):
        for svtype, svstart, arg in signal.svtypes_starts_lens:
            min_mapq = min(signal.mapq, all_signals[max(0, signal_i - 1)].mapq)
            if not config.dev_keep_lowqual_splits and min_mapq < config.mapq:
                continue
            if svtype != "NOSV":
                svlen = arg
                # print(type(svlen), type(svstart))
                yield Signal(read_id=signal.read_id,
                           read_qname=signal.read_qname,
                           contig=signal.contig,
                           ref_start=svstart,
                           ref_end=svstart,
                           qry_start=signal.qry_start,
                           qry_end=signal.qry_end,
                           strand=signal.strand,
                           mapq=signal.mapq,
                           nm=signal.nm,
                           source=signal.source,
                           svtype=svtype,
                           svlen=svlen,
                           seq=signal.seq if svtype == "INS" else None)



class SignalProvider:
    def __init__(self, config, read_id_offset):
        self.config = config

        self.signaltab = {}
        self.signalcounts = {}

        for svtype in SV.TYPES:
            self.signaltab[svtype] = {}
            self.signalcounts[svtype] = 0

        self.coverage_bins_fwd = {}
        self.coverage_bins_rev = {}
        self.coverage_bins_min_bin = None

        self.read_id = read_id_offset
        self.read_count = 0

        self.contig = None
        self.start = None
        self.end = None

    def record_signal(self, ld, pos_signaltab):
        signaltab_svtype = self.signaltab[ld.svtype]
        if pos_signaltab in signaltab_svtype:
            signaltab_svtype[pos_signaltab].append(ld)
            signal_count = len(signaltab_svtype[pos_signaltab])
            if signal_count > self.config.consensus_max_reads_bin:
                ld.seq = None
        else:
            signaltab_svtype[pos_signaltab] = [ld]
            signal_count = 1
        self.signalcounts[ld.svtype] += 1

    def build_signaltab(self, regions: list[Region], bam):

        assert (self.contig is None)
        assert (self.start is None)
        assert (self.end is None)

        externals = []
        ld_binsize = self.config.cluster_binsize

        for region in regions:
            self.contig = region.contig
            self.start = region.start if self.start is None else min(region.start, self.start)
            self.end = region.end if self.end is None else max(region.start, self.end)
            self.coverage_bins_min_bin = int(self.start / self.config.coverage_binsize) * self.config.coverage_binsize

            for ld in self.iter_region(bam, region):
                ld_contig, ld_ref_start = ld.contig, ld.ref_start

                if region.contig == ld_contig and region.start <= ld_ref_start < region.end:
                    pos_signaltab = int(ld_ref_start / ld_binsize) * ld_binsize
                    self.record_signal(ld, pos_signaltab)
                else:
                    externals.append(ld)

        return externals

    def iter_region(self, bam, region: Region):
        binsize = self.config.cluster_binsize
        coverage_binsize = self.config.coverage_binsize
        coverage_shift_bins = self.config.coverage_shift_bins
        coverage_shift_min_aln_len = self.config.coverage_shift_bins_min_aln_length
        long_ins_threshold = self.config.long_ins_length * 0.5
        qc_nm = self.config.qc_nm_measure
        advanced_tags = qc_nm
        mapq_min = self.config.mapq
        alen_min = self.config.min_alignment_length
        nm_sum = 0
        nm_count = 0
        trace_read = self.config.dev_trace_read

        for read in bam.fetch(region.contig, region.start, region.end, until_eof=False):
            if trace_read is not False:
                if trace_read == read.query_name:
                    print(f"[DEV_TRACE_READ] [0b/4] [SignalProvider.iter_region] [{region}] [{read.query_name}] has been fetched and is entering pre-filtering")

            if read.reference_start < region.start or read.reference_start >= region.end:
                continue

            self.read_id += 1
            self.read_count += 1

            alen = read.query_alignment_length
            # print(read.mapping_quality, mapq_min, alen, alen_min)
            if read.mapping_quality < mapq_min or read.is_secondary or alen < alen_min:
                continue

            has_sa = read.has_tag("SA")
            use_clips = self.config.detect_large_ins and not read.is_supplementary and not has_sa

            nm = -1
            curr_read_id = self.read_id
            if advanced_tags:
                if qc_nm:
                    if read.has_tag("NM"):
                        nm_raw = read.get_tag("NM")
                        nm_ratio = read.get_tag("NM") / float(read.query_alignment_length + 1)
                        ins_sum, del_sum = get_cigar_indels(curr_read_id, read, region.contig, self.config, use_clips, read_nm=nm)
                        nm_adj = (nm_raw - (ins_sum + del_sum))
                        nm_adj_ratio = nm_adj / float(read.query_alignment_length + 1)
                        nm = nm_adj_ratio
                        nm_sum += nm
                        nm_count += 1

            if trace_read is not False:
                if trace_read == read.query_name:
                    print(f"[DEV_TRACE_READ] [0b/4] [SignalProvider.iter_region] [{region}] [{read.query_name}] passed pre-filtering (whole-read), begin to extract signals")

            # Extract small indels
            for signal in read_iterindels(curr_read_id, read, region.contig, self.config, use_clips, read_nm=nm):
                if trace_read is not False:
                    if trace_read == read.query_name:
                        print(f"[DEV_TRACE_READ] [1/4] [preprocess.read_iterindels] [{region}] [{read.query_name}] new signal: {signal}")
                yield signal

            # Extract read splits
            if has_sa:
                if not read.is_supplementary:
                    if trace_read is not False:
                        if trace_read == read.query_name:
                            print(f"[DEV_TRACE_READ] [1/4] [preprocess.read_itersplits] [{region}] [{read.query_name}] is entering read_itersplits")
                    for signal in read_itersplits(curr_read_id, read, region.contig, self.config, read_nm=nm):
                        if trace_read is not False:
                            if trace_read == read.query_name:
                                print(f"[DEV_TRACE_READ] [1/4] [preprocess.read_itersplits] [{region}] [{read.query_name}] new signal: {signal}")
                        yield signal

            # Record in coverage table
            read_end = read.reference_start + read.reference_length
            assert (read_end == read.reference_end)
            # assert(read_end>=read.reference_start)
            if read.is_reverse:
                target_tab = self.coverage_bins_rev
            else:
                target_tab = self.coverage_bins_fwd
            covr_start_bin = (int(read.reference_start / coverage_binsize) + coverage_shift_bins * (alen >= coverage_shift_min_aln_len)) * coverage_binsize
            covr_end_bin = (int(read_end / coverage_binsize) - coverage_shift_bins * (alen >= coverage_shift_min_aln_len)) * coverage_binsize

            if covr_end_bin > covr_start_bin:
                self.coverage_bins_min_bin = min(self.coverage_bins_min_bin, covr_start_bin)
                target_tab[covr_start_bin] = target_tab[covr_start_bin] + 1 if covr_start_bin in target_tab else 1

                if read_end <= self.end:
                    target_tab[covr_end_bin] = target_tab[covr_end_bin] - 1 if covr_end_bin in target_tab else -1

        self.config.average_regional_nm = nm_sum / float(max(1, nm_count))
        self.config.qc_nm_threshold = self.config.average_regional_nm