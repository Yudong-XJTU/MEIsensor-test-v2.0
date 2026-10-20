#!/usr/bin/env python3
#
# MEIsensor post-calling annotation:
#   - EN motif: L1 endonuclease cleavage consensus (5'-TTTT/AA-3') at the insertion site
#   - SVAN: MEI annotation by calling SVAN (https://github.com/REPBIO-LAB/SVAN) on the MEIsensor VCF
#

import glob
import logging
import os
import shutil
import subprocess
import sys

import utils

log = logging.getLogger('MEIsensor.annotate')

#
# EN motif annotation
#
# L1 retrotransposition starts with a nick of the genomic target by the L1 ORF2p
# endonuclease (EN). The nick occurs at a T-rich/A-rich junction approaching the
# consensus 5'-TTTT/AA-3' on the nicked strand (Feng et al. 1996; Cost et al. 2002;
# Sen et al. 2007). Because the consensus on the opposite strand is its reverse
# complement, both nicks appear on the top strand as a T-tract followed by an
# A-tract; a long T-tract marks a plus-strand nick, a long A-tract a minus-strand nick.

EN_T_MIN = 3   # min. T-tract length on the nicked strand
EN_A_MIN = 2   # min. A-tract length on the nicked strand
EN_TRACT_MAX = 15  # cap for the reported tract lengths
EN_T_CANON = 4  # canonical 5'-TTTT/AA-3'
EN_A_CANON = 2


def _tract_len(seq, idx, base, backwards):
    """
    Length of the run of `base` starting at seq[idx], extending forward or backward.
    """
    n = 0
    while 0 <= idx < len(seq) and seq[idx] == base and n < EN_TRACT_MAX:
        n += 1
        idx += -1 if backwards else 1
    return n


def find_en_motif_in_seq(seq, ins_idx):
    """
    Best EN cleavage candidate in a flanking sequence window.

    Input:
        1. seq: upper-case reference sequence around the insertion point
        2. ins_idx: index in seq of the base preceding the insertion point; the
           insertion boundary is between seq[ins_idx] and seq[ins_idx+1]

    Output:
        (motif, strand, offset) of the best candidate, or None:
        - motif: T-tract/A-tract as found on the nicked strand, e.g. "TTTT/AA"
        - strand: strand carrying the nick ("+" or "-")
        - offset: signed distance (bp) between the cleavage junction and the
          insertion boundary; negative = junction upstream
    """
    best = None
    for j in range(len(seq) - 1):
        # top-strand T-tract ending at j, A-tract starting at j+1
        t = _tract_len(seq, j, 'T', backwards=True)
        a = _tract_len(seq, j + 1, 'A', backwards=False)
        if t < 1 or a < 1:
            continue

        candidates = []
        if t >= EN_T_MIN and a >= EN_A_MIN:
            # plus-strand nick: consensus read 5'->3' on the top strand
            candidates.append((abs(t - EN_T_CANON) + abs(a - EN_A_CANON), '+', 'T' * t, 'A' * a))
        if t >= EN_A_MIN and a >= EN_T_MIN:
            # minus-strand nick: bottom strand 5'-TTTT/AA-3' appears as top-strand TT/AAAA,
            # so the bottom-strand motif is the reverse complement (T-tract of length a)
            candidates.append((abs(a - EN_T_CANON) + abs(t - EN_A_CANON), '-', 'T' * a, 'A' * t))

        for canon_dist, strand, t_seq, a_seq in candidates:
            offset = j - ins_idx
            key = (abs(offset), canon_dist, 0 if strand == '+' else 1)
            if best is None or key < best[0]:
                best = (key, f"{t_seq}/{a_seq}", strand, offset)

    if best is None:
        return None
    return best[1], best[2], best[3]


def find_en_motif(reference_handle, contig, pos, window=25):
    """
    Annotate the L1 EN cleavage consensus near an insertion point.

    Input:
        1. reference_handle: open pysam.FastaFile
        2. contig: reference contig of the call
        3. pos: 1-based VCF POS of the base preceding the insertion
        4. window: search this many bp up/downstream of the insertion boundary

    Output:
        (motif, strand, offset) or None if no consensus-like T/A junction is found
    """
    if reference_handle is None:
        return None

    try:
        start = max(0, pos - 1 - window)
        seq = reference_handle.fetch(contig, start, pos + window).upper()
    except (KeyError, ValueError):
        return None

    if len(seq) == 0:
        return None

    ins_idx = pos - 1 - start  # seq index of the base preceding the insertion
    if not (0 <= ins_idx < len(seq) - 1):
        return None

    return find_en_motif_in_seq(seq, ins_idx)


#
# SVAN annotation
#
# SVAN (Structural Variants ANnotator, https://github.com/REPBIO-LAB/SVAN) annotates
# sequence-resolved insertions: mobile element family (FAM_N), strand, TSD, poly(A),
# conformation, transductions, nested MEIs, ... MEIsensor drives SVAN's own insertion
# annotation code (scripts/ins2fasta.py + SVAN-INS.py) unmodified:
#
#   1. prepare a SVAN input VCF (ALT holds the pure inserted sequence)
#   2. SVAN scripts/ins2fasta.py  -> insertions_seq.fa
#   3. trf insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs -> ins_trf.out
#   4. SVAN-INS.py <vcf> <trf.out> <vntr.bed> <exons.bed> <repeats.bed> <consensus.fa> <ref.fa> <out> -o <dir> -t <threads>

TRF_PARAMS = ["2", "7", "7", "80", "10", "10", "500", "-h", "-d", "-ngs"]


def resolve_svan_install(svan_path):
    """
    Locate SVAN-INS.py, scripts/ins2fasta.py and the GAPI module directory.

    Input:
        1. svan_path: SVAN repository root or the path of SVAN-INS.py itself

    Output:
        (svan_root, svan_ins_py, ins2fasta_py or None, gapi_dir or None)
    """
    svan_path = os.path.abspath(svan_path)
    if os.path.isdir(svan_path):
        svan_ins_py = os.path.join(svan_path, "SVAN-INS.py")
    else:
        svan_ins_py = svan_path
        svan_path = os.path.dirname(svan_path)

    if not os.path.isfile(svan_ins_py):
        utils.fatal_error_main(f"SVAN-INS.py not found at '{svan_ins_py}'. "
                               f"Pass the SVAN repository root (or SVAN-INS.py itself) via --svan-path.")

    ins2fasta_py = os.path.join(svan_path, "scripts", "ins2fasta.py")
    if not os.path.isfile(ins2fasta_py):
        ins2fasta_py = None

    gapi_dir = os.path.join(svan_path, "GAPI")
    if not os.path.isdir(gapi_dir):
        gapi_dir = None

    return svan_path, svan_ins_py, ins2fasta_py, gapi_dir


def resolve_svan_resources(config):
    """
    Resolve the VNTR/EXONS/REPEATS bed files and the consensus FASTA required by SVAN.
    Explicit --svan-* arguments win; otherwise files are searched in --svan-resources.

    Output:
        (vntr, exons, repeats, consensus)
    """
    res_dir = config.svan_resources

    def find_explicit(path, label):
        if path is None:
            return None
        if not os.path.isfile(path):
            utils.fatal_error_main(f"SVAN {label} file not found: '{path}'")
        return os.path.abspath(path)

    def find_in_dir(pattern, label):
        if res_dir is None:
            return None
        if not os.path.isdir(res_dir):
            utils.fatal_error_main(f"--svan-resources directory not found: '{res_dir}'")
        # ignore bwa/faidx index files that share the fasta's name
        index_exts = (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai", ".mmi")
        hits = sorted(h for h in glob.glob(os.path.join(res_dir, pattern))
                      if not h.endswith(index_exts))
        if len(hits) == 0:
            return None
        if len(hits) > 1:
            log.warning(f"Multiple candidate {label} files in '{res_dir}': "
                        f"{[os.path.basename(h) for h in hits]}; using '{hits[0]}'.")
        return os.path.abspath(hits[0])

    vntr = find_explicit(config.svan_vntr, "VNTR") or find_in_dir("*VNTR*.bed", "VNTR")
    exons = find_explicit(config.svan_exons, "EXONS") or find_in_dir("*EXONS*.bed", "EXONS")
    repeats = find_explicit(config.svan_repeats, "REPEATS") or find_in_dir("*REPEATS*.bed", "REPEATS")
    consensus = find_explicit(config.svan_consensus, "consensus") \
        or find_in_dir("*CONSENSUS*.fa*", "consensus") or find_in_dir("*consensus*.fa*", "consensus")

    missing = [label for label, path in (("VNTR bed", vntr), ("EXONS bed", exons),
                                         ("REPEATS bed", repeats), ("consensus FASTA", consensus))
               if path is None]
    if missing:
        utils.fatal_error_main(
            f"Missing SVAN resource files: {', '.join(missing)}. Provide them via --svan-vntr/--svan-exons/"
            f"--svan-repeats/--svan-consensus, or point --svan-resources at a directory containing them "
            f"(the hg38/chm13 bundles from https://zenodo.org/records/15229020 plus the retroelement "
            f"consensus FASTA).")

    return vntr, exons, repeats, consensus


def prepare_svan_input_vcf(vcf_path, out_path):
    """
    Convert a MEIsensor VCF into the sequence-resolved VCF SVAN expects:
    ALT is rewritten to hold only the inserted sequence (MEIsensor prepends
    the anchor reference base to ALT). Everything else is copied verbatim.

    Output:
        number of insertion records written
    """
    if vcf_path.endswith(".gz") or vcf_path.endswith(".bgz"):
        import pysam
        vcf_handle = pysam.BGZFile(vcf_path, "rb")
    else:
        vcf_handle = open(vcf_path, "r")

    n_records = 0
    with open(out_path, "w") as out:
        for line in vcf_handle:
            if isinstance(line, bytes):
                line = line.decode("utf-8")

            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            ref, alt = fields[3], fields[4]
            if alt.startswith("<"):
                utils.fatal_error_main(
                    "SVAN annotation needs sequence-resolved insertions in the ALT field. "
                    "Re-run MEIsensor without --symbolic to enable SVAN annotation.")

            # strip the leading anchor base(s) MEIsensor prepends to ALT
            if len(alt) > len(ref) and alt[:len(ref)] == ref:
                alt = alt[len(ref):]

            if len(alt) == 0:
                continue

            fields[4] = alt
            out.write("\t".join(fields) + "\n")
            n_records += 1

    vcf_handle.close()
    return n_records


def _write_insertions_fasta(svan_input_vcf, fasta_path):
    """
    Fallback for SVAN's scripts/ins2fasta.py: write one FASTA record per insertion,
    using the same id scheme (ID_CHROM:POS) SVAN-INS.py matches TRF hits against.
    """
    records = []
    with open(svan_input_vcf, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            records.append((f"{fields[2]}_{fields[0]}:{fields[1]}", fields[4]))

    with open(fasta_path, "w") as out:
        for seq_id, seq in records:
            out.write(f">{seq_id}\n{seq}\n")


def _run(command, env=None, stdout_path=None, cwd=None, label=""):
    """
    Run a subprocess, optionally redirecting stdout to a file; abort with the
    captured stderr tail on failure.
    """
    log.info(f"[SVAN] Running: {' '.join(str(c) for c in command)}")
    stderr_file = subprocess.PIPE if stdout_path else subprocess.PIPE
    stdout = open(stdout_path, "wb") if stdout_path else subprocess.PIPE
    try:
        proc = subprocess.run(command, env=env, cwd=cwd, stdout=stdout, stderr=stderr_file,
                              check=False, shell=False)
    except FileNotFoundError:
        utils.fatal_error_main(f"Could not execute {label or command[0]}. "
                               f"Make sure it is installed and on the PATH.")
    finally:
        if stdout_path:
            stdout.close()

    if proc.returncode != 0:
        stderr_tail = ""
        if proc.stderr is not None:
            stderr_tail = proc.stderr.decode("utf-8", errors="replace")[-2000:]
        utils.fatal_error_main(f"{label or command[0]} failed with exit code {proc.returncode}:\n{stderr_tail}")

    return proc


def _check_tools(config):
    missing = []
    for tool in [config.svan_trf, "bwa", "minimap2", "paftools.js"]:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        utils.fatal_error_main(
            f"SVAN annotation requires these executables on the PATH: trf, bwa, minimap2, paftools.js "
            f"(missing: {', '.join(missing)}). See https://github.com/REPBIO-LAB/SVAN for installation "
            f"notes; use --svan-trf to point at a non-standard trf location.")


def _read_vcf_lines(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        import pysam
        with pysam.BGZFile(path, "rb") as h:
            return [l.decode("utf-8") if isinstance(l, bytes) else l for l in h]
    with open(path, "r") as h:
        return h.readlines()


def _write_vcf_lines(path, lines):
    if path.endswith(".gz") or path.endswith(".bgz"):
        import pysam
        with pysam.BGZFile(path, "w") as h:
            for l in lines:
                h.write(l.encode() if isinstance(l, str) else l)
        try:
            pysam.tabix_index(path, preset="vcf", force=True)
        except Exception:
            log.exception(f"[SVAN] Re-indexing of {path} after merge failed.")
    else:
        with open(path, "w") as h:
            h.writelines(l if isinstance(l, str) else l.decode("utf-8") for l in lines)


def merge_svan_info(main_vcf_path, svan_vcf_path):
    """
    Merge SVAN's annotations back into the main MEIsensor VCF (rewritten in place).

    The INFO fields SVAN added (FAM_N, TSD_LEN, ...) are appended to each record's
    INFO column with an SVAN_ prefix, so SVAN's annotations never overwrite
    MEIsensor's own fields (both tools define e.g. STRAND -- SVAN's is kept as
    SVAN_STRAND). Fields SVAN merely passed through unchanged (SVTYPE, SUPPORT...)
    are not duplicated. Matching is by (CHROM, POS, ID); records SVAN dropped are
    left untouched. Header definitions for the merged fields are derived from
    SVAN's own header and inserted before #CHROM.

    Output:
        number of records that received SVAN annotations
    """
    # 1. SVAN output: INFO header definitions + per-record added INFO items
    svan_header_defs = {}  # field id -> raw ##INFO line
    svan_records = {}      # (chrom, pos, id) -> list of INFO items in output order
    for line in _read_vcf_lines(svan_vcf_path):
        if line.startswith("##INFO=<ID="):
            field_id = line.split("ID=", 1)[1].split(",", 1)[0]
            svan_header_defs[field_id] = line.rstrip("\n")
        elif line.startswith("#"):
            continue
        else:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            key = (fields[0], fields[1], fields[2])
            svan_records[key] = [it for it in fields[7].split(";") if it]

    # 2. rewrite the main VCF, appending the SVAN-added fields
    #    (records first: the merged field set is needed to build the header)
    main_lines = _read_vcf_lines(main_vcf_path)
    header_lines = []
    record_lines = []
    merged_header_fields = set()
    n_annotated = 0

    for line in main_lines:
        if line.startswith("#") or not line.strip():
            header_lines.append(line)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            record_lines.append(line)
            continue

        key = (fields[0], fields[1], fields[2])
        if key not in svan_records:
            record_lines.append(line)
            continue

        # SVAN passes the original INFO through; keep only fields it added or
        # changed (same name but different value, e.g. its own STRAND)
        existing = {}
        for it in fields[7].split(";"):
            if it:
                k, _, v = it.partition("=")
                existing[k] = v
        added = [it for it in svan_records[key]
                 if it.partition("=")[0] not in existing
                 or existing[it.partition("=")[0]] != it.partition("=")[2]]
        if not added:
            record_lines.append(line)
            continue

        fields[7] = fields[7] + ";" + ";".join("SVAN_" + it for it in added)
        merged_header_fields.update(it.split("=", 1)[0] for it in added)
        n_annotated += 1
        record_lines.append("\t".join(fields) + "\n")

    # header definitions for the merged fields, derived from SVAN's own ##INFO lines
    # (Number/Type/Description kept, ID renamed to SVAN_<ID>)
    if merged_header_fields:
        new_defs = [svan_header_defs[f].replace("ID=" + f, "ID=SVAN_" + f, 1) + "\n"
                    for f in sorted(merged_header_fields) if f in svan_header_defs]
        chrom_idx = next((i for i, ln in enumerate(header_lines) if ln.startswith("#CHROM")),
                         len(header_lines))
        header_lines[chrom_idx:chrom_idx] = new_defs

    _write_vcf_lines(main_vcf_path, header_lines + record_lines)
    return n_annotated


def run_svan_annotation(config):
    """
    Annotate the MEIsensor output VCF with SVAN by driving SVAN's own MEI code,
    then merge SVAN's annotations back into the MEIsensor VCF INFO (SVAN_* fields).
    SVAN's own annotated VCF is kept at config.svan_out.
    """
    if config.vcf is None:
        log.warning("[SVAN] No VCF output configured; skipping SVAN annotation.")
        return False

    if config.reference is None:
        utils.fatal_error_main("SVAN annotation requires --reference.")

    _check_tools(config)

    svan_root, svan_ins_py, ins2fasta_py, gapi_dir = resolve_svan_install(config.svan_path)
    vntr, exons, repeats, consensus = resolve_svan_resources(config)

    # bwa indexes for reference and consensus must be built before running SVAN
    for fasta, label in [(config.reference, "reference"), (consensus, "consensus")]:
        if not os.path.exists(fasta + ".bwt"):
            log.warning(f"[SVAN] No bwa index found next to the {label} FASTA ({fasta}); "
                        f"SVAN requires 'bwa index' to have been run on it beforehand.")

    out_dir = os.path.dirname(os.path.abspath(config.svan_out)) or "."
    os.makedirs(out_dir, exist_ok=True)
    work_dir = os.path.join(out_dir, os.path.splitext(os.path.basename(config.svan_out))[0] + ".svan_work")
    os.makedirs(work_dir, exist_ok=True)

    log.info(f"[SVAN] Annotating {config.vcf} -> {config.svan_out} (work dir: {work_dir})")

    # 1. sequence-resolved VCF for SVAN
    svan_input_vcf = os.path.join(work_dir, "meisensor.svan_input.vcf")
    n_records = prepare_svan_input_vcf(config.vcf, svan_input_vcf)
    if n_records == 0:
        log.info("[SVAN] No MEI calls to annotate; skipping SVAN.")
        return False
    log.info(f"[SVAN] Prepared {n_records} sequence-resolved insertions: {svan_input_vcf}")

    env = os.environ.copy()
    # SVAN-INS.py imports its modules flat ("import formats") while the GAPI modules
    # import each other package-style ("from GAPI import ..."), so both the GAPI
    # directory and the repository root have to be importable.
    path_entries = [gapi_dir, svan_root]
    env["PYTHONPATH"] = os.pathsep.join(p for p in path_entries if p is not None) \
        + os.pathsep + env.get("PYTHONPATH", "")
    if gapi_dir is None:
        log.warning("[SVAN] GAPI directory not found next to SVAN-INS.py; assuming the SVAN "
                    "python modules are already importable (e.g. via PYTHONPATH).")

    # 2. insertion sequences as FASTA (SVAN's own script when available)
    fasta_path = os.path.join(work_dir, "insertions_seq.fa")
    if ins2fasta_py is not None:
        _run([sys.executable, ins2fasta_py, svan_input_vcf, work_dir], env=env, cwd=work_dir,
             label="SVAN ins2fasta.py")
        fasta_path = os.path.join(work_dir, "insertions_seq.fa")
    else:
        log.warning("[SVAN] scripts/ins2fasta.py not found in the SVAN repository; writing "
                    "insertions_seq.fa with an equivalent built-in routine.")
        _write_insertions_fasta(svan_input_vcf, fasta_path)

    # 3. TRF on the insertion sequences
    trf_out = os.path.join(work_dir, "ins_trf.out")
    _run([os.path.abspath(shutil.which(config.svan_trf) or config.svan_trf), fasta_path] + TRF_PARAMS,
         env=env, stdout_path=trf_out, cwd=work_dir, label="TRF")

    # 4. SVAN-INS.py: annotate
    #    (SVAN appends '.vcf' to the fileName argument, so pass a stem)
    svan_out_base = os.path.basename(config.svan_out)
    svan_out_stem = svan_out_base[:-4] if svan_out_base.lower().endswith(".vcf") else svan_out_base
    _run([sys.executable, svan_ins_py,
          svan_input_vcf, trf_out, vntr, exons, repeats, consensus, config.reference,
          svan_out_stem,
          "-o", work_dir, "-t", str(config.threads)],
         env=env, cwd=work_dir, label="SVAN-INS.py")

    annotated = os.path.join(work_dir, svan_out_stem + ".vcf")
    if not os.path.exists(annotated):
        annotated = os.path.join(work_dir, svan_out_base)
    if not os.path.exists(annotated):
        # last resort: any VCF SVAN produced in the work dir, other than our input
        candidates = [p for p in glob.glob(os.path.join(work_dir, "*.vcf*"))
                      if os.path.basename(p) != os.path.basename(svan_input_vcf)]
        if candidates:
            annotated = max(candidates, key=os.path.getmtime)
    if not os.path.exists(annotated):
        utils.fatal_error_main(f"SVAN reported success but the annotated VCF was not found at "
                               f"'{os.path.join(work_dir, svan_out_stem + '.vcf')}'.")

    shutil.copyfile(annotated, config.svan_out)

    # 5. merge SVAN's annotations into the MEIsensor VCF INFO
    n_merged = merge_svan_info(config.vcf, config.svan_out)
    log.info(f"[SVAN] Merged SVAN annotations into {config.vcf} ({n_merged}/{n_records} records annotated)")

    if not config.svan_keep_tmp:
        svan_tmp = os.path.join(work_dir, "tmp")
        if os.path.isdir(svan_tmp):
            shutil.rmtree(svan_tmp, ignore_errors=True)

    log.info(f"[SVAN] Annotation done: {config.svan_out} (SVAN's own output; annotations also merged "
             f"into {config.vcf})")
    return True
