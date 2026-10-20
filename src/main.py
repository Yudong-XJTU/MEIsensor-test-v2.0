#!/usr/bin/env python3
#
# MEIsensor
import logging
import logging.config
import multiprocessing
import threading
import time
from collections import deque
from dataclasses import dataclass

import sys

if not sys.version_info >= (3, 10):
    print(f"Error: MEIsensor must be run with Python version 3.10 or above (detected Python version: {sys.version_info.major}.{sys.version_info.minor}). Exiting")
    exit(1)

import math
import os

import pysam

from config import Config
import vcf
import parallel_task
import utils


@dataclass
class MemoryUsage:
    """
    A snapshot of memory usage for this process and all its workers
    """
    pid: int
    mem: float  # RSS in KiB


class ResourceMonitor:
    _pid: int
    _parent: float
    _unknown: float
    _workers: dict[int, MemoryUsage]
    _filename = None

    WORKER_MEMORY_LIMIT = 1024 * 1024 * 2  # in KiB

    def __init__(self, config: Config):

        self._slots = config.threads
        self._update_interval = config.dev_monitor_memory or 30.0
        self._pid = os.getpid()

        try:
            from psutil import Process
        except ImportError:
            logging.getLogger('MEIsensor.dependencies').warning(f'Dependency psutil not available - resource monitoring is disabled.')
            self._process = None
            self._workers = {}
        else:
            self._process = Process(self._pid)

            self._running = True

            if config.dev_monitor_memory:
                self._filename = filename = f'memory-{config.run_id}.csv'
                try:
                    self._file = open(filename, "w")
                    self._file.write(self._generate_header())
                except:
                    logging.getLogger('MEIsensor.ResourceMonitor').exception(f'Unable to write memory usage to {filename}')
                else:
                    logging.getLogger('MEIsensor.ResourceMonitor').info(f'Logging memory usage to {filename}')
            else:
                self._file = None

            self.update()
            self._update_thread = threading.Thread(target=self._run_update, daemon=True)
            self._update_thread.start()

    @property
    def parent(self) -> float:
        return self._parent

    @property
    def workers(self) -> dict[int, float]:
        """
        dict pid -> memory
        """
        return {
            obj.pid: obj.mem for obj in self._workers.values() if obj is not None
        }

    @property
    def total(self) -> float:
        return self.parent + sum(self.workers.values())

    def stop(self):
        """
        Stop memory monitoring thread. This may take up to update_interval to complete.
        """
        self._running = False

    def _generate_header(self):
        return f'''total,pidp,memp,{','.join(f'pid{i},mem{i}' for i in range(self._slots))}\n'''

    def update(self):
        self._parent = self._process.memory_info().rss / 1024.0
        self._workers = {
            child.pid: MemoryUsage(child.pid, child.memory_info().rss / 1024.0) for child in self._process.children()
        }

        if self._file:
            self._file.write(f'''{','.join(self.generate())}\n''')
            self._file.flush()

        # print(f'Memory usage: {self.total / 1024:.2f} KiB')

    def _run_update(self):
        while self._running:
            self.update()
            time.sleep(self._update_interval)

        if self._file:
            self._file.close()

    def generate(self):
        yield str(self.total)
        yield str(self._pid)
        yield str(self.parent)
        for k, v in self.workers.items():
            yield str(k)
            yield str(v)

    @property
    def filename(self) -> str:
        return self._filename

    def __call__(self, worker_id: int, worker_pid: int) -> bool:
        """
        Call the resource monitor to get a recycle hint. Returns True if the given process should recycle.
        """
        if worker_pid in self._workers:
            return self._workers[worker_pid].mem > self.WORKER_MEMORY_LIMIT

        return False


def Main(processes: list[parallel_task.MEIsensorWorker]):
    # Determine MEIsensor run mode
    '''
    config = Config(
        # input='sample/HG00513/alignment/HG00513.hifi.minimap2.sorted.bam',
        input='../sample/test_alignment/HG00512.hifi.minimap2.sorted.matched_reads.bam',
        reference='../sample/reference/GRCh38.fasta',
        vcf='../sample/test_alignment/test_result/test.MEI.vcf',
    )
    '''
    config = Config()
    input_ext = config.input.split(".")[-1].lower()
    if input_ext not in ["bam"]:
        sys.exit(f"Error: Input file must be a .bam file. Provided file: {config.input}")

    # needed for running on osx
    if sys.platform == "darwin":
        multiprocessing.set_start_method("fork")

    if "bam" in input_ext or "cram" in input_ext:
        if input_ext.count("bam") + input_ext.count("cram") > 1:
            utils.fatal_error_main(f"Please specify max 1 .bam//.cram file as input (got {input_ext.count('bam')})")
        config.input = config.input

        if config.genotype_vcf is not None:
            config.mode = "genotype_vcf"
        else:
            config.mode = "call_sample"

        config.input_is_cram = False
        if "bam" in input_ext:
            config.input_mode = r"rb"
        elif "cram" in input_ext:
            config.input_mode = r"rc"
            config.input_is_cram = True

    log = logging.getLogger('MEIsensor.main')
    if config.dev_debug_log:
        logging.getLogger().setLevel(logging.DEBUG)
    if config.dev_progress_log:
        logging.getLogger('MEIsensor.progress').setLevel(logging.INFO)

    if config.mode == "call_sample":
        if config.sample_id is None:
            # config.sample_id,_=os.path.splitext(os.path.basename(config.input))
            config.sample_ids_vcf = [(0, "SAMPLE")]
        else:
            config.sample_ids_vcf = [(0, config.sample_id)]
    log.info(f"Run mode: {config.mode}")

    rkwargs = {}  # result kwargs

    monitor = ResourceMonitor(config)
    if monitor and monitor.filename is not None:
        logging.getLogger('MEIsensor.resources').info(f'Writing memory usage to {monitor.filename}')


    #
    # call_sample/genotype_vcf: Open .bam file for reading
    #
    if config.mode == "call_sample" or config.mode == "genotype_vcf":
        log.info(f"Input: {config.input}")
        bam_in = pysam.AlignmentFile(config.input, config.input_mode)
        try:
            has_index = bam_in.check_index()
            if not has_index:
                raise ValueError
        except ValueError:
            utils.fatal_error_main(f"Unable to load index for input file '{config.input}'. Please verify that your input file is sorted + indexed and that the index .bai file is valid and in the right location.")

    #
    # genotype_vcf: Read SVs from VCF to be genotyped
    #
    if config.mode == "genotype_vcf":
        path, ext = os.path.splitext(config.genotype_vcf)
        ext = ext.lower()
        if ext == ".gz":
            vcf_in_handle = pysam.BGZFile(config.genotype_vcf, "rb")
        elif ext == ".vcf":
            vcf_in_handle = open(config.genotype_vcf, "r")
        else:
            utils.fatal_error_main("Expected a .vcf or .vcf.gz file for genotyping using --genotype-vcf")
        vcf_in = vcf.VCF(config, vcf_in_handle)

        genotype_lineindex_order = []
        genotype_lineindex_svs = {}
        genotype_contig_svs = {}
        for svcall in vcf_in.read_svs_iter():
            if svcall.contig not in genotype_contig_svs:
                genotype_contig_svs[svcall.contig] = []
            assert (svcall.raw_vcf_line_index not in genotype_lineindex_svs)
            genotype_lineindex_order.append(svcall.raw_vcf_line_index)
            genotype_lineindex_svs[svcall.raw_vcf_line_index] = svcall
            genotype_contig_svs[svcall.contig].append(svcall)
        rkwargs['genotype_lineindex_order'] = genotype_lineindex_order
        log.info(f"Input: {config.genotype_vcf} ({len(genotype_lineindex_svs)} SVs to genotype)")

    #
    # Open output files
    #
    vcf_out = None
    if config.vcf is not None:

        vcf_output_info = []
        if config.sort:
            vcf_output_info.append("sorted")
        if config.vcf_output_bgz:
            vcf_output_info.append("bgzipped")
            vcf_output_info.append("tabix-indexed")

        if len(vcf_output_info) == 0:
            vcf_output_info_str = ""
        else:
            vcf_output_info_str = f"({', '.join(vcf_output_info)})"

        if os.path.exists(config.vcf) and not config.allow_overwrite:
            utils.fatal_error_main(f"Output file '{config.vcf}' already exists! Use --allow-overwrite to ignore this check and overwrite.")

        if config.vcf_output_bgz:
            if not config.sort:
                utils.fatal_error_main(".gz (bgzip) output is only supported with sorting enabled")
            vcf_handle = pysam.BGZFile(config.vcf, "w")
        else:
            vcf_handle = open(config.vcf, "w")

        vcf_out = vcf.VCF(config, vcf_handle)

        if config.mode == "call_sample":
            if config.reference is not None:
                log.info(f"Reference: {config.reference}")
            vcf_out.open_reference()

        log.info(f"Output: {config.vcf} {vcf_output_info_str}")
    #
    # Plan multiprocessing tasks
    #
    task_id = 0
    tasks = deque()
    contigs = []
    contig_tasks_intervals = {}

    if config.mode == "call_sample" or config.mode == "genotype_vcf":
        #
        # Process .bam header
        #
        task_classes = {
            'call_sample': parallel_task.CallTask,
            'genotype_vcf': parallel_task.GenotypeTask,
        }

        total_mapped = bam_in.mapped
        if (config.threads == 1 and not config.low_memory) or config.task_count_multiplier == 0:
            task_max_reads = total_mapped
        else:
            task_max_reads = max(1, math.floor(total_mapped / (config.threads * config.task_count_multiplier)))

        if total_mapped == 0:
            # Total mapped returns 0 for CRAM files
            config.task_read_id_offset_mult = 10 ** 9
        else:
            # BAM file
            config.task_read_id_offset_mult = 10 ** math.ceil(math.log(total_mapped) + 1)

        contig_lengths = []
        for contig in bam_in.get_index_statistics():
            if task_max_reads == 0:
                task_count = 1
            else:
                task_count = max(1, math.ceil(contig.mapped / float(task_max_reads)))
            contig_str = str(contig.contig)

            if config.contig and contig_str not in config.contig:
                continue

            if config.regions_by_contig and contig_str not in config.regions_by_contig:
                continue

            contigs.append(contig_str)
            contig_length = bam_in.get_reference_length(contig_str)
            contig_lengths.append((contig_str, contig_length))
            task_length = math.floor(contig_length / float(task_count))
            startpos = 0

            while startpos < contig_length - 1:
                endpos = min(contig_length - 1, startpos + task_length)
                if config.genotype_vcf is not None:
                    if contig_str in genotype_contig_svs:
                        genotype_svs = [target_sv for target_sv in genotype_contig_svs[contig_str] if target_sv.pos >= startpos and target_sv.pos < endpos]
                    else:
                        genotype_svs = []
                else:
                    genotype_svs = None

                task = task_classes[config.mode](
                    id=task_id,
                    contig=contig_str,
                    start=startpos,
                    end=endpos,
                    assigned_process_id=None,
                    genotype_svs=genotype_svs,
                    sv_id=0,
                    config=config,
                    regions=config.regions_by_contig.get(contig_str),
                )
                tasks.append(task)
                if contig_str not in contig_tasks_intervals:
                    contig_tasks_intervals[contig_str] = []
                contig_tasks_intervals[contig_str].append((task.start, task.end, task))
                startpos += task_length
                task_id += 1
        config.contig_lengths = contig_lengths

    if config.mode != "genotype_vcf" and config.vcf is not None:
        vcf_out.write_header(contig_lengths)
    elif config.mode == "genotype_vcf":
        vcf_out.rewrite_header_genotype(vcf_in.header_str)

    #
    # Start workers
    #
    if config.threads:
        for pnum in range(config.threads):
            processes.append(parallel_task.MEIsensorWorker(process_id=pnum, config=config, tasks=tasks, recycle_hint=monitor))
    else:
        processes.append(parallel_task.MEIsensorParentWorker(config=config, tasks=tasks))

    log.info("")
    if config.mode == "call_sample" or config.mode == "genotype_vcf":
        if config.input_is_cram:
            # CRAM file
            log.info(f"Scanning alignments (progress display disabled for CRAM input)")
        else:
            log.info(f"Scanning {total_mapped} aligned reads...")
    log.info("")

    #
    # Distribute analysis tasks to workers and collect results
    #
    analysis_start_time = time.monotonic()

    for p in processes:
        p.start()

    finished_tasks: list[parallel_task.Task] = []

    while any([p.run_parent() for p in processes if p.running]):
        time.sleep(0.1)

    for p in processes:
        p.finalize()
        finished_tasks.extend(p.finished_tasks)

    finished_tasks.sort(key=lambda task: task.id)

    for t in finished_tasks:
        t.result.emit(vcf_out=vcf_out, **rkwargs)

    if config.vcf is not None:
        vcf_out.close()
        if config.vcf_output_bgz:
            vcf_index_start_time = time.time()
            log.info(f"Indexing {config.vcf}...")
            try:
                pysam.tabix_index(config.vcf, preset="vcf", force=True)
            except:
                log.exception(f'Indexing failed.')
            else:
                log.info(f"Indexing done in {time.time() - vcf_index_start_time:.2f}s.")

    log.info(f"Analysis finished in {time.monotonic() - analysis_start_time:.2f}s.")
    log.info("")

    if config.mode == "call_sample" and config.vcf is not None:
        log.info(f"Wrote {vcf_out.call_count} MEI calls to {config.vcf} {vcf_output_info_str}")

    if monitor:
        log.debug(f'Stopping resource monitor.')
        monitor.stop()


if __name__ == "__main__":
    processes = []

    try:
        logging.config.dictConfig({
            'version': 1,
            'formatters': {
                'default': {
                    'format': '[%(asctime)s] %(levelname)-7s %(name)s: %(message)s',
                    'datefmt': '%Y-%m-%d %H:%M:%S'
                }
            },
            'handlers': {
                'console': {
                    'class': 'logging.StreamHandler',
                    'formatter': 'default',
                    'stream': 'ext://sys.stdout',
                }
            },
            'loggers': {
                'MEIsensor.progress': {
                    'level': logging.WARNING,
                },
                'MEIsensor.vcf': {
                    'level': logging.INFO,
                }
            },
            'root': {
                'level': logging.INFO,
                'handlers': ['console'],
            },
            'disable_existing_loggers': False,
        })
    except (ValueError, TypeError, AttributeError, ImportError):
        logging.exception(f'Error configuring loggers.')

    try:
        Main(processes)
    except (utils.MEIsensorExit, SystemExit) as exit_code:
        if len(processes):
            # Allow time for child process error messages to propagate
            print("[MEIsensor] Shutting down workers")
            time.sleep(10)
        for proc in processes:
            try:
                proc.process.terminate()
            except:
                pass

        for proc in processes:
            try:
                proc.process.join()
            except:
                pass
        exit(exit_code.code)
    except:
        logging.getLogger('MEIsensor.main').exception(f'Unhandled error while running MEIsensor.')
        exit(1)
