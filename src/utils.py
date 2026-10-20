#!/usr/bin/env python3
#
# util

import logging
import statistics
import sys
from dataclasses import dataclass
from typing import Optional


log = logging.getLogger(__name__)


@dataclass
class Region:
    """
    A region from a bed file.

    Coordinates are 0-based, half open
    """
    contig: str
    start: int
    end: int

    @classmethod
    def from_bed_line(cls, line: str) -> Optional['Region']:
        try:
            contig, start, end = line.split('\t')[:3]
            return Region(
                contig=contig,
                start=int(start),
                end=int(end),
            )
        except ValueError as ex:
            log.warning(f'Invalid region line: {line} ({ex})')
            return None

    def __str__(self) -> str:
        return f'{self.contig}:{self.start}-{self.end}'


class MEIsensorExit(Exception):
    code: int = 1


def stdev(nums):
    nums = list(nums)
    return statistics.stdev(nums) if len(nums) > 1 else 0


def median(nums):
    return int(statistics.median(nums))


def median_noavg(nums):
    nums = sorted(list(nums))
    mid = int(len(nums) / 2)
    return nums[mid]


def median_modes(nums):
    max_count = 0
    counts = {}
    for n in nums:
        if n not in counts:
            counts[n] = 1
        else:
            counts[n] += 1
        max_count = max(max_count, counts[n])
    return median_noavg(k for k, n in counts.items() if max_count - n < 3)


def mean(nums):
    nums = list(nums)
    return sum(nums) / len(nums)


def mean_or_none(nums):
    nums = list(nums)
    if len(nums) == 0:
        return None
    else:
        return sum(nums) / len(nums)


def mean_or_none_round(nums):
    r = mean_or_none(nums)
    if r is None:
        return r
    else:
        return round(r)


def trim(nums, pct=25):
    nums = sorted(list(nums))
    trim_n = int(len(nums) / float(100.0) * pct)
    if trim_n > 0:
        return nums[trim_n:-trim_n]
    else:
        return nums


def most_common(nums):
    counts = {}
    for n in nums:
        if n not in counts:
            counts[n] = 1
        else:
            counts[n] += 1
    return sorted(((counts[n], n) for n in counts), reverse=True)


def error(msg):
    sys.stderr.write("[MEIsensor] Error: " + msg + "\n")
    sys.stderr.flush()


def fatal_error(msg):
    error(msg + " (fatal, exiting)")
    exit(1)


def fatal_error_main(msg):
    error(msg + " (fatal, exiting)")
    raise MEIsensorExit


center = median_modes
