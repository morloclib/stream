#!/usr/bin/env python3
"""Time morloc's stream ingest (stream-cpp) against native C++.

    python3 bench.py [REPS]

Builds both programs if missing, writes the inputs to data/ on first use,
runs every configuration REPS times (default 3) and appends the medians to
results.tsv. Both programs parse with fasta.hpp; native.cpp writes the same
batches as zstd frames without morloc. Each morloc output is read back and
its entry and base counts checked against the native program's.
"""

import datetime
import os
import random
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "out")
ZSTD = os.environ.get("ZSTD_PREFIX", "/env/.pixi/envs/default")
THREADS = os.cpu_count() or 1

# name -> (entries, bases per entry, bases per line)
INPUTS = {
    "reads": (1_000_000, 150, 150),
    "genome": (25, 4_000_000, 60),
    "tiny": (1, 10, 60),
}


def write_fasta(path, entries, length, width):
    rng = random.Random(1)
    with open(path, "w") as f:
        for i in range(entries):
            f.write(f">seq{i} synthetic\n")
            s = "".join(rng.choices("ACGT", k=length))
            for j in range(0, length, width):
                f.write(s[j : j + width])
                f.write("\n")


def build():
    os.chdir(HERE)
    if not os.path.exists("ingest-native"):
        subprocess.run(
            ["g++", "-O2", "-std=c++17", f"-I{ZSTD}/include", "native.cpp", "-o", "ingest-native",
             f"-L{ZSTD}/lib", f"-Wl,-rpath,{ZSTD}/lib", "-lzstd", "-pthread"],
            check=True)
    if not os.path.exists("ingest-morloc"):
        subprocess.run(["morloc", "make", "-o", "ingest-morloc", "main.loc"], check=True)


def run(cmd):
    """Wall seconds, CPU seconds, largest process's max RSS (MiB), stdout.

    The rusage covers the command and the descendants it waited for."""
    with tempfile.TemporaryFile("w+") as out, tempfile.TemporaryFile("w+") as err:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=out, stderr=err)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.perf_counter() - t0
        p.returncode = os.waitstatus_to_exitcode(status)
        out.seek(0)
        err.seek(0)
        if p.returncode:
            raise SystemExit(f"{' '.join(cmd)} failed: {err.read()}")
        return wall, ru.ru_utime + ru.ru_stime, ru.ru_maxrss / 1024, out.read()


def counts_of(stdout):
    return tuple(int(x) for x in stdout.replace("[", " ").replace("]", " ").replace(",", " ").split())


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    build()
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    configs = [
        ("native", "0", "1"),
        ("native", "6", str(THREADS)),
        ("morloc", "0", None),
        ("morloc", "3", None),
    ]
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    rows = []
    for name, (entries, length, width) in INPUTS.items():
        src = os.path.join(DATA, name + ".fa")
        if not os.path.exists(src):
            write_fasta(src, entries, length, width)
        mb = os.path.getsize(src) / 2**20
        expected = None
        for prog, level, threads in configs:
            dest = os.path.join(OUT, f"{name}-{prog}-{level}")
            walls, cpus, rss = [], [], []
            for _ in range(reps):
                if os.path.exists(dest):
                    os.remove(dest)
                if prog == "native":
                    w, c, r, out = run(["./ingest-native", level, threads, src, dest])
                    expected = counts_of(out)
                else:
                    w, c, r, _ = run(["./ingest-morloc", "ingestFa", "--compression-level", level, src, dest])
                walls.append(w)
                cpus.append(c)
                rss.append(r)
            if prog == "morloc":
                got = counts_of(run(["./ingest-morloc", "countFa", dest])[3])
                if got != expected:
                    raise SystemExit(f"{name} level {level}: morloc wrote {got}, native read {expected}")
            label = f"{prog} level {level}" + (f" threads {threads}" if threads else "")
            wall = statistics.median(walls)
            rows.append((stamp, name, f"{mb:.1f}", label, f"{wall:.3f}", f"{statistics.median(cpus):.3f}",
                         f"{max(rss):.0f}", f"{mb / wall:.0f}", f"{os.path.getsize(dest) / 2**20:.1f}"))
            print("\t".join(rows[-1][1:]), flush=True)
    path = os.path.join(HERE, "results.tsv")
    new = not os.path.exists(path)
    with open(path, "a") as f:
        if new:
            f.write("date\tinput\tinput_mib\tprogram\twall_s\tcpu_s\tmax_rss_mib\tmib_per_s\toutput_mib\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
