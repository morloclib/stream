# File helpers for the stream test suite.
import os
import tempfile

_dir = None


def tmp_path(name):
    global _dir
    if _dir is None:
        _dir = tempfile.mkdtemp(prefix="stream-test-")
    return os.path.join(_dir, name)


def truncate_file(path, nbytes):
    size = os.path.getsize(path)
    with open(path, "r+b") as f:
        f.truncate(max(0, size - nbytes))


def tmp_dir(name):
    d = tmp_path(name)
    os.makedirs(d, exist_ok=True)
    return d


def files_in(d):
    return len(os.listdir(d))


def scrambled(n, width, bits, seed):
    """n records of `width` words with repeats, spanning [0, 2^bits)."""
    import random
    rng = random.Random(seed)
    top = 1 << bits
    edges = [0, 1, (1 << 32) - 1, 1 << 32, top - 1] + ([1 << 63, (1 << 63) - 1] if bits > 63 else [])
    pool = [rng.randrange(top) for _ in range(max(1, n // 3))] + edges
    return [rng.choice(pool) for _ in range(n * width)]


def sorted_records(xs, width):
    recs = sorted(tuple(xs[i : i + width]) for i in range(0, len(xs), width))
    return [x for r in recs for x in r]


def pieces(k, xs):
    return [xs[i : i + k] for i in range(0, len(xs), k)]
