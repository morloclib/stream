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


def scrambled(n, width, seed):
    """n records of `width` words with repeats, spanning the u64 range."""
    import random
    rng = random.Random(seed)
    pool = [rng.randrange(1 << 64) for _ in range(max(1, n // 3))] + [0, (1 << 64) - 1, 1 << 63, (1 << 63) - 1]
    return [rng.choice(pool) for _ in range(n * width)]


def sorted_records(xs, width):
    recs = sorted(tuple(xs[i : i + width]) for i in range(0, len(xs), width))
    return [x for r in recs for x in r]


def pieces(k, xs):
    return [xs[i : i + k] for i in range(0, len(xs), k)]
