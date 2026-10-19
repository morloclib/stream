// Native ingest: the parser of fasta.hpp, with batches serialized as
// length-prefixed records and gathered into 16 MiB blocks; each block is cut
// into 2 MiB frames compressed with zstd on `threads` threads while parsing
// goes on. This is the work a morloc stream file write does (a 16 MiB write
// buffer per sub-packet; level 3 is zstd 6 in 2 MiB frames), without morloc.
//
//   native LEVEL THREADS INPUT OUTPUT
//
// LEVEL is a zstd level, 0 for none. Prints "entries<TAB>bases".

#include "fasta.hpp"

#include <zstd.h>

#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <mutex>
#include <thread>

constexpr size_t FRAME = size_t(2) << 20;
constexpr size_t BLOCK = size_t(16) << 20;

static void put_u64(std::string& out, uint64_t x) {
    out.append(reinterpret_cast<const char*>(&x), sizeof x);
}

static std::string serialize(const std::vector<BenchFa>& batch) {
    std::string out;
    size_t n = 8;
    for (const auto& e : batch) n += 16 + e.header.size() + e.bases.size();
    out.reserve(n);
    put_u64(out, batch.size());
    for (const auto& e : batch) {
        put_u64(out, e.header.size());
        out += e.header;
        put_u64(out, e.bases.size());
        out += e.bases;
    }
    return out;
}

static std::string compress_block(const std::string& raw, int level, size_t threads) {
    if (level == 0) return raw;
    const size_t nframes = (raw.size() + FRAME - 1) / FRAME;
    std::vector<std::string> frames(nframes);
    auto work = [&](size_t first) {
        ZSTD_CCtx* cx = ZSTD_createCCtx();
        for (size_t i = first; i < nframes; i += threads) {
            const size_t at = i * FRAME;
            const size_t len = std::min(FRAME, raw.size() - at);
            std::string& out = frames[i];
            out.resize(ZSTD_compressBound(len));
            size_t k = ZSTD_compressCCtx(cx, out.data(), out.size(), raw.data() + at, len, level);
            if (ZSTD_isError(k)) throw std::runtime_error(ZSTD_getErrorName(k));
            out.resize(k);
        }
        ZSTD_freeCCtx(cx);
    };
    std::vector<std::thread> ts;
    for (size_t t = 1; t < std::min(threads, nframes); ++t) ts.emplace_back(work, t);
    work(0);
    for (auto& t : ts) t.join();
    std::string out;
    put_u64(out, nframes);
    for (const auto& f : frames) {
        put_u64(out, f.size());
        out += f;
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: native LEVEL THREADS INPUT OUTPUT\n";
        return 2;
    }
    const int level = std::atoi(argv[1]);
    const size_t threads = std::max(1, std::atoi(argv[2]));
    std::FILE* out = std::fopen(argv[4], "wb");
    if (!out) {
        std::cerr << "cannot open " << argv[4] << "\n";
        return 1;
    }

    // One compressor thread behind the parser, holding at most two batches.
    std::mutex m;
    std::condition_variable cv;
    std::deque<std::string> queue;
    bool done = false;
    bool failed = false;
    std::thread writer([&]() {
        std::string pending;
        auto put = [&]() {
            std::string block = compress_block(pending, level, threads);
            pending.clear();
            if (std::fwrite(block.data(), 1, block.size(), out) != block.size()) {
                std::lock_guard<std::mutex> g(m);
                failed = true;
            }
        };
        for (;;) {
            std::string raw;
            {
                std::unique_lock<std::mutex> g(m);
                cv.wait(g, [&] { return !queue.empty() || done; });
                if (queue.empty()) break;
                raw = std::move(queue.front());
                queue.pop_front();
                cv.notify_all();
            }
            pending += raw;
            if (pending.size() >= BLOCK) put();
        }
        if (!pending.empty()) put();
    });

    uint64_t entries = 0, bases = 0;
    bench_parse_fasta(argv[3], [&](const std::vector<BenchFa>& batch) {
        entries += batch.size();
        for (const auto& e : batch) bases += e.bases.size();
        std::string raw = serialize(batch);
        std::unique_lock<std::mutex> g(m);
        cv.wait(g, [&] { return queue.size() < 2; });
        queue.push_back(std::move(raw));
        cv.notify_all();
    });
    {
        std::lock_guard<std::mutex> g(m);
        done = true;
        cv.notify_all();
    }
    writer.join();
    if (std::fclose(out) != 0 || failed) {
        std::cerr << "write error in " << argv[4] << "\n";
        return 1;
    }
    std::cout << entries << "\t" << bases << "\n";
    return 0;
}
