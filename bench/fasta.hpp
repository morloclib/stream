#ifndef __STREAM_BENCH_FASTA_HPP__
#define __STREAM_BENCH_FASTA_HPP__

// The FASTA parser both ingest benchmarks share: the morloc reader below and
// the native program in native.cpp.

#include <cstdio>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct BenchFa {
    std::string header;
    std::string bases;
};

// A batch is handed on at this many entries or bases, whichever comes first.
constexpr size_t BENCH_BATCH_ENTRIES = 4096;
constexpr size_t BENCH_BATCH_BASES = size_t(8) << 20;

// Reads a FASTA file one batch at a time.
class BenchFaReader {
  public:
    explicit BenchFaReader(const std::string& path) : path_(path), buf_(size_t(1) << 20) {
        f_ = std::fopen(path.c_str(), "rb");
        if (!f_) throw std::runtime_error("cannot open " + path);
    }
    BenchFaReader(const BenchFaReader&) = delete;
    BenchFaReader& operator=(const BenchFaReader&) = delete;
    ~BenchFaReader() {
        if (f_) std::fclose(f_);
    }

    // The next batch; empty at the end of the file.
    std::vector<BenchFa> next() {
        std::vector<BenchFa> batch;
        if (pending_) {
            batch.push_back(std::move(*pending_));
            pending_.reset();
        }
        size_t bases = batch.empty() ? 0 : batch.back().bases.size();
        std::string line;
        while (line_(line)) {
            const char* s = line.data();
            size_t n = line.size();
            if (n > 0 && s[n - 1] == '\r') --n;
            if (n == 0) continue;
            if (s[0] == '>') {
                BenchFa e{std::string(s + 1, n - 1), std::string()};
                if (!batch.empty() && (batch.size() >= BENCH_BATCH_ENTRIES || bases >= BENCH_BATCH_BASES)) {
                    pending_ = std::move(e);
                    return batch;
                }
                batch.push_back(std::move(e));
            } else {
                if (batch.empty()) throw std::runtime_error(path_ + ": sequence before the first '>' header");
                batch.back().bases.append(s, n);
                bases += n;
            }
        }
        return batch;
    }

  private:
    // The next line without its newline; false at the end of the file.
    bool line_(std::string& out) {
        out.clear();
        for (;;) {
            if (at_ == len_) {
                if (eof_) return !out.empty();
                len_ = std::fread(buf_.data(), 1, buf_.size(), f_);
                at_ = 0;
                if (len_ == 0) {
                    if (std::ferror(f_)) throw std::runtime_error("read error in " + path_);
                    eof_ = true;
                    return !out.empty();
                }
            }
            const char* p = buf_.data() + at_;
            const char* nl = static_cast<const char*>(std::memchr(p, '\n', len_ - at_));
            if (!nl) {
                out.append(p, len_ - at_);
                at_ = len_;
                continue;
            }
            out.append(p, nl - p);
            at_ = (nl - buf_.data()) + 1;
            return true;
        }
    }

    std::string path_;
    std::FILE* f_ = nullptr;
    std::vector<char> buf_;
    size_t at_ = 0, len_ = 0;
    bool eof_ = false;
    std::optional<BenchFa> pending_;
};

// Calls `emit(batch)` for each batch of entries in the file, in order.
template <class Emit>
void bench_parse_fasta(const std::string& path, Emit emit) {
    BenchFaReader r(path);
    for (auto batch = r.next(); !batch.empty(); batch = r.next()) emit(batch);
}

template <class Sink>
void bench_read_fasta(const std::string& path, Sink sink) {
    bench_parse_fasta(path, [&](const std::vector<BenchFa>& batch) { sink(batch); });
}

inline std::tuple<int, int> bench_count_batch(const std::tuple<int, int>& acc, const std::vector<BenchFa>& batch) {
    int bases = std::get<1>(acc);
    for (const auto& e : batch) bases += static_cast<int>(e.bases.size());
    return {std::get<0>(acc) + static_cast<int>(batch.size()), bases};
}

#endif
