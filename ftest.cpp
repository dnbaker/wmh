#include "fbmh.h"
#include <random>
#include <chrono>

#define _T()  std::chrono::high_resolution_clock::now()
int main() {
    size_t _M = 100, _N = 1000000;
    if(char *s = std::getenv("_M")) _M = std::strtoull(s, nullptr, 10);
    if(char *s = std::getenv("_N")) _N = std::strtoull(s, nullptr, 10);
    std::vector<uint64_t> data(_N);
    std::vector<double> weights(_N);
    std::poisson_distribution<size_t> dst(50);
    wy::WyRand<uint64_t> rng;
    for(auto &i: weights) i = dst(rng);
    for(const auto s: {S_PMH1, S_PMH1A, S_PMH2, S_BMH2, S_BMH1}) {
        for(size_t i = 0; i < 10; ++i) {
            size_t to_use = (i + 1) * (_N / 10);
            auto t = _T();
            auto samples = wmh::minhash(weights.data(), data.data(), to_use, _M, s);
            auto tt = _T();
            std::fprintf(stderr, "%s\t%zu\t%g\n", wmh::mh2str(s),
                         to_use, std::chrono::duration<double, std::milli>(tt - t).count());
        }
    }
}
