#include "bmh.h"
#include <cassert>
#include <cinttypes>
#include <chrono>

int main(int argc, char **argv) {
    if(argc > 3) {
        std::fprintf(stderr, "usage: %s <optional: n, # items> <optional: m, # sigs>\n", argv[0]);
        std::exit(EXIT_FAILURE);
    }
    size_t n = argc < 2 ? 1000: std::atoi(argv[1]);
    size_t m = argc < 3 ? 100: std::atoi(argv[2]);
    bmh_t<float> bmh(m), bmh2(m), bmh3(m), bmh4(m), bmh7(m);
    pmh1_t pmh1(m), pmh2(m);
    pmh1a pmh3(m), pmh4(m);
    assert(bmh.m() == m);
    assert(bmh2.m() == m);
    assert(std::equal(bmh.hvals_.data(), bmh.hvals_.data() + m, bmh2.hvals_.data()));
    auto start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < n; ++i) {
        if(i % 16 == 0) bmh7.update_1(i, 1);
        bmh.update_1(i, 1);
        bmh2.update_1(i, 1);
        bmh3.update_1(i, 1);
        bmh4.update_2(i, 1);
        //if(i % 100 == 0) std::fprintf(stderr, "Processed %zu/%zu\r\n", i, n);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "Updates for 5 BMH: %gs\n", static_cast<double>((stop - start).count() / 1000) * .001);
    start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < n; ++i) {
        pmh1.update(i, 1);
        pmh2.update(i, 1);
        pmh3.update(i, 1);
        pmh4.update(i, 1);
        //if(i % 100 == 0) std::fprintf(stderr, "Processed %zu/%zu\r\n", i, n);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::fprintf(stderr, "Updates for 2PMH: %gs\n", static_cast<double>((stop - start).count() / 1000) * .001);
    bmh4.finalize_2();
    auto s1 = bmh.to_sigs(), s2 = bmh2.to_sigs();
    auto s4 = bmh4.to_sigs();
    auto s7 = bmh7.to_sigs();
    auto ss1 = pmh1.to_sigs();
    auto ss2 = pmh2.to_sigs();
    pmh3.finalize();
    pmh4.finalize();
    auto ss3 = pmh3.to_sigs(), ss4 = pmh4.to_sigs();
    size_t nmatch = 0;
    for(size_t i = 0; i < m; ++i) {
        nmatch += s7[i] == s1[i];
        //std::fprintf(stderr, "%" PRIu64 ":%" PRIu64 ":%" PRIu64 "\n", s1[i], s2[i], s3[i]);
        //std::fprintf(stderr, "[%" PRIu64 ":%" PRIu64 ":%" PRIu64 "]\n", s4[i], s5[i], s6[i]);
    }
    std::fprintf(stderr, "Expected 1/16 signatures to match. Matching: %zu/%zu\n", nmatch, m);
    assert(std::equal(s1.begin(), s1.end(), s2.begin()));
    assert(std::equal(s1.begin(), s1.end(), s4.begin()));
    //assert(!std::equal(s1.begin(), s1.end(), s3.begin()));
}
