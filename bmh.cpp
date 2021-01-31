#include "bmh.h"
#include <cassert>
#include <cinttypes>

int main() {
    wd_t<float> wf;
    wd_t<double> wd;
    std::fprintf(stderr, "vals: %zu, %u\n", size_t(wd.ft2it()), int(wf.ft2it()));
    poisson_process_t<double, size_t> p(10, 1.3);
    size_t m = 20, n = 1000;
    bmh_t bmh(m), bmh2(m), bmh3(m);
    assert(bmh.m() == m);
    assert(bmh2.m() == m);
    assert(std::equal(bmh.hvals_.data(), bmh.hvals_.data() + m, bmh2.hvals_.data()));
    for(size_t i = 0; i < n; ++i) {
        bmh.update_1(i, 1);
        bmh2.update_1(i, 1);
        bmh3.update_1(i, 1);
    }
    auto s1 = bmh.to_sigs(), s2 = bmh2.to_sigs(), s3 = bmh3.to_sigs();
    for(size_t i = 0; i < m; ++i) {
        std::fprintf(stderr, "%" PRIu64 ":%" PRIu64 ":%" PRIu64 "\n", s1[i], s2[i], s3[i]);
    }
    assert(std::equal(s1.begin(), s1.end(), s2.begin()));
    //assert(!std::equal(s1.begin(), s1.end(), s3.begin()));
}
