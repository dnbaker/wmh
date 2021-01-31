#ifndef FBMH_H__
#define FBMH_H__
#include "bmh.h"
namespace wmh {

template<typename BMH, typename FT, typename IT>
std::vector<uint64_t> minwise_det(const FT *weights, const IT *indices, size_t n, size_t m) {
    BMH h(m);
    for(size_t i = 0; i < n; ++i) {
        h.add(indices[i], weights[i]);
    }
    h.finalize();
    return h.to_sigs();
}

template<typename FT, typename IT>
std::vector<uint64_t> minhash(const FT *weights, const IT *indices, size_t n, size_t m, Sketchers sketcher, bool use_double=sizeof(FT) > 4) {
    std::vector<uint64_t> ret;
    switch(sketcher) {
        case S_BMH1: {
            if(use_double) ret = minwise_det<BagMinHash1<double>, FT, IT>(weights, indices, n, m);
            else ret = minwise_det<BagMinHash1<float>, FT, IT>(weights, indices, n, m);
        }
        break;

        case S_BMH2:
        if(use_double) ret = minwise_det<BagMinHash2<double>, FT, IT>(weights, indices, n, m);
        else ret = minwise_det<BagMinHash2<float>, FT, IT>(weights, indices, n, m);
        break;

        case S_PMH1:
        if(use_double) ret = minwise_det<pmh1_t<double>>(weights, indices, n, m); 
        else ret = minwise_det<pmh1_t<float>>(weights, indices, n, m); 
        break;
        case S_PMH1A: {
        ret = use_double ? minwise_det<pmh1a_t<double>>(weights, indices, n, m)
                         : minwise_det<pmh1a_t<float>>(weights, indices, n, m);
        }
        break;
        case S_PMH2: {
        ret = use_double ? minwise_det<pmh2_t<double>>(weights, indices, n, m)
                         : minwise_det<pmh2_t<float>>(weights, indices, n, m);
        }
        break;
        default: break;
    }
    return ret;
}
} // namespace wmh

#endif
