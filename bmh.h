#ifndef BAGMINHASH_H__
#define BAGMINHASH_H__
#include <stdexcept>
#include <cassert>
#include "aesctr/wy.h"
#include <queue>
#include <div.h>

template<typename FT>
struct mvt_t {
    std::vector<FT> data_;
    mvt_t(size_t m, const FT maxv=std::numeric_limits<FT>::max()): data_((m << 1) - 1, maxv)
    {
        assert(getm() == m);
    }


    FT *data() {return data_.data();}
    const FT *data() const {return data_.data();}
    // Check size and max
    size_t getm() const {return (data_.size() >> 1) + 1;}
    FT max() const {return data_.back();}
    FT operator[](size_t i) const {return data_[i];}
    void reset() {
        std::fill(data_.begin(), data_.end(), std::numeric_limits<FT>::max());
    }

    bool update(size_t index, FT x) {
        const auto sz = data_.size();
        const auto mv = getm();
        if(x < data_[index]) {
            do {
                data_[index] = x;
                index = mv + (index >> 1);
                if(index >= sz) break;
                size_t lhi = (index - mv) << 1;
                size_t rhi = lhi + 1;
                x = std::max(data_[lhi], data_[rhi]);
            } while(x < data_[index]);
            return true;
        }
        return false;
    }
};

template <typename T>
class MaxValueTracker {
    const uint32_t m;
    std::vector<T> values;

public:
    MaxValueTracker(uint32_t _m, const T infinity=std::numeric_limits<T>::max()) : m(_m), values((_m << 1) - 1, infinity) {}

    uint64_t getm() const {return this->m;}
    T *data() {return values.data();}
    const T *data() const {return values.data();}
    void update(uint32_t idx, T value) {
        assert(idx < m);
        while(value < values[idx]) {
            values[idx] = value;
            idx = m + (idx >> 1);
            if (idx >= values.size()) break;
            uint32_t leftChildIdx = (idx - m) << 1;
            uint32_t rightChildIdx = leftChildIdx + 1;
            value = std::max(values[leftChildIdx], values[rightChildIdx]);
        }
    }

    const T& max() const {
        return values.back();
    }

    const T& operator[](uint32_t idx) const {
        return values[idx];
    }
};

template<typename FT>
using DefIT = std::conditional_t<sizeof(FT) == 4, uint32_t,
               std::conditional_t<sizeof(FT) == 8, uint64_t,
               std::conditional_t<sizeof(FT) == 2, uint16_t,
               std::conditional_t<sizeof(FT) == 1, uint8_t,
               void>>>>;

template<typename FT>
struct wd_t {
    using IT = DefIT<FT>;
    static_assert(std::is_integral_v<IT>, "Sanity check");

    union ITFTU {
        IT i_; FT f_;
        ITFTU(): i_(0) {}
    };
    static constexpr IT ft2it(FT val=std::numeric_limits<FT>::max()) {
        ITFTU tmp;
        tmp.f_ = val;
        return tmp.i_;
    }
    static constexpr FT it2ft(IT val) {
        ITFTU tmp;
        tmp.i_ = val;
        return tmp.f_;
    }
    template<typename OIT, typename=std::enable_if_t<std::is_integral_v<OIT>>>
    static constexpr FT cvt(OIT val) {return it2ft(val);}
    template<typename OFT, typename=std::enable_if_t<std::is_floating_point_v<OFT>>>
    static constexpr IT cvt(OFT val) {return ft2it(val);}
    static constexpr FT maxv = std::numeric_limits<FT>::max();
    static const IT maxi = cvt(maxv);
    using IntType = IT;
};

template<typename FT, typename IT=DefIT<FT>>
struct poisson_process_t {
    static_assert(std::is_arithmetic_v<FT>, "Must be arithmetic");
    static_assert(std::is_integral_v<IT>, "Must be intgral");
    // Algorithm 4
    FT x_, weight_, minp_, maxq_;
    IT idx_ = std::numeric_limits<IT>::max();
    uint64_t wyv_; // RNG state
    using wd = wd_t<FT>;
public:
    poisson_process_t(FT x, FT w, FT p, FT q, uint64_t seed):
        x_(x), weight_(w), minp_(p), maxq_(q),
        wyv_(seed)
        {
        assert(minp_ < maxq_);
    }
    poisson_process_t& operator=(poisson_process_t &&) = default;
    poisson_process_t& operator=(const poisson_process_t &) = default;
    poisson_process_t(const poisson_process_t &o) = default;
    poisson_process_t(poisson_process_t &&o) = default;
    poisson_process_t(IT id, FT w): x_(0.), weight_(w), minp_(0.), maxq_(std::numeric_limits<FT>::max()), wyv_(id) {

    }
    IT widxmax() const {
        return wd::cvt(maxq_);
    }
    IT widxmin() const {
        return wd::cvt(minp_);
    }
    bool partially_relevant() const {
        return wd::cvt(widxmin() + 1) <= weight_;
    }
    bool fully_relevant() const {
        return wd::cvt(widxmax()) <= weight_;
    }
    bool can_split() const {
        return widxmax() > widxmin() + 1;
    }
    // Note: > is reversed, for use in pq
    bool operator>(const poisson_process_t &o) const {return x_ < o.x_;}
    bool operator<(const poisson_process_t &o) const {return x_ > o.x_;}
    template<typename OIT>
    void step(const schism::Schismatic<OIT> &fastmod) {
        // Top 52-bits as U01 for exponential with weight of q - p,
        // bottom logm bits for index
        uint64_t xi = wy::wyhash64_stateless(&wyv_);
        x_ += -std::log(static_cast<double>(xi >> 12) * 0x1p-52) / (maxq_ - minp_);
        assert(fastmod.mod(static_cast<IT>(xi)) == (static_cast<IT>(xi) % fastmod.d()) || !std::fprintf(stderr, "lhs: %zu. rhs: %zu. xi: %u. d(): %u\n", size_t(fastmod.mod(static_cast<IT>(xi))), size_t(static_cast<IT>(xi) % fastmod.d()), static_cast<IT>(xi), int(fastmod.d())));
        idx_ = fastmod.mod(xi);
    }
    void step(size_t m) {
        uint64_t xi = wy::wyhash64_stateless(&wyv_);
        x_ += -std::log(static_cast<double>(xi >> 12) * 0x1p-52) / (maxq_ - minp_);
        idx_ = xi % m;
    }
    poisson_process_t split() {
        uint64_t midpoint = (uint64_t(widxmin()) + widxmax()) / 2;
        double midval = wd::cvt(midpoint);
        uint64_t xval = wd::cvt(x_) ^ wy::wyhash64_stateless(&midpoint);
        const double p = (midval - minp_) / (maxq_ - minp_);
        const double rv = static_cast<double>((wy::wyhash64_stateless(&xval) >> 12) * 0x1p-52);
        if(rv < p) {
            auto oldmaxq = maxq_;
            maxq_ = midval;
            return poisson_process_t(x_, weight_, midval, oldmaxq, xval);
        } else {
            auto oldminp = minp_;
            minp_ = midval;
            return poisson_process_t(x_, weight_, oldminp, midval, xval);
        }
    }
};


template<typename FT>
static inline uint64_t reg2sig(FT v) {
    uint64_t t = 0;
    std::memcpy(&t, &v, std::min(sizeof(v), sizeof(uint64_t)));
    t ^= 0xcb1eb4b41a93fe67uLL;
    return wy::wyhash64_stateless(&t);
}

template<typename FT=float>
struct bmh_t {
    using wd = wd_t<FT>;
    using IT = typename wd::IntType;
    using PoissonP = poisson_process_t<FT, IT>;

    struct pq_t: public std::priority_queue<PoissonP, std::vector<PoissonP>> {
        auto &getc() {return this->c;}
        const auto &getc() const {return this->c;}
        void clear() {
            this->c.clear();
            assert(this->size() == 0);
        }
    };
    pq_t heap_;
    mvt_t<FT> hvals_;
    schism::Schismatic<IT> div_;
    auto m() const {return hvals_.getm();}

    bmh_t(size_t m): hvals_(m), div_(m) {
        heap_.getc().reserve(m);
    }
    void update_2(IT id, FT w) {
        auto &tmp = heap_.getc();
        if(w <= 0.) return;
        PoissonP p(id, w);
        p.step(div_);
        if(p.fully_relevant()) hvals_.update(p.idx_, p.x_);
        const size_t offset = tmp.size();
        size_t mainiternum = 0, subin  =0;
        while(p.x_ < hvals_.max()) {
            ++mainiternum;
            //std::fprintf(stderr, "x: %g. max: %g\n", p.x_, hvals_.max());
            while(p.can_split() && p.partially_relevant()) {
                ++subin;
                //std::fprintf(stderr, "min %g max %g, splitting!\n", p.minp_, p.maxq_);
                auto pp = p.split();
                if(p.fully_relevant())
                    hvals_.update(p.idx_, p.x_);
                if(pp.partially_relevant()) {
                    pp.step(div_);
                    if(pp.fully_relevant()) hvals_.update(pp.idx_, pp.x_);
                    if(pp.partially_relevant()) heap_.push(std::move(pp));
                    tmp.emplace_back(std::move(pp));
                    std::push_heap(tmp.begin() + offset, tmp.end());
                }
                //std::fprintf(stderr, "Finishing subloop at %zu/%zu\n", mainiternum, subin);
            }
            if(p.fully_relevant()) {
                p.step(div_);
                hvals_.update(p.idx_, p.x_);
                if(p.x_ <= hvals_.max()) {
                    tmp.emplace_back(std::move(p));
                    std::push_heap(tmp.begin() + offset, tmp.end());
                }
            }
            if(tmp.size() == offset) break;
            std::pop_heap(tmp.begin() + offset, tmp.end());
            p = std::move(tmp.back());
            tmp.pop_back();
        }
        auto bit = tmp.begin() + offset;
        for(;bit != tmp.begin() && tmp.front().x_ > hvals_.max();--bit) {
            std::pop_heap(tmp.begin(), bit, std::greater<>());
        }
        for(auto hit = tmp.begin() + offset; hit != tmp.end(); ++hit) {
            if(hit->x_ <= hvals_.max()) {
                *bit++ = std::move(*hit);
                std::push_heap(tmp.begin(), bit, std::greater<>());
            }
        }
        tmp.erase(bit, tmp.end());
    }
    void finalize_2() {
        auto &tmp = heap_.getc();
        std::make_heap(tmp.begin(), tmp.end());
        while(tmp.size()) {
            std::pop_heap(tmp.begin(), tmp.end());
            auto p = std::move(tmp.back());
            tmp.pop_back();
            if(p.x_ > hvals_.max()) break;
            while(p.can_split() && p.partially_relevant()) {
                auto pp = p.split();
                if(p.fully_relevant()) hvals_.update(p.idx_, p.x_);
                if(pp.partially_relevant()) {
                    pp.step(div_);
                    if(pp.fully_relevant()) hvals_.update(pp.idx_, pp.x_);
                    if(pp.x_ <= hvals_.max()) {
                        tmp.emplace_back(std::move(pp));
                        std::push_heap(tmp.begin(), tmp.end());
                    }
                }
            }
            if(p.fully_relevant()) {
                p.step(div_);
                hvals_.update(p.idx_, p.x_);
                if(p.x_ <= hvals_.max()) {
                    tmp.emplace_back(std::move(p));
                    std::push_heap(tmp.begin(), tmp.end());
                }
            }
        }
    }
    void update_1(IT id, FT w) {
        if(w <= 0.) return;
        PoissonP p(id, w);
        p.step(div_);
        if(p.fully_relevant()) hvals_.update(p.idx_, p.x_);
        size_t mainiternum = 0, subin  =0;
        while(p.x_ < hvals_.max()) {
            ++mainiternum;
            //std::fprintf(stderr, "x: %g. max: %g\n", p.x_, hvals_.max());
            while(p.can_split() && p.partially_relevant()) {
                ++subin;
                //std::fprintf(stderr, "min %g max %g, splitting!\n", p.minp_, p.maxq_);
                auto pp = p.split();
                if(p.fully_relevant())
                    hvals_.update(p.idx_, p.x_);
                if(pp.partially_relevant()) {
                    pp.step(div_);
                    if(pp.fully_relevant()) hvals_.update(pp.idx_, pp.x_);
                    if(pp.partially_relevant()) heap_.push(std::move(pp));
                }
                //std::fprintf(stderr, "Finishing subloop at %zu/%zu\n", mainiternum, subin);
            }
            if(p.fully_relevant()) {
                p.step(div_);
                hvals_.update(p.idx_, p.x_);
                if(p.x_ <= hvals_.max()) heap_.push(std::move(p));
            }
            if(heap_.empty()) break;
            p = std::move(heap_.top());
            heap_.pop();
            //std::fprintf(stderr, "heap size: %zu\n", heap_.size());
        }
        heap_.clear();
    }
    template<typename IT=uint64_t>
    std::vector<IT> to_sigs() const {
        std::vector<IT> ret(m());
        std::transform(hvals_.data(), hvals_.data() + m(), ret.begin(), reg2sig<FT>);
        return ret;
    }
};


template<typename FT=float>
struct pmh1_t {
    using wd = wd_t<FT>;
    using IT = typename wd::IntType;

    mvt_t<FT> hvals_;
    schism::Schismatic<IT> div_;
    std::vector<FT> res_;
    pmh1_t(size_t m): hvals_(m), div_(m), res_(m) {}

    void update(const IT id, const FT w) {
        if(w <= 0.) return;
        const FT wi = 1. / w;
        uint64_t hi = id;
        uint64_t xi = wy::wyhash64_stateless(&hi);
        auto hv = -std::log((xi >> 12) * 0x1p-52) * wi;
        while(hv < hvals_.max()) {
            //if(++iternum % 100 == 0) std::fprintf(stderr, "Inner loop %zu with hv = %g and max = %g\n", iternum, hv, hvals_.max());
            auto idx = div_.mod(xi);
            if(hvals_.update(idx, hv)) {
                res_[idx] = id;
                if(hv >= hvals_.max()) break;
            }
            xi = wy::wyhash64_stateless(&hi);
            hv += -std::log((xi >> 12) * 0x1p-52) * wi;
        }
    }
    size_t m() const {return res_.size();}
    template<typename IT=uint64_t>
    std::vector<IT> to_sigs() const {
        std::vector<IT> ret(m());
        std::transform(res_.data(), res_.data() + m(), ret.begin(), reg2sig<FT>);
        return ret;
    }
};
template<typename FT=float>
struct pmh1a_t {
    using wd = wd_t<FT>;
    using IT = typename wd::IntType;
    struct BufEl {
        uint64_t id;
        double hv, wi;
        uint64_t rv;
    };

    mvt_t<FT> hvals_;
    std::vector<BufEl> buffer_;
    schism::Schismatic<IT> div_;
    std::vector<IT> res_;
    pmh1a_t(size_t m): hvals_(m), div_(m), res_(m) {}

    void update(const IT id, const FT w) {
        if(w <= 0.) return;
        const FT wi = 1. / w;
        uint64_t hi = id;
        uint64_t xi = wy::wyhash64_stateless(&hi);
        auto hv = -std::log((xi >> 12) * 0x1p-52) * wi;
        if(hv >= hvals_.max()) return;
        auto idx = div_.mod(xi);
        if(hvals_.update(idx, hv)) {
            res_[idx] = id;
            if(hv >= hvals_.max()) return;
        }
        buffer_.push_back(BufEl{id, hv, wi, hi});
    }
    void finalize() {
      while(!buffer_.empty()) {
            auto wit = buffer_.begin();
            for(auto oit = buffer_.begin(); oit != buffer_.end(); ++oit) {
                const auto id = oit->id;
                auto &h = oit->hv;
                auto &rng = oit->rv;
                auto wi = oit->wi;
                if(h >= hvals_.max()) continue;
                h += wi * -std::log((wy::wyhash64_stateless(&rng) >> 12) * 0x1p-52);
                if(h >= hvals_.max()) continue;
                auto k = div_.mod(rng);
                if(hvals_.update(k, h)) {
                    res_[k] = id;
                    if(h >= hvals_.max()) continue;
                }
                *wit++ = {id, h, wi, rng};
            }
            buffer_.erase(wit, buffer_.end());
        }
    }
    size_t m() const {return res_.size();}
    template<typename IT=uint64_t>
    std::vector<IT> to_sigs(bool randomize=false) const {
        std::vector<IT> ret(m());
        std::copy(res_.begin(), res_.end(), ret.begin());
        if(randomize)
            std::transform(ret.begin(), ret.end(), ret.begin(),
                           [](uint64_t x) {return wy::wyhash64_stateless(&x);});
        return ret;
    }
};

#endif
