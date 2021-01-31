#ifndef BAGMINHASH_H__
#define BAGMINHASH_H__
#include <stdexcept>
#include <cassert>
#include "aesctr/wy.h"
#include <queue>

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


    void update(size_t index, FT x) {
        const auto sz = data_.size();
        const auto mv = getm();
        while(x < data_[index]) {
            data_[index] = x;
            index = mv + (index >> 1);
            if(index >= sz) break;
            size_t lhi = (index - mv) << 1;
            size_t rhi = lhi + 1;
            x = std::max(data_[lhi], data_[rhi]);
        }
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
    static constexpr FT maxv = std::numeric_limits<FT>::max();
    static const IT maxi = ft2it(maxv);
    using IntType = IT;
};

template<typename FT=double, typename IT=DefIT<FT>>
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
        return wd::ft2it(maxq_);
    }
    IT widxmin() const {
        return wd::ft2it(minp_);
    }
    bool partially_relevant() const {
        return wd::it2ft(widxmin() + 1) <= weight_;
    }
    bool fully_relevant() const {
        return wd::it2ft(widxmax()) <= weight_;
    }
    bool can_split() const {
        //std::fprintf(stderr, "indices: %zu %zu\n", size_t(rhi), size_t(lhi));
        return widxmax() > widxmin() + 1;
    }
    // Note: > is reversed, for use in pq
    bool operator>(const poisson_process_t &o) const {return x_ < o.x_;}
    bool operator<(const poisson_process_t &o) const {return x_ > o.x_;}
    void step(size_t m) {
        // Top 52-bits as U01 for exponential with weight of q - p,
        // bottom logm bits for index
        uint64_t xi = wy::wyhash64_stateless(&wyv_);
        x_ += -std::log(static_cast<double>(xi >> 12) * 0x1p-52) / (maxq_ - minp_);
        idx_ = xi % m;
    }
    poisson_process_t split() {
        auto wmin = widxmin(), wmax = widxmax();
        uint64_t midpoint = (uint64_t(wmin) + wmax) / 2;
        //std::fprintf(stderr, "%zu-%zu are splitting at %zu\n", size_t(wmin), size_t(wmax), size_t(midpoint));
        double midval = wd::it2ft(midpoint);
        uint64_t xval = wd::ft2it(x_) ^ wy::wyhash64_stateless(&midpoint);
        const double p = (midval - minp_) / (maxq_ - minp_);
        const double rv = (wy::wyhash64_stateless(&xval) >> 12) * 0x1p-52;
        //std::fprintf(stderr, "splitting %g->%g at %g with prob %g\n", minp_, maxq_, midval, p);
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



template<typename FT=float>
struct bmh_t {
    using wd = wd_t<FT>;
    using IT = typename wd::IntType;
    using PoissonP = poisson_process_t<FT, IT>;

    size_t maxspace_;
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
    auto m() const {return hvals_.getm();}

    bmh_t(size_t m): hvals_(m) {
        heap_.getc().reserve(m);
    }
    void update_1(IT id, FT w) {
        if(w <= 0.) return;
        PoissonP p(id, w);
        const auto _m = m();
        p.step(_m);
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
                    pp.step(_m);
                    if(pp.fully_relevant()) hvals_.update(pp.idx_, pp.x_);
                    if(pp.partially_relevant()) heap_.push(std::move(pp));
                }
                //std::fprintf(stderr, "Finishing subloop at %zu/%zu\n", mainiternum, subin);
            }
            if(p.fully_relevant()) {
                p.step(_m);
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
    static inline uint64_t reg2sig(FT v) {
        uint64_t t = 0;
        std::memcpy(&t, &v, std::min(sizeof(v), sizeof(uint64_t)));
        return wy::wyhash64_stateless(&t);
    }
    template<typename IT=uint64_t>
    std::vector<IT> to_sigs() const {
        std::vector<IT> ret(m());
        for(size_t i = 0; i < m(); ++i) std::fprintf(stderr, "sig %zu is %g\n", i, hvals_[i]);
        //std::transform(&hvals_.data_[0], &hvals_.data_[mv], ret.begin(), [](auto x) {uint64_t xv = 0; std::memcpy(&x, &xv, sizeof(x)); return wy::wyhash64_stateless(&xv);});
        std::transform(hvals_.data(), hvals_.data() + m(), ret.begin(), reg2sig);
        return ret;
    }
};

#endif
