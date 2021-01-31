#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "../fbmh.h"
#include <unordered_map>

namespace py = pybind11;

static std::string standardize_dtype(std::string x) {
    static const std::unordered_map<std::string, std::string> map {
    {"<f8", "d"},
    {"f8", "d"},
    {"<f4", "f"},
    {"f4", "f"},
    {"<u4", "I"},
    {"u4", "I"},
    {"<i4", "I"},
    {"i4", "I"},
    {"<u8", "L"},
    {"u8", "L"},
    {"<i8", "L"},
    {"i8", "L"},
    {"<u2", "H"},
    {"u2", "H"},
    {"<i2", "H"},
    {"i2", "H"},
    {"<u1", "B"},
    {"u1", "B"},
    {"<i1", "B"},
    {"i1", "B"}
    };  
    if(auto it = map.find(x); it != map.end()) x = it->second;
    return x;
}

Sketchers obj2sk(py::object obj) {
    if(py::isinstance<py::int_>(obj))
        return (Sketchers)obj.cast<int>();
    static const std::unordered_map<std::string, Sketchers> map {
        {"bmh", Sketchers::S_BMH1},
        {"bagminhash", Sketchers::S_BMH1},
        {"bagminhash1", Sketchers::S_BMH1},
        {"bagminhash2", Sketchers::S_BMH2},
        {"bmh1", Sketchers::S_BMH1},
        {"bmh2", Sketchers::S_BMH2},
        {"pmh2", Sketchers::S_PMH2},
        {"pmh1", Sketchers::S_PMH1},
        {"pmh1a", Sketchers::S_PMH1A},
        {"pmh", Sketchers::S_PMH2}
    };
    if(py::isinstance<py::str>(obj)) {
        auto os = obj.cast<std::string>();
        std::transform(os.begin(), os.end(), os.begin(), [](auto x) {return std::tolower(x);});
        auto it = map.find(os);
        if(it != map.end()) return it->second;
    }
    throw std::runtime_error("Could not decode to Sketches type");
}

template<typename FT>
py::object perform_hash(py::array_t<FT, py::array::c_style | py::array::forcecast> weights, py::array ids, size_t m, py::object sto) {
    auto stype = obj2sk(sto);
    py::object oret = py::none();
    std::vector<uint64_t> ret;
    auto winf = weights.request(), iinf = ids.request();
    
    switch(standardize_dtype(iinf.format)[0]) {
        case 'f': {py::array_t<float, py::array::c_style | py::array::forcecast> iids(ids);
                   auto iidi = iids.request();
                   ret = wmh::minhash((FT *)winf.ptr, (float *)iidi.ptr, iidi.size, m, stype);}
        break;
        case 'I': case 'i': {py::array_t<uint32_t, py::array::c_style | py::array::forcecast> iids(ids);
                   auto iidi = iids.request();
                   ret = wmh::minhash((uint32_t *)winf.ptr, (uint32_t *)iidi.ptr, iidi.size, m, stype);}
        break;
        case 'L': case 'l':
        default: {
            py::array_t<uint64_t, py::array::c_style | py::array::forcecast> iids(ids);
            auto iidi = iids.request();
            ret = wmh::minhash((uint64_t *)winf.ptr, (uint64_t *)iidi.ptr, iidi.size, m, stype);
        }
    }
    py::array_t<uint64_t> npret(m);
    auto npi = npret.request();
    std::copy(ret.begin(), ret.end(), (uint64_t *)npi.ptr);
    return npret;
}
namespace py = pybind11;
PYBIND11_MODULE(wmh, m) {
    m.def("hash", perform_hash<float>, py::arg("weights"), py::arg("ids"), py::arg("m") = 250, py::arg("stype") = 0,
        "Performs minhash on two arrays: one of weights, and one of ids, sketching using a weighted minhash sketcher. Options: {bmh1, bmh2, pmh1, pmh1a}");
    m.def("hash", perform_hash<double>, py::arg("weights"), py::arg("ids"), py::arg("m") = 250, py::arg("stype") = 0,
        "Performs minhash on two arrays: one of weights, and one of ids, sketching using a weighted minhash sketcher. Options: {bmh1, bmh2, pmh1, pmh1a, pmh2}");
}
