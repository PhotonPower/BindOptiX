// pybind11 module bindoptix._core (stage Z11). The Python facade is bindoptix.fast.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cstdint>
#include <stdexcept>
#include <string>

#include "bindoptix/core.hpp"

namespace py = pybind11;
using namespace bindoptix;

template <typename T>
using Arr = py::array_t<T, py::array::c_style | py::array::forcecast>;

namespace {

template <typename T>
T* mut(Arr<T>& a, std::int64_t expected, const char* name) {
  if (a.size() != expected) {
    throw std::invalid_argument(std::string(name) + ": wrong size");
  }
  return a.mutable_data();
}

void run(py::dict prm, Arr<double> x, Arr<std::int64_t> state, Arr<std::int64_t> spot,
         Arr<double> t_event, Arr<double> e, Arr<double> t_start, Arr<std::uint64_t> pid,
         Arr<double> visit_entry, Arr<double> visit_time, Arr<double> visit_off,
         Arr<double> off_floor, Arr<double> snaps, Arr<int> poly_offset, Arr<double> poly_xy,
         Arr<double> rates, Arr<int> snap_step, const std::string& backend) {
  Params p{};
  auto g = [&](const char* k) { return prm[k].cast<double>(); };
  auto gi = [&](const char* k) { return prm[k].cast<int>(); };
  p.width = g("width");
  p.height = g("height");
  p.pressure_gradient = g("pressure_gradient");
  p.n_terms = gi("n_terms");
  p.has_flow = gi("has_flow");
  p.diffusion = g("diffusion");
  p.radius = g("radius");
  p.vsed[0] = g("vsed_x");
  p.vsed[1] = g("vsed_y");
  p.vsed[2] = g("vsed_z");
  p.swim_speed = g("swim_speed");
  p.p_tumble = g("p_tumble");
  p.tumble_cos = g("tumble_cos");
  p.floor_contact = gi("floor_contact");
  p.k_ns = g("k_ns");
  p.slip = g("slip");
  p.floor_diffusion_factor = g("floor_diffusion_factor");
  p.x_out = g("x_out");
  p.contact_gap = g("contact_gap");
  p.brownian = gi("brownian");
  p.wall_hydrodynamics = gi("wall_hydrodynamics");
  p.n_spots = gi("n_spots");
  p.seed = prm["seed"].cast<std::uint64_t>();
  p.dt = g("dt");
  p.t0 = g("t0");
  p.n_steps = gi("n_steps");
  p.step0 = prm["step0"].cast<std::int64_t>();
  p.record_visits = gi("record_visits");
  p.n_snap = gi("n_snap");
  const std::int64_t n = state.size();
  const std::int64_t nv = p.record_visits ? n * p.n_spots : 0;
  if (poly_offset.size() != p.n_spots + 1 || rates.size() != std::int64_t(p.n_steps) * p.n_spots ||
      snap_step.size() != p.n_snap) {
    throw std::invalid_argument("spot, rate or snapshot arrays have the wrong size");
  }
  p.poly_offset = poly_offset.data();
  p.poly_xy = poly_xy.data();
  p.rates = rates.data();
  p.snap_step = snap_step.data();
  State s{mut(x, 3 * n, "x"),
          mut(state, n, "state"),
          mut(spot, n, "spot"),
          mut(t_event, n, "t_event"),
          mut(e, 3 * n, "e"),
          t_start.data(),
          pid.data(),
          nv ? mut(visit_entry, nv, "visit_entry") : nullptr,
          nv ? mut(visit_time, nv, "visit_time") : nullptr,
          nv ? mut(visit_off, nv, "visit_off") : nullptr,
          p.record_visits ? mut(off_floor, n, "off_floor") : nullptr,
          mut(snaps, std::int64_t(p.n_snap) * 3 * n, "snaps"),
          n};
  if (t_start.size() != n || pid.size() != n) throw std::invalid_argument("t_start/pid size");
  if (backend == "cuda") {
    if (!cuda_available()) throw std::runtime_error("CUDA backend not available");
    py::gil_scoped_release release;
    run_cuda(p, s, poly_xy.size() / 2);
  } else if (backend == "cpu") {
    py::gil_scoped_release release;
    run_cpu(p, s);
  } else {
    throw std::invalid_argument("backend must be 'cpu' or 'cuda'");
  }
}

Arr<double> uniform_v(std::uint64_t seed, Arr<std::uint64_t> pid, std::uint64_t step,
                      std::uint64_t lane) {
  Arr<double> out(pid.size());
  auto* o = out.mutable_data();
  for (py::ssize_t i = 0; i < pid.size(); ++i) o[i] = uniform(seed, pid.data()[i], step, lane);
  return out;
}

Arr<double> normal_v(std::uint64_t seed, Arr<std::uint64_t> pid, std::uint64_t step,
                     std::uint64_t lane) {
  Arr<double> out(pid.size());
  auto* o = out.mutable_data();
  for (py::ssize_t i = 0; i < pid.size(); ++i) o[i] = normal(seed, pid.data()[i], step, lane);
  return out;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "BindOptiX compiled core (stage Z11): particle transport on CPU and CUDA";
  m.def("run", &run);
  m.def("uniform", &uniform_v);
  m.def("normal", &normal_v);
  m.def("cuda_built", []() {
#if defined(BINDOPTIX_HAVE_CUDA)
    return true;
#else
    return false;
#endif
  });
  m.def("cuda_available", &cuda_available);
}
