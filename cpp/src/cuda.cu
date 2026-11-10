// CUDA backend (stage Z11): one thread per cell, the same simulate_particle as on the CPU.
// Built with -DBINDOPTIX_ENABLE_CUDA=ON on the large machine (RTX 6000 Ada sm_89, RTX 3090 sm_86).
#include <cuda_runtime.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "bindoptix/core.hpp"

namespace bindoptix {

namespace {

void check(cudaError_t err, const char* what) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string("CUDA ") + what + ": " + cudaGetErrorString(err));
  }
}

template <typename T>
struct DeviceBuffer {
  T* ptr = nullptr;
  std::size_t count = 0;
  explicit DeviceBuffer(std::size_t n) : count(n) {
    if (n) check(cudaMalloc(&ptr, n * sizeof(T)), "malloc");
  }
  ~DeviceBuffer() {
    if (ptr) cudaFree(ptr);
  }
  void upload(const T* host) {
    if (count) check(cudaMemcpy(ptr, host, count * sizeof(T), cudaMemcpyHostToDevice), "upload");
  }
  void download(T* host) const {
    if (count) check(cudaMemcpy(host, ptr, count * sizeof(T), cudaMemcpyDeviceToHost), "download");
  }
  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;
};

__global__ void kernel(Params p, State s) {
  const std::int64_t i = static_cast<std::int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < s.n) simulate_particle(p, s, i);
}

}  // namespace

bool cuda_available() {
  int n = 0;
  return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

void run_cuda(const Params& hp, const State& hs, std::int64_t n_vertices) {
  const std::int64_t n = hs.n;
  const std::size_t ns = static_cast<std::size_t>(hp.n_spots);
  const std::size_t nv = hp.record_visits ? static_cast<std::size_t>(n) * ns : 0;
  DeviceBuffer<int> off(ns + 1), snap(hp.n_snap);
  DeviceBuffer<double> xy(2 * static_cast<std::size_t>(n_vertices));
  DeviceBuffer<double> rates(static_cast<std::size_t>(hp.n_steps) * ns);
  DeviceBuffer<double> x(3 * n), te(n), e(3 * n), ts(n), ve(nv), vt(nv), vo(nv);
  DeviceBuffer<double> of(hp.record_visits ? n : 0), sn(static_cast<std::size_t>(hp.n_snap) * 3 * n);
  DeviceBuffer<std::int64_t> st(n), sp(n);
  DeviceBuffer<std::uint64_t> pid(n);
  off.upload(hp.poly_offset);
  snap.upload(hp.snap_step);
  xy.upload(hp.poly_xy);
  rates.upload(hp.rates);
  x.upload(hs.x);
  te.upload(hs.t_event);
  e.upload(hs.e);
  ts.upload(hs.t_start);
  st.upload(hs.state);
  sp.upload(hs.spot);
  pid.upload(hs.pid);
  if (hp.record_visits) {
    ve.upload(hs.visit_entry);
    vt.upload(hs.visit_time);
    vo.upload(hs.visit_off);
    of.upload(hs.off_floor);
  }
  Params dp = hp;
  dp.poly_offset = off.ptr;
  dp.poly_xy = xy.ptr;
  dp.rates = rates.ptr;
  dp.snap_step = snap.ptr;
  State ds{x.ptr, st.ptr, sp.ptr, te.ptr, e.ptr, ts.ptr, pid.ptr,
           ve.ptr, vt.ptr, vo.ptr, of.ptr, sn.ptr, n};
  const int block = 128;
  const int grid = static_cast<int>((n + block - 1) / block);
  if (n) {
    kernel<<<grid, block>>>(dp, ds);
    check(cudaGetLastError(), "launch");
    check(cudaDeviceSynchronize(), "synchronize");
  }
  x.download(hs.x);
  te.download(hs.t_event);
  e.download(hs.e);
  st.download(hs.state);
  sp.download(hs.spot);
  sn.download(hs.snaps);
  if (hp.record_visits) {
    ve.download(hs.visit_entry);
    vt.download(hs.visit_time);
    vo.download(hs.visit_off);
    of.download(hs.off_floor);
  }
}

}  // namespace bindoptix
