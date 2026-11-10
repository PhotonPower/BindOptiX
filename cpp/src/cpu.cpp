// CPU backend: OpenMP over cells (cells are independent).
#include "bindoptix/core.hpp"

namespace bindoptix {

void run_cpu(const Params& p, const State& s) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (std::int64_t i = 0; i < s.n; ++i) {
    simulate_particle(p, s, i);
  }
}

#if !defined(BINDOPTIX_HAVE_CUDA)
bool cuda_available() { return false; }
void run_cuda(const Params&, const State&, std::int64_t) {}
#endif

}  // namespace bindoptix
