// BindOptiX compiled core: particle transport of stage Z3 (docs/roadmap.md, stage Z11).
//
// One function, simulate_particle, advances one cell through all steps; it is compiled for the
// CPU (OpenMP over cells) and for CUDA (one thread per cell). It reproduces the Python reference
// bindoptix.particles.run: same equations, same order of operations, same counter-based random
// numbers (SplitMix64 over seed, particle id, step and lane), double precision throughout.
// Results agree with the reference up to the rounding of exp/log/sin/cos (a few ulp).
#pragma once

#include <cmath>
#include <cstdint>

#if defined(__CUDACC__)
#define BX_HD __host__ __device__
#else
#define BX_HD
#endif

namespace bindoptix {

constexpr int FREE = 0, FLOOR = 1, CAPTURED = 2, STUCK = 3, EXITED = 4, WAITING = 5;
constexpr double PI = 3.14159265358979323846;

// ---- counter-based random numbers (bindoptix.rng) -------------------------------------------

BX_HD inline std::uint64_t mix(std::uint64_t x) {
  std::uint64_t z = x + 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

BX_HD inline std::uint64_t hash4(std::uint64_t seed, std::uint64_t pid, std::uint64_t step,
                                 std::uint64_t lane) {
  std::uint64_t h = mix(seed);
  h = mix(h ^ pid);
  h = mix(h ^ step);
  return mix(h ^ lane);
}

BX_HD inline double uniform(std::uint64_t seed, std::uint64_t pid, std::uint64_t step,
                            std::uint64_t lane) {
  return static_cast<double>(hash4(seed, pid, step, lane) >> 11) * (1.0 / 9007199254740992.0);
}

BX_HD inline double normal(std::uint64_t seed, std::uint64_t pid, std::uint64_t step,
                           std::uint64_t lane) {
  const double u1 = 1.0 - uniform(seed, pid, step, 2 * lane);
  const double u2 = uniform(seed, pid, step, 2 * lane + 1);
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * PI * u2);
}

// ---- parameters ---------------------------------------------------------------------------

struct Params {
  // channel flow (bindoptix.flow.DuctFlow)
  double width, height, pressure_gradient;
  int n_terms;
  int has_flow;
  // cell (bindoptix.cells.CellModel)
  double diffusion, radius, vsed[3], swim_speed, p_tumble, tumble_cos;
  // channel model (bindoptix.particles.ChannelModel)
  int floor_contact;  // 1 contact, 0 reflect
  double k_ns;        // may be +inf
  double slip, floor_diffusion_factor, x_out, contact_gap;
  int brownian, wall_hydrodynamics;
  // spots: convex counter-clockwise polygons
  int n_spots;
  const int* poly_offset;   // n_spots + 1 vertex offsets
  const double* poly_xy;    // vertices (x, y)
  const double* rates;      // n_steps * n_spots capture rates per step (may be +inf)
  // run
  std::uint64_t seed;
  double dt, t0;
  int n_steps;
  std::int64_t step0;
  int record_visits;
  int n_snap;
  const int* snap_step;  // 1-based step within this run at which each snapshot is taken
};

struct State {
  double* x;           // n * 3
  std::int64_t* state; // n
  std::int64_t* spot;  // n
  double* t_event;     // n
  double* e;           // n * 3
  const double* t_start;
  const std::uint64_t* pid;
  double* visit_entry; // n * n_spots
  double* visit_time;
  double* visit_off;
  double* off_floor;   // n
  double* snaps;       // n_snap * n * 3
  std::int64_t n;
};

// ---- channel flow ---------------------------------------------------------------------------

// The Fourier series run over odd n with k_n = nπ/h. Instead of two exponentials and a sine per
// term, powers and the sine are advanced by recurrences (A^(n+2) = A^n A², sin((n+2)θ) =
// 2 cos 2θ sin nθ − sin((n−2)θ)); this is ~20 times faster and agrees with the termwise
// evaluation of the reference to about 1e-13 relative.

BX_HD inline double velocity(const Params& p, double y, double z) {
  const double h = p.height, hw = 0.5 * p.width;
  const bool inside = std::fabs(y) <= hw && z >= 0.0 && z <= h;
  if (!inside) return 0.0;
  const double ay = std::fabs(y < -hw ? -hw : (y > hw ? hw : y));
  const double a1 = std::exp(PI * (ay - hw) / h), a2 = a1 * a1;
  const double b1 = std::exp(-2.0 * PI * ay / h), b2 = b1 * b1;
  const double c1 = std::exp(-2.0 * PI * hw / h), c2 = c1 * c1;
  const double th = PI * z / h, c2t = 2.0 * std::cos(2.0 * th);
  double an = a1, bn = b1, cn = c1, s_prev = -std::sin(th), s_n = std::sin(th), sum = 0.0;
  for (int i = 0; i < p.n_terms; ++i) {
    const double n = 2.0 * i + 1.0;
    const double ratio = an * (1.0 + bn) / (1.0 + cn);
    sum += (1.0 / (n * n * n)) * (1.0 - ratio) * s_n;
    const double s_next = c2t * s_n - s_prev;
    s_prev = s_n;
    s_n = s_next;
    an *= a2;
    bn *= b2;
    cn *= c2;
  }
  return 4.0 * h * h * p.pressure_gradient / (PI * PI * PI) * sum;
}

BX_HD inline double floor_shear(const Params& p, double y) {
  const double h = p.height, hw = 0.5 * p.width;
  if (std::fabs(y) > hw) return 0.0;
  const double ay = std::fabs(y);
  const double a1 = std::exp(PI * (ay - hw) / h), a2 = a1 * a1;
  const double b1 = std::exp(-2.0 * PI * ay / h), b2 = b1 * b1;
  const double c1 = std::exp(-2.0 * PI * hw / h), c2 = c1 * c1;
  double an = a1, bn = b1, cn = c1, sum = 0.0;
  for (int i = 0; i < p.n_terms; ++i) {
    const double n = 2.0 * i + 1.0;
    sum += (1.0 / (n * n)) * (an * (1.0 + bn) / (1.0 + cn));
    an *= a2;
    bn *= b2;
    cn *= c2;
  }
  return 0.5 * p.pressure_gradient * h - 4.0 * h * p.pressure_gradient / (PI * PI) * sum;
}

// ---- wall hydrodynamics (bindoptix.wall) ------------------------------------------------------

BX_HD inline double lambda_perp(double gap, double a) {
  const double h = gap > 0.0 ? gap : 0.0;
  return (6.0 * h * h + 2.0 * a * h) / (6.0 * h * h + 9.0 * a * h + 2.0 * a * a);
}

BX_HD inline double dlambda_perp(double gap, double a) {
  const double h = gap > 0.0 ? gap : 0.0;
  const double num = 6.0 * h * h + 2.0 * a * h;
  const double den = 6.0 * h * h + 9.0 * a * h + 2.0 * a * a;
  return ((12.0 * h + 2.0 * a) * den - num * (12.0 * h + 9.0 * a)) / (den * den);
}

BX_HD inline double lambda_par(double gap, double a) {
  const double h = gap > 0.0 ? gap : 0.0;
  const double r = a / (a + h);
  const double lam = 1.0 - 9.0 / 16.0 * r + std::pow(r, 3) / 8.0 - 45.0 / 256.0 * std::pow(r, 4) -
                     std::pow(r, 5) / 16.0;
  return lam < 0.0 ? 0.0 : (lam > 1.0 ? 1.0 : lam);
}

// ---- spots ----------------------------------------------------------------------------------

BX_HD inline bool in_spot(const Params& p, int k, double x, double y) {
  const int v0 = p.poly_offset[k], v1 = p.poly_offset[k + 1], nv = v1 - v0;
  for (int i = 0; i < nv; ++i) {
    const double ax = p.poly_xy[2 * (v0 + i)], ay = p.poly_xy[2 * (v0 + i) + 1];
    const int j = v0 + (i + 1) % nv;
    const double bx = p.poly_xy[2 * j], by = p.poly_xy[2 * j + 1];
    if ((bx - ax) * (y - ay) - (by - ay) * (x - ax) < 0.0) return false;
  }
  return true;
}

// ---- tumbling (bindoptix.particles._turn) -----------------------------------------------------

BX_HD inline void turn(double* e, double cos_a, double phi) {
  const double hx = std::fabs(e[2]) < 0.9 ? 0.0 : 1.0;
  const double hz = std::fabs(e[2]) < 0.9 ? 1.0 : 0.0;
  // p1 = e × helper, helper = (hx, 0, hz)
  double p1[3] = {e[1] * hz - e[2] * 0.0, e[2] * hx - e[0] * hz, e[0] * 0.0 - e[1] * hx};
  const double n1 = std::sqrt(p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]);
  for (double& v : p1) v /= n1;
  const double p2[3] = {e[1] * p1[2] - e[2] * p1[1], e[2] * p1[0] - e[0] * p1[2],
                        e[0] * p1[1] - e[1] * p1[0]};
  const double sa = std::sqrt(1.0 - cos_a * cos_a > 0.0 ? 1.0 - cos_a * cos_a : 0.0);
  const double c = std::cos(phi), s = std::sin(phi);
  double o[3];
  for (int i = 0; i < 3; ++i) o[i] = cos_a * e[i] + sa * (c * p1[i] + s * p2[i]);
  const double no = std::sqrt(o[0] * o[0] + o[1] * o[1] + o[2] * o[2]);
  for (int i = 0; i < 3; ++i) e[i] = o[i] / no;
}

// ---- one cell through all steps --------------------------------------------------------------

BX_HD inline void simulate_particle(const Params& p, const State& s, std::int64_t i) {
  const double d = p.brownian ? p.diffusion : 0.0;
  const double sig = std::sqrt(2.0 * d * p.dt);
  const double sig_floor = std::sqrt(2.0 * d * p.floor_diffusion_factor * p.dt);
  const double w2 = 0.5 * p.width, h = p.height;
  double* x = s.x + 3 * i;
  double* e = s.e + 3 * i;
  std::int64_t st = s.state[i];
  const std::uint64_t pid = s.pid[i];
  const std::uint64_t seed = p.seed;
  int snap_k = 0;
  for (int k = 1; k <= p.n_steps; ++k) {
    const std::uint64_t step = static_cast<std::uint64_t>(p.step0 + k);
    const double t = p.t0 + k * p.dt;
    if (st == WAITING && s.t_start[i] <= t - p.dt + 1e-9 * p.dt) {
      st = FREE;
      if (p.floor_contact && x[2] <= 0.0) {
        st = FLOOR;
        x[2] = 0.0;
      }
    }
    if (st == FREE) {
      const double u = p.has_flow ? velocity(p, x[1], x[2]) : 0.0;
      double drift[3] = {u, 0.0, 0.0};
      const double z0 = x[2];
      double scale[3] = {1.0, 1.0, 1.0};
      double d_perp = d;
      if (p.wall_hydrodynamics) {
        const double lp = lambda_perp(z0, p.radius), ls = lambda_par(z0, p.radius);
        drift[0] += p.vsed[0] * ls;
        drift[1] += p.vsed[1] * ls;
        drift[2] += p.vsed[2] * lp + d * dlambda_perp(z0, p.radius);
        for (int c = 0; c < 3; ++c) drift[c] += p.swim_speed * e[c];
        scale[0] = scale[1] = std::sqrt(ls);
        scale[2] = std::sqrt(lp);
        d_perp = d * lp;
      } else {
        for (int c = 0; c < 3; ++c) drift[c] += p.vsed[c] + p.swim_speed * e[c];
      }
      double q[3];
      for (int c = 0; c < 3; ++c)
        q[c] = x[c] + drift[c] * p.dt + sig * scale[c] * normal(seed, pid, step, c);
      if (q[1] > w2) q[1] = 2.0 * w2 - q[1];
      if (q[1] < -w2) q[1] = -2.0 * w2 - q[1];
      if (q[2] > h) q[2] = 2.0 * h - q[2];
      if (!p.floor_contact) {
        q[2] = std::fabs(q[2]);
      } else {
        const double zc = p.contact_gap;
        bool hit = q[2] <= zc;
        if (d > 0.0 && !hit) {
          const double g0 = z0 - zc > 0.0 ? z0 - zc : 0.0;
          const double g1 = q[2] - zc > 0.0 ? q[2] - zc : 0.0;
          const double pb = d_perp > 0.0 ? std::exp(-g0 * g1 / (d_perp * p.dt)) : 0.0;
          hit = uniform(seed, pid, step, 10) < pb;
        }
        if (hit) {
          q[2] = zc;
          st = FLOOR;
        }
      }
      x[0] = q[0];
      x[1] = q[1];
      x[2] = q[2];
    }
    if (st == FLOOR) {
      const double slide = p.slip * p.radius * floor_shear(p, x[1]);
      x[0] += slide * p.dt + sig_floor * normal(seed, pid, step, 20);
      x[1] += sig_floor * normal(seed, pid, step, 21);
      x[1] = x[1] < -w2 ? -w2 : (x[1] > w2 ? w2 : x[1]);
      double r = p.k_ns;
      int which = -1;
      for (int sp = 0; sp < p.n_spots; ++sp) {
        if (in_spot(p, sp, x[0], x[1])) {
          which = sp;
          r = p.rates[static_cast<std::int64_t>(k - 1) * p.n_spots + sp];
          break;
        }
      }
      if (p.record_visits) {
        if (which < 0) {
          s.off_floor[i] += p.dt;
        } else {
          const std::int64_t at = i * p.n_spots + which;
          if (std::isnan(s.visit_entry[at])) {
            s.visit_entry[at] = t;
            s.visit_off[at] = s.off_floor[i];
          }
          s.visit_time[at] += p.dt;
        }
      }
      const double prob = std::isinf(r) ? 1.0 : -std::expm1(-r * p.dt);
      const bool caught = uniform(seed, pid, step, 22) < prob;
      if (caught) {
        st = which >= 0 ? CAPTURED : STUCK;
        if (which >= 0) s.spot[i] = which;
        s.t_event[i] = t;
      } else if (p.swim_speed > 0.0 && e[2] > 0.0) {
        st = FREE;
      }
    }
    if (p.swim_speed > 0.0 && (st == FREE || st == FLOOR)) {
      if (uniform(seed, pid, step, 30) < p.p_tumble) {
        turn(e, p.tumble_cos, 2.0 * PI * uniform(seed, pid, step, 31));
      }
    }
    if ((st == FREE || st == FLOOR) && x[0] > p.x_out) {
      st = EXITED;
      s.t_event[i] = t;
    }
    while (snap_k < p.n_snap && p.snap_step[snap_k] == k) {
      double* o = s.snaps + (static_cast<std::int64_t>(snap_k) * s.n + i) * 3;
      o[0] = x[0];
      o[1] = x[1];
      o[2] = x[2];
      ++snap_k;
    }
  }
  s.state[i] = st;
}

// backends
void run_cpu(const Params& p, const State& s);
bool cuda_available();
void run_cuda(const Params& p, const State& s, std::int64_t n_vertices);

}  // namespace bindoptix
