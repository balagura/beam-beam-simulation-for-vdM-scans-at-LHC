// Output data classes for B*B simulation of beam-beam effects in van der Meer
// scans at LHC.
//
// Author V. Balagura, balagura@cern.ch (Nov 2019)
//

#ifndef output_hh 
#define output_hh 1

#include <vector>
#include <string>
#include <complex>
#include <atomic>
#include <memory> // for unique_ptr
#include "gzstream.hh"
#include "bb.hh"

using namespace std;

// -------------------- Output Data classes --------------------
//
// Integrals_Per_Turn accumulate overlap integrals per turn
//
struct Integrals_Per_Turn {
  void resize(int n_turns);
  bool empty();
  void add(int i_turn, double add);
  const vector<double>& integ();
  void write(ostream& os, const string& prefix);
  static string config_name() { return "integrals.per.turn"; }
protected:
  vector<double> double_integ;
  bool double_converted;
  // Unfortunately, at the moment g++ 8.3 does not support atomic<double> +=
  // operations. Therefore, the double accumulators are expressed via
  // atomic<long long int> (by converting double's to long long int's and back). This
  // is much faster than using mutex locks in the threads.
  //
  // In addition, one can not use vector of atomic<long long int> because such
  // atomics do not have a copy constructor. So, dynamic array created by new()
  // is used instead:
  unique_ptr<atomic<long long int>[] > atomic_integ;
  // To not loose precision in conversions double<->long long int, the added
  // doubles are multiplied by a large double constant before converting to
  // long long int. Then, when all threads finish, they are divided by the same
  // constant.  The constant should be large enough to keep precision but
  // small enough so that the accumulator does not exceed the long long int
  // range. The accumulator is the sum of the terms kicker_density(x,y) *
  // point_weight[i]. Since the weights point_weight[i] are normalized, they
  // are all less than 1. So, sum(kicker_density(x,y) * point_weight[i]) <
  // sum(kicker_density(x,y)) = 1. Therefore, the maximal long long int value,
  // LLONG_MAX, is chosen as a constant.
};

//
// Avr_XY_Per_Turn - accumulate x+iy complex coordinate of the kicked bunch center-of-mass
struct Avr_XY_Per_Turn {
  void resize(int n_turns);
  bool empty();
  //
  void init(double max_xy_radius, int n_points); // sets atomic_avr_xy_converter
  void add(int i_turn, complex<double> add);
  void write(ostream& os, const string& prefix);
  static string config_name() { return "avr.xy.per.turn"; }
  Avr_XY_Per_Turn() : size(0), atomic_avr_xy_converter(0) {}
protected:
  unique_ptr<atomic<long long int>[]> atomic_avr_x, atomic_avr_y;
  int size;
  double atomic_avr_xy_converter;
};

//
// Integrals_per_particle and Avr_XY_Per_Particle - overlap integrals and average
// x+iy coordinates accumulated by one macroparticle during groups of many
// turns (NO_BB, ADIABATIC, STABILIZATION, BB). Ie. the kicked bunch density
// is represented in this case by a delta-function concentrated at only one
// particle which makes, however, many betatron oscillations (many turns).
//
template <class T>
// T = double for Integrals_Per_Particle, complex<double> for Avr_XY_Per_Particle
struct Per_Particle_Base : protected array<vector<T>, PHASES> { // [phase][particle]
  void resize(int n_points);
  bool empty();
  void add(int phase, int particle, T add);
  void normalize_by_n_turns(int n_turns, int phase, int particle); // parallelized, so per particle
  array<T, PHASES> weighted_average(const vector<double>& particle_weights);
  void write(ostream& os, const string& prefix);
  static string config_name();
};
struct Integrals_Per_Particle : public Per_Particle_Base<double> {
  static string config_name() { return "integrals.per.particle"; }
};
struct Avr_XY_Per_Particle : public Per_Particle_Base<complex<double> > {
  static string config_name() { return "avr.xy.per.particle"; }
};

// Points store raw 4D coordinates for selected turns in
// Points[turn][particle][0/1 for x/y]
//
// This vector will be filled by all threads in parallel and when all finish,
// will be written to file. Unfortunately, direct writing to the same file by
// many threads would block them, this dictates storing large v_points in
// memory.
struct Points : protected vector<vector<array<complex<double>, 2> > > { // [turn][particle][0/1 for x/y]
  void resize(int n_points,
	      int n_turns,
	      int select_one_turn_out_of);
  using vector<vector<array<complex<double>, 2> > >::empty;
  void add(int turn, int particle, complex<double> xZ, complex<double> yZ);
  void write(ostream& os, const string& prefix);
  static string config_name() { return "points"; }
protected:
  int n_selected_turns;
};

//
// Create a vector of the classes above (denoted generically as class Data)
// and a boolean flag write_it controlling whether the corresponding data
// should be written out.
//
template<class Data>
struct Output {
  Output() : write_it(false) {}
  vector<Data> data; // [ip]
  bool write_it;
  ogzstream file;
  static string config_name() { return Data::config_name(); }
};
//
// Main class: std::tuple<...> of everything above.
//
// It is made as a tuple to implement some rudimentary looping over types for
// initialization and writing to files. Unfortunately, the support of
// metaprogramming in C++ is very limited. Here, the looping is performed by
// the cryptic line
//
// apply([&f](auto& ... x) { (..., f(x)); }, *((Outputs_tuple*)this));
//
// using std::apply(function, tuple) and a fold expression both appearing in
// C++17. The "function" is a lambda [&f](auto& ... x) which makes "folding"
// (... operation pack) where "operation" is just the comma-operator "," and
// the pack contains f() applied in this way to all members of the tuple.
// The function f() should be applicable to all elements of the tuple.
//
// In principle, it should be easy for the compiler to call f() sequentially
// for all tuple's elements, as it knows all the types, positions in memory
// etc., but expressing it in C++ is cumbersome.
//
// In spite of that, here, it is made generic with apply+fold in order to be
// flexible in activating / deactivating the corresponding algorithms and the
// following outputs using configuration switches, and, possibly, for adding
// new output types (or removing existing).
//
// In addition to the generic initialization + output, for the specific
// algorithms one needs to access the structures in the tuple individually
// (knowing its type). They can be done using a helper member function like
//
// get<Integrals_Per_Particle>(ip)
//
// It returns the data of Integrals_Per_Particle type for the interaction point
// "ip" (counting from zero). Integrals_Per_Particle should be a unique class
// name in the tuple (otherwise, if the tuple contains two structures of the
// same type, the compiler will produce an error, see the documentation on
// std::get(std::tuple) which is used internally).
//
// To check whether the particular class and the corresponding algorithm has
// been chosen in the configuration, one can use
//
// is_active<Integrals_Per_Particles>(ip)
//
// "ip" is needed here since the algorithms are configurable per interaction
// point.
//
typedef tuple<Output<Integrals_Per_Turn>,
	      Output<Avr_XY_Per_Turn>,
	      Output<Integrals_Per_Particle>,
	      Output<Avr_XY_Per_Particle>,
	      Output<Points> > Outputs_tuple;
struct Outputs : public Outputs_tuple {
  Outputs(int n_ip,
	  int n_points,
	  double max_xy_radius,
	  const string& output_options, // white-space separated list of options
	  const string& output_dir); // for double <-> llong int converter in Avr_XY_Per_Turn
  void reset(int n_ip,
	     int n_points,
	     int n_turns,
	     int select_one_turn_out_of);
  void write(int step);
  // access to tuple individual components
  template<class T> T&   get      (int ip) { return std::get<Output<T> >(*this).data[ip]; }
  // should we process the data of type "T"?
  template<class T> bool is_active(int ip) { return !std::get<Output<T> >(*this).data[ip].empty(); }
};

// ------------------------------------------------------------
//
// Per_Particle_Base<T> = array<vector<T>, PHASES>  [phase][particle]
//
template <class T>
void Per_Particle_Base<T>::resize(int n_points) {
  for (int phase=0; phase<PHASES; ++phase) (*this)[phase].assign(n_points, 0.); // 0+i0 for complex
}
template <class T>
bool Per_Particle_Base<T>::empty() { return (*this)[0].empty(); }
template <class T>
// initialized by zeros as any vector of doubles
void Per_Particle_Base<T>::add(int phase, int particle, T add) {
  (*this)[phase][particle] += add;
}
template <class T>
void Per_Particle_Base<T>::normalize_by_n_turns(int n_turns, int phase, int particle) {
  // parallelized, so per particle
  (*this)[phase][particle] /= n_turns;
}
template <class T>
array<T, PHASES> Per_Particle_Base<T>::weighted_average(const vector<double>& particle_weights) {
  array<T, PHASES> avr;
  fill(avr.begin(), avr.end(), 0.);
  for (int phase=0; phase<PHASES; ++phase) {
    const vector<T>& v = (*this)[phase];
    T& a = avr[phase];
    for (size_t i=0; i<v.size(); ++i) a += particle_weights[i] * v[i];
  }
  return avr;
}

// operator<<(ostream&, complex<double>) is already overloaded in C++ standard
// and prints (real,imag) instead of "real imag" needed here. Use "dump"
// below instead:
inline void dump(ostream& os, double x)          { os << x; }
inline void dump(ostream& os, complex<double> x) { os << real(x) << " " << imag(x); }
inline bool isnan(complex<double> x) { return isnan(real(x)) || isnan(imag(x)); }
//
template <class T>
void Per_Particle_Base<T>::write(ostream& os, const string& prefix) {
  for (int phase=0; phase<PHASES; ++phase) {
    for (size_t i = 0; i < (*this)[phase].size(); ++i) {
      T x = (*this)[phase][i];
      if (!isnan(x)) { // if in some phase n_turn=0, normalization to such n_turns gives nan
	os << prefix << " " << phase << " "  << i << " ";
	dump(os, x);
	os << '\n';
      }
    }
  }
}

#endif
