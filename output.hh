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
#include <gzstream.h>
#include <config.hh>
#include <n.hh>

using namespace std;

// -------------------- Output Data classes --------------------
//
// Integrals accumulate overlap integrals per turn
//
struct Integrals {
  void resize(const N& n);
  void resize(int N_turns);
  bool empty();
  void add(int i_turn, double add);
  const vector<double>& integ();
  void write(ostream& os, const string& prefix);
  static string config_name() { return "integrals"; }
  Integrals() : size(0) {}
protected:
  vector<double> double_integ;
  // Unfortunately, at the moment g++ 8.3 does not support atomic<double> +=
  // operations. Therefore, the double accumulators are expressed via
  // atomic<long long int> (by converting double's to long long int's and back). This
  // is much faster than using mutex locks in the threads.
  //
  // In addition, one can not use vector of atomic<long long int> because such
  // atomics do not have a copy constructor. So, dynamic array created by new()
  // is used instead:
  unique_ptr<atomic<long long int>[] > atomic_integ;
  int size;
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
// Avr_Z - accumulate x+iy complex coordinate of the kicked bunch center-of-mass
struct Avr_Z {
  void resize(const N& n);
  void resize(int N_turns);
  bool empty();
  //
  void init(double max_xy_radius, int n_points); // sets atomic_avr_xy_converter
  void add(int i_turn, complex<double> add);
  void write(ostream& os, const string& prefix);
  static string config_name() { return "centers"; }
  Avr_Z() : size(0), atomic_avr_xy_converter(0) {}
protected:
  unique_ptr<atomic<long long int>[]> atomic_avr_x, atomic_avr_y;
  int size;
  double atomic_avr_xy_converter;
};

//
// Integrals_Per_Circle and Avr_Z_Per_Circle - overlap integrals and average
// x+iy coordinates accumulated by one macroparticle during groups of many
// turns (NO_BB, ADIABATIC, STABILIZATION, BB). Ie. the kicked bunch density
// is represented in this case by a delta-function concentrated at only one
// particle which makes, however, many betatron oscillations (many turns).
//
template <class T>
// T = double for Integrals_Per_Circle, complex<double> for Avr_Z_Per_Circle
struct Per_Circle_Base : protected array<vector<T>, N::PHASES> { // [phase][circle]
  void resize(const N& n);
  void resize(int N_circles);
  bool empty();
  void add(int phase, int circle, T add);
  void normalize_by_n_turns(int N_turns, int phase, int circle); // parallelized, so per circle
  array<T, N::PHASES> weighted_average(const vector<double>& circle_weights);
  void write(ostream& os, const string& prefix);
  static string config_name();
};
struct Integrals_Per_Circle : public Per_Circle_Base<double> {
  static string config_name() { return "integrals.per.circle"; }
};
struct Avr_Z_Per_Circle : public Per_Circle_Base<complex<double> > {
  static string config_name() { return "centers.per.circle"; }
};

// Points store raw 4D coordinates for selected turns in
// Points[turn][circle][0/1 for x/y]
//
// This vector will be filled by all threads in parallel and when all finish,
// will be written to file. Unfortunately, direct writing to the same file by
// many threads would block them, this dictates storing large v_points in
// memory.
struct Points : protected vector<vector<array<complex<double>, 2> > > { // [turn][circle][0/1 for x/y]
  void resize(const N& n);
  void resize(int N_turns_stored,
	      int N_circles,
	      int N_selected_turns);
  using vector<vector<array<complex<double>, 2> > >::empty;
  void add(int turn, int circle, complex<double> xZ, complex<double> yZ);
  void write(ostream& os, const string& prefix);
  static string config_name() { return "points"; }
protected:
  int n_selected_turns;
};

//
// Create a vector of the classes above (denoted generically as class Data)
// and add two boolean flags: whether the corresponding data should be written
// out (bool write_i, per interaction point) and, if not, should it still be
// processed (bool always_fill). The latter is intended for
// Integrals_Per_Circle class which is used to calculate the final beam-beam
// corrections and, therefore, should be processed in any case.
//
template<class Data>
struct Output {
  struct Data_Write_It : public Data {
    bool write_it;
    Data_Write_It() : write_it(false) {}
  };
  vector<Data_Write_It> data; // [ip]
  ogzstream file;
  static bool always_fill() { return false; }
  static string config_name() { return Data::config_name(); }
};
// partial specialization for Integrals_Per_Circle class
template<> inline bool Output<Integrals_Per_Circle>::always_fill() { return true; }

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
// get<Integrals_Per_Circle>(ip)
//
// It returns the data of Integrals_Per_Circle type for the interaction point
// "ip" (counting from zero). Integrals_Per_Circle should be a unique class
// name in the tuple (otherwise, if the tuple contains two structures of the
// same type, the compiler will produce an error, see the documentation on
// std::get(std::tuple) which is used internally).
//
// To check whether the particular class and the corresponding algorithm has
// been chosen in the configuration, one can use
//
// is_active<Integrals_Per_Circles>(ip)
//
// "ip" is needed here since the algorithms are configurable per interaction
// point.
//
typedef tuple<Output<Integrals>,
	      Output<Avr_Z>,
	      Output<Integrals_Per_Circle>,
	      Output<Avr_Z_Per_Circle>,
	      Output<Points> > Outputs_tuple;
struct Outputs : public Outputs_tuple {
  Outputs(const N& n, Config& c, string output_dir,
	  double max_xy_radius); // for double <-> llong int converter in Avr_Z
  void write(int step, const vector<complex<double> >& z_kicker);
  // access to tuple individual components
  template<class T> T&   get      (int ip) { return std::get<Output<T> >(*this).data[ip]; }
  // should we process the data of type "T"?
  template<class T> bool is_active(int ip) { return !std::get<Output<T> >(*this).data[ip].empty(); }
};

// ------------------------------------------------------------
//
// Per_Circle_Base<T> = array<vector<T>, N::PHASES>  [phase][circle]
//
template <class T>
void Per_Circle_Base<T>::resize(const N& n) { resize(n.points()); }
template <class T>
void Per_Circle_Base<T>::resize(int N_circles) {
  for (int phase=0; phase<N::PHASES; ++phase) (*this)[phase].resize(N_circles);
}
template <class T>
bool Per_Circle_Base<T>::empty() { return (*this)[0].empty(); }
template <class T>
// initialized by zeros as any vector of doubles
void Per_Circle_Base<T>::add(int phase, int circle, T add) {
  (*this)[phase][circle] += add;
}
template <class T>
void Per_Circle_Base<T>::normalize_by_n_turns(int N_turns, int phase, int circle) { // parallelized, so per circle
  (*this)[phase][circle] /= N_turns;
}
template <class T>
array<T, N::PHASES> Per_Circle_Base<T>::weighted_average(const vector<double>& circle_weights) {
  array<T, N::PHASES> avr;
  fill(avr.begin(), avr.end(), 0.);
  for (int phase=0; phase<N::PHASES; ++phase) {
    const vector<T>& v = (*this)[phase];
    T& a = avr[phase];
    for (size_t i=0; i<v.size(); ++i) a += circle_weights[i] * v[i];
  }
  return avr;
}
template <class T>
void Per_Circle_Base<T>::write(ostream& os, const string& prefix) {
  for (int phase=0; phase<N::PHASES; ++phase) {
    for (size_t i = 0; i < (*this)[phase].size(); ++i) {
      os << prefix << " "
	 << phase << " "
	 << i << " "
	 << (*this)[phase][i] << '\n';
    }
  }
}

#endif
