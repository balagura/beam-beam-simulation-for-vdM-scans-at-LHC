// Output data classes for B*B simulation of beam-beam effects in van der Meer
// scans at LHC.
//
// Author V. Balagura, balagura@cern.ch (Nov 2019)
//

#include <algorithm>
#include <numeric>
#include <output.hh>

using namespace std;

// -------------------- Output Data classes --------------------
//
// Integrals_Per_Turn accumulate overlap integrals per turn
//

void Integrals_Per_Turn::resize(const N& n) { resize(n.turns()); }
void Integrals_Per_Turn::resize(int N_turns) {
  if (empty()) {
    atomic_integ = unique_ptr<atomic<long long int>[] >(new atomic<long long int>[N_turns]);
    double_integ.resize(N_turns);
  }
  double_converted = false;
  for (int i=0; i<N_turns; ++i) {
    // side note in case one wants to construct vector<atomic>: contrary to other vectors,
    // atomic vectors in current C++ are not initialized by zeros:
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0883r0.pdf
    // Here, array created by new() is used instead, so initialization by
    // zeros is anyway necessary
    atomic_integ[i].store(0);
  }
}
bool Integrals_Per_Turn::empty() { return double_integ.empty(); }
void Integrals_Per_Turn::add(int i_turn, double add) {
  atomic_integ[i_turn] += (long long int)(add * LLONG_MAX);
}
const vector<double>& Integrals_Per_Turn::integ() {
  if (!double_converted) {
    // transform atomic long long int's to double's
    for (int i=0; i<double_integ.size(); ++i) double_integ[i] = double(atomic_integ[i]) / LLONG_MAX;
  }
  return double_integ;
}
void Integrals_Per_Turn::write(ostream& os, const string& prefix) {
  const vector<double>& v = integ();
  for (size_t i_turn = 0; i_turn < v.size(); ++i_turn) {
    os << prefix << " "
       << i_turn << " "
       << v[i_turn] << '\n';
  }
}

//
// Avr_XY_Per_Turn - accumulate x+iy complex coordinate of the kicked bunch center-of-mass
void Avr_XY_Per_Turn::resize(const N& n) { resize(n.turns()); }
void Avr_XY_Per_Turn::resize(int N_turns) {
  if (empty()) {
    atomic_avr_x = unique_ptr<atomic<long long int>[] >(new atomic<long long int>[N_turns]);
    atomic_avr_y = unique_ptr<atomic<long long int>[] >(new atomic<long long int>[N_turns]);
    size = N_turns;
  }
  for (size_t i=0; i<N_turns; ++i) {
    atomic_avr_x[i].store(0);
    atomic_avr_y[i].store(0);
  }
}
bool Avr_XY_Per_Turn::empty() { return size == 0; }
//
void Avr_XY_Per_Turn::init(double max_xy_radius, int n_points) {
  atomic_avr_xy_converter = LLONG_MAX / max(max_xy_radius, 1.) / 100 / n_points;
  // As the maximal x * w[i], y * w[i] weighted particle's coordinate at a
  // given turn take the maximal initially simulated radius (note, w[i]<1)
  // times 100 for safety (initial circles are distorted by beam-beam, so
  // the maximal distance from 0,0 can grow a little, but not by 100). Any
  // accumulated sum of n coordinates can not exceed this value times n
  // (very conservatively since the coordinates tend to compensate each
  // other in the sum, so that the overall sum gives the average bunch
  // center). Ie. (sum of n coordinates) * atomic_avr_xy_converter = (sum of
  // n coordinates) * (LLONG_MAX / max_xy_radius / 100 / n_points) <
  // LLONG_MAX. If max_xy_radius < 1 - set it to 1.
}
void Avr_XY_Per_Turn::add(int i_turn, complex<double> add) {
  atomic_avr_x[i_turn] += (long long int)(real(add) * atomic_avr_xy_converter);
  atomic_avr_y[i_turn] += (long long int)(imag(add) * atomic_avr_xy_converter);
}
void Avr_XY_Per_Turn::write(ostream& os, const string& prefix) {
  for (size_t i_turn = 0; i_turn < size; ++i_turn) {
    // transform atomic long long int's to double's
    complex<double> avr_z(atomic_avr_x[i_turn] / atomic_avr_xy_converter,
			  atomic_avr_y[i_turn] / atomic_avr_xy_converter);
    os << prefix << " "
       << i_turn << " "
       << real(avr_z) << " " << imag(avr_z) << '\n';
  }
}

// Store raw 4D coordinates for selected turns in
// Points[turn][particle][0/1 for x/y]
//
// This vector will be filled by all threads in parallel and when all finish,
// will be written to file. Unfortunately, direct writing to the same file by
// many threads would block them, this dictates storing large v_points in
// memory.
// Points = vector<vector<array<complex<double>, 2> > >  [turn][particle][0/1 for x/y]
void Points::resize(const N& n) { resize(n.turns_stored(), n.points(), n.select_turns()); }
void Points::resize(int N_turns_stored,
		    int N_particles,
		    int N_selected_turns) {
  // no need to initialize, resize() is sufficient
  vector<vector<array<complex<double>, 2> > >::resize(N_turns_stored,
						      vector<array<complex<double>, 2> >(N_particles));
  n_selected_turns = N_selected_turns;
}
void Points::add(int turn, int particle, complex<double> xZ, complex<double> yZ) {
  array<complex<double>, 2 >& a = (*this)[turn][particle];
  a[0] = xZ;
  a[1] = yZ;
}
void Points::write(ostream& os, const string& prefix) {
  for (int turn = 0; turn < size(); ++turn) {
    for (int particle = 0; particle < (*this)[turn].size(); ++particle) {
      array<complex<double>, 2 >& a = (*this)[turn][particle];
      os << prefix << " "
	 << n_selected_turns * (turn+1) - 1 << " "
	 << particle << " "
	// since the complex coordinated eg. in X-X' plane is defined as zX =
	// X - iX', ie. with the minus sign, to print out X' one needs to
	// invert the sign of imag(zX):	
	 << real(a[0]) << " " << -imag(a[0]) << " "
	 << real(a[1]) << " " << -imag(a[1]) << "\n";
    }
  }
}

Outputs::Outputs(const N& n, Config& c, string output_dir,
		 double max_xy_radius) { // for double <-> llong int converter in Avr_XY_Per_Turn
  map<string, pair<vector<int>, bool> > selected;
  // selected[name] = (ips, boolean), boolean should be then set to true if
  // "name" coincides with one of output Data config_name()'s, otherwise
  // wrong option in config file will be signalled
  for (int ip=0; ip<n.ip(); ++ip) {
    string key = "output." + to_string(ip+1); // config IP number starts from one
    if (c.defined(key)) {
      for (string s : c.vs(key)) {
	pair<vector<int>, bool>& p = selected[s];
	p.first.push_back(ip);
	p.second = false;
      }
    }
  }
  auto f = [&n, &selected, &output_dir](auto&& x) {
	     x.data.resize(n.ip());
	     //
	     string s = x.config_name();
	     auto it = selected.find(x.config_name());
	     if (it != selected.end()) {
	       it->second.second = true; // flag this name as matched
	       for (int ip: it->second.first) {
		 x.data[ip].write_it = true;
		 // initialize IPs specified in it->second.first == vector<int>
		 x.data[ip].resize(n);
	       }
	     }
	     if (x.always_fill()) { // ensure all IPs are resize'd()
	       for (int ip=0; ip<n.ip(); ++ip) {
		 if (x.data[ip].empty()) x.data[ip].resize(n);
	       }
	     }
	     
	     if (any_of(x.data.begin(), x.data.end(),
			[](auto& e) { return e.write_it; })) {
	       string file_name = x.config_name();
	       replace(file_name.begin(), file_name.end(), '.', '_');
	       x.file.open((output_dir + "/" + file_name + ".txt.gz").c_str());
	       x.file << scientific; // change to scientific format for doubles
	     }
	   };
  apply([&f](auto& ... x) { (..., f(x)); }, *((Outputs_tuple*)this));
  //
  for (const auto& s: selected) {
    if (!s.second.second) { // no match found with any config_name()
      cerr << "output option " << s.first << " does not exist. Terminating\n";
      exit(1);
    }
  }
  // Avr_XY_Per_Turn, unfortunately, needs a special initialization of
  // double<->atomic<llong int> converter
  for (int ip=0; ip<n.ip(); ++ip) get<Avr_XY_Per_Turn>(ip).init(max_xy_radius, n.points());
}
void Outputs::reset(const N& n) { // for double <-> llong int converter in Avr_XY_Per_Turn
  auto f = [&n](auto&& x) {
	     for (int ip=0; ip<n.ip(); ++ip)
	       if (x.data[ip].write_it || x.always_fill()) x.data[ip].resize(n);
	   };
  apply([&f](auto& ... x) { (..., f(x)); }, *((Outputs_tuple*)this));
}
void Outputs::write(int step) {
  auto f = [&step](auto&& x) {
	     for(size_t ip=0; ip<x.data.size(); ++ip) {
	       if (x.data[ip].write_it) {
		 string prefix =
		   to_string(step) + " " + to_string(ip) + " ";
		 x.data[ip].write(x.file, prefix);
	       }
	     }
	   };
  apply([&f](auto& ... x) { (..., f(x));}, *((Outputs_tuple*)this));
}
