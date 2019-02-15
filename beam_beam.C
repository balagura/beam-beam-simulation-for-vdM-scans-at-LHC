// Author V. Balagura, balagura@cern.ch (Jan 2019)

#include <vector>
#include <array>
#include <map>
#include <string>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <config.hh>
#include <boost/filesystem.hpp> // to mkdir output subdirectory
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp> // file_sink()
#include <boost/iostreams/filter/gzip.hpp> // gzip_compressor()

using namespace std;
using namespace std::complex_literals;

int main(int argc, char** argv) {
  // -------------------- Across 2D X-Y plane with weights --------------------
  // Our integral is 4-dimensional: over X,X',Y,Y'. But, one can integrate over
  // X-X', Y-Y' "circles" by averaging over turns, as one point "rotates" and
  // "swipes" those circles. So, one can integrate only over "circle" radii (rX,
  // rY). This is done below using 2-dimensional rX, rY Gaussian distribution,
  // it is covered by the grid of points inside the circle with R = 5 sigma. The
  // 2D Gaussian weight is assigned to every point, then the integral is
  // calculated as the sum of (2nd bunch densities * weight) over the sample of
  // "particles" starting from this grid with random phases and "turning" many
  // times.
  string config_file;
  if (argc != 2) {
    cout << "Usage: <program name> <configuration file>\n";
    return 0;
  }
  ifstream conf_file(argv[1]);
  Config c(conf_file);
  boost::filesystem::path output_dir;
  if (c.defined("output.directory")) {
    output_dir = c.s("output.directory");
  } else {
    string inp(argv[1]); // drop .txt from the end of configuration file
    size_t n = inp.size();
    if (inp.substr(n-4) == ".txt") {
      output_dir = inp.substr(0, n-4);
    } else {
      cerr << "The name of the configuration file must have the "
	   << "form *.txt (" << inp << " is given)\n";
      return 1;
    }
  }
  if (boost::filesystem::exists(output_dir)) {
    cerr << "Output directory " << output_dir << " already exists\n";
    return 1;
  }
  if (!boost::filesystem::create_directories(output_dir)) {
    cerr << "Can not create output directory " << output_dir << endl;
    return 1;
  }
  cout << "Output directory: " << output_dir << endl;
  double alpha = 1. / 137.035;
  double hbar = 0.197327e-15; // in Gev * m
  double beta0 = 1;
  
  // Qx,y are normally given with 2-3 digits after comma. So, after 1000 turns
  // the points return to their original positions (as 1000*Qx,y becomes
  // integer). To avoid this and add extra randomness, add "negligible"
  // irrational numbers exp(-9) = 0.0001234098 and sqrt(2)/1e4 =
  // 0.0001414214. This also makes Qx/Qy irrational,and opens otherwise closed
  // Lissajous figures of betatron oscillations in X-Y plane.  In reality Qx,y
  // should be irrational.
  double Qx = c["Qx"] + sqrt(2) / 1e4;
  double Qy = c["Qy"] + exp(-9);
  complex<double> zQx = exp(2i * M_PI * Qx);
  complex<double> zQy = exp(2i * M_PI * Qy);

  int N_points = c("N.points");
  int N_points_in_step; // <= N_points
  int N_no_bb = c("N.no.beam.beam.turns");
  int N_turns = c("N.turns");
  double N1 = c["N1"];
  double N2 = c["N2"];
  int N_steps = c.vd("x2").size();
  vector<complex<double> > z2(N_steps);
  for (int step = 0; step < N_steps; ++step) {
    z2[step] = complex<double>(c.vd("x2")[step], c.vd("y2")[step]);
  }
  double sig1 = c["sig1.x"];
  double sig2 = c["sig2.x"];
  double sig1_sq = sig1 * sig1;
  double sig2_sq = sig2 * sig2;
  double two_sig1_sq = 2 * sig1_sq;
  double two_sig2_sq = 2 * sig2_sq;
  double sig_sq_limit = 25 * sig2_sq;
  double k_int = 0.5 / M_PI / sig2_sq;
  // kick
  double k = 2 * c("Z1") * c("Z2") * alpha * hbar * N2 /
    beta0 / c["p"] * c["beta"] * 1e12; // in um^2
  bool is_kick_const = c.s("constant_kick") == "yes";
  // print Config when all type conversions are resolved
  cout << c;
  // store xZ, yZ once per "select_turns" turns
  int select_turns = c("select.one.turn.out.of");
  int N_turns_stored = N_turns / select_turns;
  struct Summary {
    vector<double> integ;
    vector<complex<double> > avr_z;
    double int0_analytic, int0, int0_to_analytic, int0_rel_err,
      int_to_int0_correction;
    complex<double> z, analytic_kick, approx_analytic_z;      
    Summary(int n_turns) : integ(n_turns), avr_z(n_turns) {}
    Summary() {}
  };
  vector<Summary> summary(N_steps, Summary(N_turns));
  vector<complex<double> > z(N_points), kick(N_points);
  vector<double> r2(N_points), e(N_points);
  double vdm_sig_sq = sig1_sq + sig2_sq;
  double two_vdm_sig_sq = 2 * vdm_sig_sq;


  boost::iostreams::filtering_ostream rx_ry_weights, points;
  bool output_rx_ry_weights = find(c.vs("output").begin(),
				   c.vs("output").end(),
				   "rx.ry.weights") != c.vs("output").end();
  bool output_points = find(c.vs("output").begin(),
			    c.vs("output").end(),
			    "points") != c.vs("output").end();
  if (output_rx_ry_weights) {
    rx_ry_weights.push(boost::iostreams::gzip_compressor());
    string file_name = (output_dir / "rx_ry_weights.txt.gz").string();
    rx_ry_weights.push(boost::iostreams::file_sink(file_name));
  }
  if (output_points) {
    points.push(boost::iostreams::gzip_compressor());
    string file_name = (output_dir / "points.txt.gz").string();
    points.push(boost::iostreams::file_sink(file_name));
  }
  //
  // sqrt(N.points) will be used to sample rX and rY coordinates in the interval
  // [0, 5*sig]
  int N_sqrt = int(sqrt(N_points));
  double rXY_bin = 5 * sig1 / N_sqrt;
  vector<double> rXY_bins(N_sqrt);
  for (int i=0; i<N_sqrt; ++i) rXY_bins[i] = (i + 0.5) * rXY_bin;
  // only N_points_in_step <= N_points will be used in rx, ry, w:
  vector<double> rx(N_points), ry(N_points), w(N_points);
  double w_fact = rXY_bin / sig1_sq;
  w_fact *= w_fact;
  // main loop
  for (int step = 0; step < N_steps; ++step) {
    cout << "Step " << step+1 << endl;
    // initialization
    N_points_in_step = 0;
    {
      for (int ix=0; ix<N_sqrt; ++ix) {
	for (int iy=0; iy<N_sqrt; ++iy) {
	  complex<double> z1(rXY_bins[ix], rXY_bins[iy]);
	  // select only the intersection of two 5*sig circles around bunches
	  // from 5sig X 5sig rectangle:
	  if (norm(z1)            < sig_sq_limit &&
	      norm(z1 - z2[step]) < sig_sq_limit) {
	    rx[N_points_in_step] = real(z1);
	    ry[N_points_in_step] = imag(z1);
	    ++N_points_in_step;
	  }
	}
      }
    }
    // distribution of radius length for two-dimensional Gaussian:
    //
    // exp(-r^2/2/sig^2) * r / sig^2 dr.
    //
    // So, the weight for two two-dimensional Gaussians for dr = rX, rY grid bin:
    {
      for (int i=0; i<N_points_in_step; ++i)
	w[i] = exp(-(rx[i] * rx[i] + ry[i] * ry[i])
		   / two_sig1_sq) * rx[i] * ry[i] * w_fact;
      double w_sum = accumulate(w.begin(),
				w.begin() + N_points_in_step, 0.);
      for (int i = 0; i < N_points_in_step; ++i) w[i] = w[i] / w_sum;    
    }
    if (output_rx_ry_weights) {
      for (size_t p = 0; p < N_points_in_step; ++p) {
	rx_ry_weights << step << " "
		      << p << " "
		      << rx[p] << " "
		      << ry[p] << " "
		      << w[p]<< '\n';
      }
    }    
    // store selected points in a[ selected_turn ][ point ][ 0/1 for x/y]
    vector<vector<array<complex<double>, 2> > >
      a(N_turns_stored, vector<array<complex<double>, 2> >
	(N_points_in_step));
    // For every (X,Y) pair obtain X-X', Y-Y' 4-dimensional
    // coordinates by randomly rotating X and Y radii
    vector<complex<double> > xZ(N_points_in_step), yZ(N_points_in_step);
    if (c.defined("seed")) {
      srand48(c("seed"));
    } else {
      srand48(time(NULL));
    }
    cout << "First random number: " << drand48() << endl;
    for (int i = 0; i < N_points_in_step; ++i) {
      xZ[i] = rx[i] * exp(2i * M_PI * drand48());
      yZ[i] = ry[i] * exp(2i * M_PI * drand48());
    }
    complex<double> kick_const;
    if (is_kick_const) {
      double z2_norm = norm(z2[step]);
      if (z2_norm == 0.) {
	kick_const = 0;
      } else {
	kick_const = k * (-z2[step]) / z2_norm * (1 - exp(-z2_norm / (two_sig2_sq + two_sig1_sq)));
      }
    }
    for (int i_turn = 0; ; ++i_turn) {
      if ((i_turn + 1) % select_turns == 0) {
	// last i_turn = N_turn-1 is written, first i_turn=0 - not
	// (assuming N_turn is a multiple of select_runs)
	int j_turn = (i_turn + 1) / select_turns - 1;
	// last j_turn = N_turn/select_runs - 1
	for (int point = 0; point < N_points_in_step; ++point) {
	  a[j_turn][point][0] = xZ[point];
	  a[j_turn][point][1] = yZ[point];
	}
	cout << "  Turn " << i_turn+1 << endl;
      }
      // always store integ, avr_z
      // summary[step].avr_z[i_turn] are .integ[i_turn] are initialized by zeros
      // as any vector of doubles
      for (int i = 0; i < N_points_in_step; ++i) {
	z[i] = complex<double>( real( xZ[i] ), real( yZ[i] ));
	summary[step].avr_z[i_turn] += z[i] * w[i];
	z[i] = z[i] - z2[step]; // w.r.t. the center of the 2nd bunch
	r2[i] = norm(z[i]);
	e[i] = exp(- r2[i] / two_sig2_sq);
	// overlap integral is calculated as a sum of (2nd bunch profile
	// density * weight) over the sample of 1st bunch "particles"
	summary[step].integ[i_turn] += e[i] * w[i];
      }
      summary[step].integ[i_turn] *= k_int;
      if (i_turn >= N_turns-1) break;
      //
      if (i_turn < N_no_bb) {
	// first N_no_bb times turn without kick to measure the initial,
	// possibly biased, integral; then calculate how much the kick
	// changes it. Taking the ratio to the initial integral (instead
	// of the exact integral) should cancel the bias at least
	// partially and improve precision.
	for (size_t i = 0; i < N_points_in_step; ++i) {
	  xZ[i] *= zQx;
	  yZ[i] *= zQy;
	}
      } else if (is_kick_const) {
	for (size_t i = 0; i < N_points_in_step; ++i) {
	  // The momentum kick is subtracted below because of the minus sign in
	  // the definition of eg. zX = X - iX'. The kick is added to X',
	  // so that i*kick is subtracted from zX.
	  xZ[i] = (xZ[i] - 1i * real( kick_const )) * zQx;
	  yZ[i] = (yZ[i] - 1i * imag( kick_const )) * zQy;
	}
      } else {
	for (size_t i = 0; i < N_points_in_step; ++i) {
	  if (r2[i] == 0.) {
	    xZ[i] *= zQx;  // kick = 0
	    yZ[i] *= zQy;
	  } else {
	    kick[i] = k * z[i] / r2[i]  * (1 - e[i]);
	    // The momentum kick is subtracted below because of the minus sign in
	    // the definition of eg. zX = X - iX'. The kick is added to X',
	    // so that i*kick is subtracted from zX.
	    xZ[i] = (xZ[i] - 1i * real( kick[i] )) * zQx;
	    yZ[i] = (yZ[i] - 1i * imag( kick[i] )) * zQy;
	  }
	}
      }
    }
    if (output_points) {
      const char* xyC = "xy";
      for (size_t xy = 0; xy < 2; ++xy) {
	for (size_t j_turn = 0; j_turn < N_turns_stored; ++j_turn) {
	  for (size_t point = 0; point < N_points_in_step; ++point) {
	    points << step << " " << xyC[xy] << " "
		   << j_turn * select_turns << " "
		   << point << " "
		   << real(a[j_turn][point][xy]) << " "
		   << imag(a[j_turn][point][xy]) << '\n';
	  }
	}
      }
    }
    summary[step].int0_analytic = 1. / M_PI / two_vdm_sig_sq *
      exp(-norm(z2[step]) / two_vdm_sig_sq);
    summary[step].int0 = accumulate(summary[step].integ.begin(),
				    summary[step].integ.begin()+N_no_bb, 0.) /
      double(N_no_bb);
    summary[step].int0_to_analytic =
      summary[step].int0 / summary[step].int0_analytic;
    double sd_int0 = 0;
    for (auto i = summary[step].integ.begin();
	 i != summary[step].integ.begin() + N_no_bb; ++i) {
      double dx = *i - summary[step].int0;
      sd_int0 += dx * dx;
    }    
    summary[step].int0_rel_err =
      sqrt(sd_int0 / (N_no_bb - 1.) / double(N_no_bb)) / summary[step].int0;
    // assume that after 2000 turns the distribution stabilizes
    summary[step].int_to_int0_correction =
      accumulate(summary[step].integ.begin() + N_no_bb + 2000,
		 summary[step].integ.end(), 0.)
      / (summary[step].integ.size() - 2000 - N_no_bb) / summary[step].int0;
    summary[step].z =
      accumulate(summary[step].avr_z.begin() + N_no_bb + 2000,
		 summary[step].avr_z.end(), complex<double>(0, 0))
      / complex<double>(summary[step].avr_z.size() - 2000 - N_no_bb, 0);
    {
      double r2 = norm( z2[step] );
      if (r2 == 0) {
	summary[step].analytic_kick = 0;
      } else {
	// the sign of this kick is not inverted (to print it out in summary),
	// ie. it is directly for X', not for z which is defined with minus as
	// X - iX';
	// minus sign is since the kick is opposite to the direction to z2
	summary[step].analytic_kick = - k / 2 * z2[step] / r2 * (1 - exp(-r2 / two_vdm_sig_sq));
      }
    }
    // estimation as Re(-i*kick / (1 - exp(2*pi*Q*i))), minus acording to our
    // z=X-iX':
    summary[step].approx_analytic_z =
      complex<double>(real(summary[step].analytic_kick) / tan(M_PI * Qx),
		      imag(summary[step].analytic_kick) / tan(M_PI * Qy));
  }
  cout << "Final results are written to " << output_dir << " directory\n"
       << "The summary appears in \"summary.txt\".\n"
       << "The 4th column \"int_to_int0_correction\" in this file contains the\n"
       << "final correction, i.e. the ratio of the luminosities with and\n"
       << "without beam-beam electromagnetic interaction, for the scan step and\n"
       << "the corresponding beam offset specified in the first three columns.\n"
       << "The columns 5-8 refer to no-beam-beam case: : analytic and numeric\n"
       << "integrals and their ratio (\"int0_analytitc\", \"int0\" and\n"
       << "\"int0_to_analytic\"), and the estimated statistical error of\n"
       << "\"int0\" (\"int0_err\").\n"
       << "The (9,10) column pair gives the numerically calculated x-y center of\n"
       << "of the first bunch, while (11,12) and (13,14) pairs give the\n"
       << "analitically calculated shifts, the first in the coordinates\n"
       << "beta*(dx/dz, dy/dz) (which should be precise) and the second -\n"
       << "an estimation in x,y under the assumption that its ratio to the\n"
       << "first pair is -1/tan(pi*Q).\n";
  {
    string file_name = (output_dir / "config.txt").string();
    ofstream config_out(file_name);
    config_out << c;    
  }
  {
    string file_name = (output_dir / "summary.txt").string();
    ofstream fsummary(file_name);
    for (int step = 0; step < N_steps; ++step) {
      fsummary << step << " "
	       << real(z2[step]) << " "
	       << imag(z2[step]) << " "
	       << summary[step].int_to_int0_correction << " "
	       << summary[step].int0_analytic << " "
	       << summary[step].int0 << " "
	       << summary[step].int0_to_analytic << " "
	       << summary[step].int0_rel_err << " "
 	       << real(summary[step].z) << " "
	       << imag(summary[step].z) << " "
	       << real(summary[step].analytic_kick) << " "
	       << imag(summary[step].analytic_kick) << " "
	       << real(summary[step].approx_analytic_z) << " "
	       << imag(summary[step].approx_analytic_z) << "\n";
    }
  }
  if (find(c.vs("output").begin(),
	   c.vs("output").end(), "integrals") != c.vs("output").end()) {
    boost::iostreams::filtering_ostream integrals;
    integrals.push(boost::iostreams::gzip_compressor());
    string file_name = (output_dir / "integrals.txt.gz").string();
    integrals.push(boost::iostreams::file_sink(file_name));
    for (int step = 0; step < N_steps; ++step) {
      for (size_t turn = 0; turn < summary[step].integ.size(); ++turn) {
	integrals << step << " "
		  << real(z2[step]) << " "
		  << imag(z2[step]) << " "
		  << turn << " "
		  << summary[step].integ[turn] << '\n';
      }
    }
  }
  if (find(c.vs("output").begin(),
	   c.vs("output").end(), "centers") != c.vs("output").end()) {
    boost::iostreams::filtering_ostream centers;
    centers.push(boost::iostreams::gzip_compressor());
    string file_name = (output_dir / "centers.txt.gz").string();
    centers.push(boost::iostreams::file_sink(file_name));
    for (int step = 0; step < N_steps; ++step) {
      for (size_t turn = 0; turn < summary[step].avr_z.size(); ++turn) {
	centers << step << " "
		<< real(z2[step]) << " "
		<< imag(z2[step]) << " "
		<< turn << " "
		<< real(summary[step].avr_z[turn]) << " "
		<< imag(summary[step].avr_z[turn]) << '\n';
      }
    }
  }
  return(0);
}
