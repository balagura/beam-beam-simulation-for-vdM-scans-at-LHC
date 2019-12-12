//
// B*B simulation of beam-beam effects in van der Meer scans at LHC.
// Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------
//
// Reads configuration from istream
//
#include "config.hh"
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>

using namespace std;

/*
  There are six types of configurable data:
  
     vector<long int>  vector<double>  vector<string>
            long int          double          string
	    
  In the configuration file there might be a line
  key = 1
  
  How to interpret it, what can be a type corresponding to "key"? The data "1"
  might be long, double, string, or even a vector of those, in total 6
  possibilities. In the beginning, the "minimal" possible type is assigned:
  long int.  If in the program the user refers to "1" later as
  config.d("key"), ie. as a double, the "key" will be moved from long int's to
  double's. If then config.vd("key") appears, it will be moved further to
  vector<double>'s. The conversions can be only in one direction, ie. to the
  right and up, not to the left or down. If "key" is "upgraded" to
  vector<double>, it can not be downgraded back to long int.

  To allow such up/right-conversions one needs C++ type manipulations inside
  the 3x2 matrix of types above. The idea is simple, but expressing it in C++
  is a pain due to a lack of effective meta-programming tools in the language.

  The following performs
  
  1) right-conversions with Right-converter<From, To> template-specialized
  class and with right_convert<From, To> helper function,
  
  2) up+right conversions with function-overloaded convert(From from, To)
  where To type is sent as an argument just for overloading but To's object is
  empty. In this way a compiler knows that To's object is not needed and can
  drop this dependency.
*/
//
// right conversion, includes identity conversion From=To needed for
// up-conversion From -> vector<From>
//
template <class To> struct Right_Converter {
  template <class From> To operator()(From a) { return a; }
};
template<> struct Right_Converter<string> {
  template <class From> string operator()(From a) { return to_string(a); }
  string operator()(const string& a) { return a; }
};
template<class From, class To> To right_convert(From a) { return Right_Converter<To>()(a); }
//
// The same can be done explicitly (but then it is more difficult to add
// further types):
//
//long int right_converter(long int      a, to_type<long int>) { return a; }
//double   right_converter(long int      a, to_type<double>  ) { return a; }
//string   right_converter(long int      a, to_type<string>  ) { return to_string(a); }
//double   right_converter(double        a, to_type<double>  ) { return a; }
//string   right_converter(double        a, to_type<string>  ) { return to_string(a); }
//string   right_converter(const string& a, to_type<string>  ) { return a; }
//
// up+right conversion
//
template<class From, class To> To convert(From from, To) {
  return right_convert<From, To>(from);
}
template<class From, class To> vector<To> convert(From from, vector<To>) {
  return vector<To>(1, right_convert<From, To>(from));
}
template<class From, class To> vector<To> convert(const vector<From>& from, vector<To>) {
  vector<To> v(from.size());
  transform(from.begin(), from.end(), v.begin(),
	    [](auto x) { return right_convert<From, To>(x); });
  return v;
}
//
// move data from map<string, From> to map<string, To>
//
template<class From, class To>
typename map<string, To>::iterator get_from(map<string, From>& from,
					    map<string, To>& to,
					    const string& name) {
    auto j = from.find(name);
    if (j == from.end()) return to.end();
    auto i = to.insert(make_pair(name, convert(j->second, To()))).first;
    from.erase(j);
    return i;
}
//
// use get_from(from, to, name) explicitly
//
const vector<string>& Config::vs(const string& name) {
  // the longest search: first in vectors, then scalars of strings, then in
  // vectors and scalars of doubles and then in vectors+scalaras of integers
  auto& m = mvs;
  auto i = m.find(name);
  if (i == m.end()) i = get_from(ms,  m, name);
  if (i == m.end()) i = get_from(mvd, m, name);
  if (i == m.end()) i = get_from(md,  m, name);
  if (i == m.end()) i = get_from(mvl, m, name);
  if (i == m.end()) i = get_from(ml,  m, name);
  if (i == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else return i->second;
}
const vector<double>& Config::vd(const string& name) {
  auto& m = mvd;
  auto i = m.find(name);
  if (i == m.end()) i = get_from(md,  m, name);
  if (i == m.end()) i = get_from(mvl, m, name);
  if (i == m.end()) i = get_from(ml,  m, name);
  if (i == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else return i->second;
}
const vector<long int>& Config::vl(const string& name) {
  auto& m = mvl;
  auto i = m.find(name);
  if (i == m.end()) i = get_from(ml, m, name);
  if (i == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else return i->second;
}
const string& Config::s(const string& name) {
  auto& m = ms;
  auto i = m.find(name);
  if (i == m.end()) i = get_from(md, m, name);
  if (i == m.end()) i = get_from(ml, m, name);
  if (i == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else return i->second;
}
double Config::d(const string& name) {
  auto& m = md;
  auto i = m.find(name);
  if (i == m.end()) i = get_from(ml, m, name);
  if (i == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else return i->second;
}
long int Config::l(const string& name) const {
  auto& m = ml;
  auto i = m.find(name);
  if (i == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else return i->second;
}
bool Config::defined(const std::string& name) {
  return
    ml .find(name) != ml.end()  ||
    md .find(name) != md.end()  ||
    ms .find(name) != ms.end()  ||
    mvl.find(name) != mvl.end() ||
    mvd.find(name) != mvd.end() ||
    mvs.find(name) != mvs.end();
}

void Config::read(istream& is) {
  string line;
  while (getline(is, line)) {
    istringstream iss(line);
    iss >> ws;
    int c = iss.peek();
    if (c == '#' || c == EOF) continue; // comments start from #, empty lines are ignored
    string key;
    if (!getline(iss, key, '=')) {
      cerr << "Line \"" << line << "\" is not recognized\n";
      exit(EXIT_FAILURE);
    }
    { // remove trailing white space
      size_t pos = key.find_first_of(" \t");
      if (pos != string::npos) key.erase(pos);
    }
    vector<string> items;
    {
      string item;
      while (iss >> item) items.push_back(item);
    }
    if (items.empty()) {
      // just ignore it
      return;
    }      
    vkeys.push_back(key);
    char* rest;
    vector<long int> ints(items.size());
    // try to convert all items to long integers, if at least one failed - to doubles
    for (vector<string>::const_iterator it = items.begin(); it != items.end(); ++it) {
      *(ints.begin() + (it - items.begin())) = strtol(it->c_str(), &rest, 10);
      if (*rest) break; // conversion to long int failed
    }
    if (*rest == 0) {
      if (ints.size() == 1) {
	ml.insert(make_pair(key, ints[0]));
      } else {
	mvl.insert(make_pair(key, vector<long int>(ints)));
      }
    } else {
      vector<double> doubles(items.size());
      for (vector<string>::const_iterator it = items.begin(); it != items.end(); ++it) {
	*(doubles.begin() + (it - items.begin())) = strtod(it->c_str(), &rest);
	if (*rest) break; // conversion to double failed
      }
      if (*rest == 0) {
	if (doubles.size() == 1) {
	  md.insert(make_pair(key, doubles[0]));
	} else {
	  mvd.insert(make_pair(key, vector<double>(doubles)));
	}
      } else {
	if (items.size() == 1) {
	  ms.insert(make_pair(key, items[0]));
	} else {
	  mvs.insert(make_pair(key, vector<string>(items)));
	}
      }
    }
  }
}
template<class T> void Config::print(ostream& os, const map<string, T>& m) const {
  for (typename map<string, T>::const_iterator i = m.begin(); i!=m.end(); ++i)
    os << i->first << " = " << i->second << endl;
}
template<class T> void Config::print(ostream& os, const map<string, vector<T> >& m) const {
  for (typename map<string, vector<T> >::const_iterator i = m.begin(); i!=m.end(); ++i) {
    os << i->first << " =";
    for (typename vector<T>::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
      os << " " << *j;
    }
    os << endl;
  }
}    
ostream& operator<<(ostream& os, const Config& c) {
  os << "# ---------- Long integers: ----------\n";
  c.print(os, c.ml);
  os << "# ---------- Doubles: ----------\n";
  c.print(os, c.md);
  os << "# ---------- Strings: ----------\n";
  c.print(os, c.ms);
  os << "# ---------- Vectors of long integers: ----------\n";
  c.print(os, c.mvl);
  os << "# ---------- Vectors of doubles: ----------\n";
  c.print(os, c.mvd);
  os << "# ---------- Vectors of strings: ----------\n";
  c.print(os, c.mvs);
  return os;
}
