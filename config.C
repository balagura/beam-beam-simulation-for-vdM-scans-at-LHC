// Author V. Balagura, balagura@cern.ch (07.01.2019)

#include <config.hh>
#include <sstream>
#include <cstdlib>
#include <cstdio>

using namespace std;

long int Config::l(const string& name) const {
  return Config::find_or_fail(ml, name)->second;
}
double Config::d(const string& name) {
  // search first in doubles and then also in integers
  map<string, double>::const_iterator i = md.find(name);
  if (i != md.end()) {
    return i->second;
  } else {
    map<string, long int>::iterator j = Config::find_or_fail(ml, name);
    // if found in integers: move it from integers to doubles and return
    i = md.insert(make_pair(name, double(j->second))).first;
    ml.erase(j);
    return i->second;
  }
}
const string& Config::s(const string& name) const {
  return Config::find_or_fail(ms, name)->second;
}
const vector<long int>& Config::vl(const string& name) {
  // search first in vectors and then in "scalars" of integers
  map<string, vector<long int> >::const_iterator i = mvl.find(name);
  if (i != mvl.end()) {
	return i->second;
  } else {
    map<string, long int>::iterator j = Config::find_or_fail(ml, name);
    i = mvl.insert(make_pair(name, vector<long int>(1, j->second))).first;
    ml.erase(j);
    return i->second;
  }
}
const vector<double>&   Config::vd(const string& name) {
  // search first in vectors+scalars of doubles and also in integers
  map<string, vector<double> >::const_iterator i = mvd.find(name);
  if (i != mvd.end()) {
    return i->second;
  } else {
    map<string, vector<long int> >::iterator j = mvl.find(name);
    if (j != mvl.end()) {
      // if found in integers: move it from integers to doubles and return
      vector<double> vd(j->second.size());
      copy(j->second.begin(), j->second.end(), vd.begin());
      i = mvd.insert(make_pair(name, vd)).first;
      mvl.erase(j);
      return i->second;
    } else {
      map<string, double>::iterator j1 = md.find(name);
      if (j1 != md.end()) {
	// if found in double scalars: move it to vector of doubles and return
	i = mvd.insert(make_pair(name, vector<double>(1, j1->second))).first;
	md.erase(j1);
	return i->second;
      } else {
	map<string, long int>::iterator j2 = Config::find_or_fail(ml, name);
	i = mvd.insert(make_pair(name, vector<double>(1, j2->second))).first;
	ml.erase(j2);
	return i->second;
      }
    }
  }
}
const vector<string>& Config::vs(const string& name) const {
  return Config::find_or_fail(mvs, name)->second;
}
bool Config::defined(const std::string& name) {
  if (ml.find(name) != ml.end()) return true;
  if (md.find(name) != md.end()) return true;
  if (ms.find(name) != ms.end()) return true;
  if (mvl.find(name) != mvl.end()) return true;
  if (mvd.find(name) != mvd.end()) return true;
  if (mvs.find(name) != mvs.end()) return true;
  return false;
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

template<class T> typename map<string, T>::const_iterator
Config::find_or_fail(const map<string, T>& m, const string& name) {
  typename map<string, T>::const_iterator it = m.find(name);
  if (it == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else
    return it;
}    
template<class T> typename map<string, T>::iterator
Config::find_or_fail(map<string, T>& m, const string& name) {
  typename map<string, T>::iterator it = m.find(name);
  if (it == m.end()) {
    cerr << "Configuration parameter " << name << " is not defined\n";
    exit(EXIT_FAILURE);
  } else
    return it;
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
