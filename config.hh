#ifndef config_hh 
#define config_hh 1

/* "Config" class reads the configuration file in the format
   
   key1 = X
     or
   key2 = X1 X2 X3 ...

   where X, X1, X2, ... can be long integers, doubles or simply "not numbers"
   (strings). The syntax of the C++ program is the following:
   Config c;
   c.read(ifstream("file_name"));

   or
   Config c(ifstream("file_name"));
   
   The long integer, double or std::string parameters can then be accessed as
   c.l("key1")
   c.d("key1")
   c.s("key1")

   respectively, for a single X. " .l/d/s" stand for the long integer, double
   and string types. Alternatively and equivalently, long int and double can be
   accessed using shorter variants
   c("key1")
   c["key1"]

   Note that the circular () and square [] brackets are used for the long int
   and the double, respectively.
   In case of multiple X1 X2 X3 ...
   c.vl("key2")
   c.vd("key2")
   c.vs("key2")

   give std::vector of the corresponding type, ".v" stands for the vector.

   If one of the double-type parameters appears as a long integer in the
   configuration file, eg.
   key = 100

   it is initially recognized as an integer, but when accessed as a double:
   c.d("key")

   its type automtically changes to double, and it can not be accessed
   as an integer anymore. The same rule holds for a conversion of a vector of
   integers to a vector of doubles, and also for a conversion of single values
   to vectors: an integer to a vector of integers or doubles, and a double to a
   vector of doubles.
   
*/
// Author V. Balagura, balagura@cern.ch (07.01.2019)

#include <iostream>
#include <vector>
#include <map>
#include <string>

struct Config {
  Config() {}
  Config(std::istream& is) { read(is); }
  void read(std::istream& is);
  long int           l(const std::string& name) const;
  long int  operator()(const std::string& name) const { return l(name); };
  double             d(const std::string& name);
  double    operator[](const std::string& name)       { return d(name); }
  const std::string& s(const std::string& name) const;
  const std::vector<long int>&      vl(const std::string& name);
  const std::vector<double>&        vd(const std::string& name);
  const std::vector<std::string>&   vs(const std::string& name) const;
  bool defined(const std::string& name);
  friend std::ostream& operator<<(std::ostream& os, const Config& c);
private:
  template<class T> static typename std::map<std::string, T>::const_iterator
  find_or_fail(const std::map<std::string, T>& m, const std::string& name);
  template<class T> static typename std::map<std::string, T>::iterator
  find_or_fail(std::map<std::string, T>& m, const std::string& name);
  template<class T> void print(std::ostream& os,
			       const std::map<std::string, T>& m) const;
  template<class T> void print(std::ostream& os, const std::map<std::string,
			       std::vector<T> >& m) const;
  std::map<std::string, long int>    ml;
  std::map<std::string, double>      md;
  std::map<std::string, std::string> ms;
  std::map<std::string, std::vector<long int> >      mvl;
  std::map<std::string, std::vector<double> >        mvd;
  std::map<std::string, std::vector<std::string> >   mvs;
};

#endif
