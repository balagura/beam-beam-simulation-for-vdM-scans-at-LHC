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

   respectively. "l/d/s" function names stand for the long integer, double and
   string types, respectively. Alternatively and equivalently, long int and
   double can be accessed using shorter variants 
   c("key1")
   c["key1"]
   
   These are C++ operators () and [] and they are equaivalent to l() and d()
   member functions, respectively.
   
   In case of multiple X1 X2 X3 ... the "vector" versions with "v" prefix
   should be used:
   c.vl("key2")
   c.vd("key2")
   c.vs("key2")

   They return std::vector of the corresponding type (long int, double or
   string, respectively).

   Therefore, there are six types of configurable data:
   
     vector<long int>  vector<double>  vector<string>
            long int          double          string
	    
   In the configuration file there might be a line
   key = 1
  
   How to interpret it, what can be a type corresponding to "key"? The data
   "1" might be long, double, string, or even a vector of those, in total 6
   possibilities. In the beginning, the "minimal" possible type is assigned:
   long int.  If in the program the user refers to "1" later as c.d("key"),
   ie. as a double, the "key" will be moved from the container of long int's
   to double's. If then config.vd("key") appears, it will be moved further to
   vector<double>'s. In the table above the conversions can only be performed
   in the directions to the right and from the bottom row to the
   top. Conversions to the left or down are forbidden, eg. if "key" is
   "upgraded" to vector<double>, it can not be downgraded back to long int.
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
  const std::string& s(const std::string& name);
  const std::vector<long int>&      vl(const std::string& name);
  const std::vector<double>&        vd(const std::string& name);
  const std::vector<std::string>&   vs(const std::string& name);
  bool defined(const std::string& name);
  // returns vector of all keys to which value(s) have been assigned:
  const std::vector<std::string>& keys() const { return vkeys; }
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
  std::vector<std::string> vkeys;
};

#endif
