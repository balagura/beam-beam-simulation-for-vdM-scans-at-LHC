// Author V. Balagura, balagura@cern.ch (07.01.2019)

#include <config.hh>
#include <fstream>

using namespace std;

template<class T> ostream& operator<<(ostream& os, const vector<T>& v) {
  for (typename vector<T>::const_iterator i = v.begin(); i != v.end(); ++i)
    os << *i << " ";
  return os;
}

int main(int agc, char** argv) {
  Config c;
  ifstream in(argv[1]);
  c.read(in);
  cout << ">>> Initial config content\n";
  cout << c;
  cout << "Is \"aaa\" defined? " << (c.defined("aaa") ? "yes" : "no") << endl; 
  cout << "Is \"N\" defined? " << (c.defined("N") ? "yes" : "no") << endl; 
  cout << ">>> Access config variables:\n";
  cout << "c.l(\"N\") = c(\"N\") = " << c.l("N") << " = " << c("N") << endl;
  cout << "c.d(\"x\") = c[\"x\"] = " << c.d("x") << " = " << c["x"] << endl;
  cout << "c.s(\"longName123\") = " << c.s("longName123") << endl;
  cout << "c.vl(\"vN\") = " << c.vl("vN") << endl;
  cout << "c.vd(\"vX\") = " << c.vd("vX") << endl;
  cout << "c.vs(\"vS\") = " << c.vs("vS") << endl;

  cout << ">>> Check conversions\n";
  
  cout << ">>> long int -> double "
       << c.l("lToD") << " -> " << c.d("lToD") << endl;
  cout << ">>> long int -> string "
       << c.l("lToS") << " -> " << c.s("lToS") << endl;
  cout << ">>> long int -> vector<long int> "
       << c.l("lToVL") << " -> " << c.vl("lToVL") << endl;
  cout << ">>> long int -> vector<double> "
       << c.l("lToVD") << " -> " << c.vd("lToVD") << endl;
  cout << ">>> long int -> vector<string> "
       << c.l("lToVS") << " -> " << c.vs("lToVS") << endl;

  cout << ">>> double -> vector<double> "
       << c.d("dToVD") << " -> " << c.vd("dToVD") << endl;
  cout << ">>> double -> string "
       << c.d("dToS") << " -> " << c.s("dToS") << endl;
  cout << ">>> double -> vector<string> "
       << c.d("dToVS") << " -> " << c.vs("dToVS") << endl;

  cout << ">>> string -> vector<string> "
       << c.s("sToVS") << " -> " << c.vs("sToVS") << endl;

  cout << ">>> vector<long int> -> vector<double> "
       << c.vl("vlToVD") << " -> " << c.vd("vlToVD") << endl;
  cout << ">>> vector<long int> -> vector<string> "
       << c.vl("vlToVS") << " -> " << c.vs("vlToVS") << endl;

  cout << ">>> vector<double> -> vector<string> "
       << c.vd("vdToVS") << " -> " << c.vs("vdToVS") << endl;

  cout << ">>> Final config content\n";
  cout << c;
  return 0;
}
