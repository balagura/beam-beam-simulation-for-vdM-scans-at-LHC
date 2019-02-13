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
  cout << ">>> Access zToDouble as double: " << c.d("zToDouble") << endl;
  cout << ">>> Access vZToDouble as a vector of doubles: "
       << c.vd("vZToDouble") << endl;
  cout << ">>> Access iToVDouble as a vector of doubles: "
       << c.vd("iToVDouble") << endl;
  cout << ">>> Access dToVDouble as a vector of doubles: "
       << c.vd("dToVDouble") << endl;
  cout << ">>> Access iToVInt as a vector of doubles: "
       << c.vl("iToVInt") << endl;
  cout << ">>> Final config content\n";
  cout << c;
  return 0;
}
