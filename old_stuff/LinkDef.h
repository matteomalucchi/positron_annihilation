#include <vector>
#ifdef __ROOTCLING__
#pragma link C++ class vector<vector <float> >+;
#pragma link C++ class vector<float>+;
#endif

//rootcling -v4 -f mydict.cxx  -rmf libmydict.rootmap -rml libmydict.so  LinkDef.h
// g++ -shared -fPIC -o libmydict.so mydict.cxx `root-config --cflags --libs`