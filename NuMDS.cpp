#include "NuMDS.h"
#include <sstream>

int main(int argc, char *argv[]) {
  stringstream ss;
  ss << argv[3];
  ss >> seed;
  ss.clear();
  ss << argv[2];
  ss >> cutoff_time;
  ss.clear();
  para1 = 60;
  file = argv[1];
  if (BuildInstance(argv[1]) != 0) {
    cerr << "Open instance file failed." << endl;
    return 1;
  }


  if (seed < 0U || seed > ~0U) {
    seed = 10;
  }

  if (cutoff_time < 0 || cutoff_time > (int)(~0U >> 1)) {
    cutoff_time = 1000;
  }

  srand(seed);

  start = chrono::steady_clock::now();
  reduction();

  start = chrono::steady_clock::now();
  ConstructDS();
  if (candidate_size <= 1)
  {
	  cout << argv[1] << " " << seed << " " << fix_must_in.size() << "\n";
	  return 0;
  }
  LocalSearch();
  int number1 = 0;
  int number2 = 0;
  for (int i = 1; i <= v_num; i++) {
	  if (record_in[i] == 0) {
		  number1++;
	  }
	  if (record_in[i] == 1) {
		  number2++;
	  }
  }

  cout << file << " " << seed << " " << best_c_size << endl;


  FreeMemory();

  return 0;
}
