#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;

int main() {

  ifstream myfile("baryonData.txt");
  ofstream myWfile("baryonDataDist.txt");
  int class_number = 0;
  double class_interval = 0.0;
  double k;
  std::vector<double> dist;
  vector<int> distCount;
  vector<double> baryon_number;
  if (myfile.is_open()) {

    while (!myfile.eof()) {
      myfile >> k;
      baryon_number.push_back(k);
    }
    myfile.close();
  }

  // cout<<baryon_number.size()<<endl;
  cout << "input class number= ";
  cin >> class_number;
  cout << endl;

  sort(baryon_number.begin(), baryon_number.end());
  // cout<<"I am okay!"<<endl;

  class_interval =
      (double)(baryon_number[baryon_number.size() - 1] - baryon_number[0]) /
      (double)class_number;

  for (int i = 0; i < class_number; i++) {
    dist.push_back((double)baryon_number[0] +
                   (double)(i + 0.5) * class_interval);
    distCount.push_back(0);
  }

  for (int i = 0; i < baryon_number.size(); i++) {
    int counter = floor((double)(baryon_number[i] - baryon_number[0]) /
                        (double)class_interval);
    distCount[counter]++;
  }
  for (int i = 0; i < dist.size(); i++) {
    myWfile << dist[i] << " " << distCount[i] << endl;
  }
  return 0;
}
