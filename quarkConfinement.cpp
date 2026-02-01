/*
Sat Mar 31 09:00:31 2012
This code simulates the fluxtube model of
A. Patel ( Nucl. Phys.
which includes the baryonic vertices with Nc links
attached at it. The model involves a cubic lattice with all the
links valued either +1 ( outward ),-1 (inward) or 0 (unoccupied).
Likewise one has quarks sitting at lattice sites valued again +1, -1
or 0. The baryonic vertices are also valued +1,-1 or 0.
However their distributed is dictated by the conservation
of Flux with is nothing but the Gauss law
Sum of link values at any site =
quark value at site - (Nc * Baryon value at site)*/

// first update: 22/02/2015 ( Zahidul Islam Jitu )

#include <boost/random.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;
using namespace boost;
/********************************************
Random Number Stuff - uses boost
********************************************/
mt19937 boost_rng(getpid() * time(NULL));
uniform_01<mt19937> number(boost_rng);
__inline double my_rand(void) { return (number()); }

// This function assigns the occupancy of links, sites ,....

double fill(void) {
  if (my_rand() > (2. / 3)) {
    return (1.);
  } else if (my_rand() < (1. / 3)) {
    return (-1.);
  } else
    return (0.);
}

/*************************************************
Parameters of the model
*************************************************/
// Number of Colors
const int Nc = 3;
// Lattice size
const int Size = 30;
// SampleSize
const int SampleSize = 8192;
const int updates = SampleSize / 4;
double E[SampleSize];
double beta;

/*******************************************
The variables
*******************************************/
double Link[Size][Size][Size][3]; // 3 directions on the lattice
double Q[Size][Size][Size];       // the quark variable
double B[Size][Size][Size];       // The baryon variable
double m, v;
// The mass of the quark and the baryonic vertex energy
int counter;
int Prev(int x) {
  int m;
  if (x != 0)
    m = x - 1;
  else {
    m = Size - 1;
    // Periodic boundary Condition
  }
  return (m);
}
/**************************************************
//
//
various Statistical averages using jackknife technology
//
**************************************************/
void multijack_average(double sample[], double *average, double *error) {
  // Given the array sample[ ] , we first calculate its average
  double sum = 0.0;
  // we sum up the samples
  for (int i = 0; i < SampleSize; i++) {
    sum += sample[i];
  }
  // Hence the average
  *average = sum / SampleSize;
  double bin_variance_squared, bin_variance, prev_variance = 0.0, jacked_sum;
  double jacked_average[SampleSize];
  double dev;
  int k;
  /*
  Make equal sized bins from the samples with z elements
  The bin size should be within the range 1 <= z <= ( SampleSize -1)
  but a quick glance at the formula
  sigma = \sqrt { \frac{k-1}{k} \sum^k_{i} ( x - x_i)^2
  shows that k >=2 to get non-zero variance. So it is meaningless to take
  binsizes greater than the half of the sample size.
  */
  for (int z = 1; z <= SampleSize / 2; z++) {
    k = SampleSize / z; // This is the number of bins produced
    bin_variance_squared = 0.0;
    // Now calculate the sample average after removing the z samples in the
    // m-th bin
    for (int m = 0; m < k; m++) {
      // To find the sum after removing the samples in the m-th bin
      // we simply remove the samples from the total sum calculated
      // previously
      jacked_sum = sum;
      for (int n = 0; n < z; n++) {
        jacked_sum -= sample[m * z + n];
      }
      // Thus the average after removing the m-th bin ( of size m ) is
      // given by
      jacked_average[m] = jacked_sum / (SampleSize - z);
      // Next we calculate the variance in the jackknifed samples by summing
      // the squares of the deviations
      dev = *average - jacked_average[m];
      bin_variance_squared += dev * dev;
    }
    // so the variance over this bin is
    bin_variance = sqrt((k - 1) * (bin_variance_squared) / k);
    // we are interested in the largest of the all these variances
    // which is found as follows
    if (bin_variance > prev_variance) {
      prev_variance = bin_variance;
    }
  }
  *error = prev_variance;
}
/**********************************************************
//
//
Functions related energy estimation
//
***************************************************************/
// Evaluate the size energy at the point (x,y,z)
double Site_Energy(int x, int y, int z) {
  double E = 0.;
  for (int d = 0; d < 3; d++) {
    E += Link[x][y][z][d] * Link[x][y][z][d];
  }
  return (E + m * (Q[x][y][z] * Q[x][y][z]) + v * (B[x][y][z] * B[x][y][z]));
}

double total_energy(void) {
  double energy = 0;
  for (int x = 0; x < Size; x++) {
    for (int y = 0; y < Size; y++) {
      for (int z = 0; z < Size; z++) {
        energy += Site_Energy(x, y, z);
      }
    }
  }
  return (energy);
}

int site_flux(int x, int y, int z) {
  int sum = Link[x][y][z][0] - Link[Prev(x)][y][z][0] + Link[x][y][z][1] -
            Link[x][Prev(y)][z][1];
  sum += Link[x][y][z][2] - Link[x][y][Prev(z)][2];
  return (sum);
}

// change of code, from void return type to bool type
bool place_vertex(int local_flux, int x, int y, int z, int *my_Q, int *my_B) {
  // we place quarks or vertices
  // quarks based on flux conservation
  // There are nine possibilities
  bool result = false; //(new code)this checks whether we should keep the new
                       //configuration based on flux conservation
  if (local_flux == 1) {
    *my_Q = 1;
    *my_B = 0;
    result = true;
  } else if (local_flux == -1) {
    *my_Q = -1;
    *my_B = 0;
    result = true;
  } else if (local_flux == 0) {
    *my_Q = 0;
    *my_B = 0;
    result = true;
  } else if (local_flux == Nc) {
    *my_Q = 0;
    *my_B = -1;
    result = true;
    // counter++;
  } else if (local_flux == Nc - 1) {
    *my_Q = -1;
    *my_B = -1;
    result = true;
    // counter--;
  } else if (local_flux == Nc + 1) {
    *my_B = -1;
    *my_Q = 1;
    result = true;
  } else if (local_flux == -Nc + 1) {
    *my_B = 1;
    *my_Q = 1;
    result = true;
  } else if (local_flux == -Nc - 1) {
    *my_B = 1;
    *my_Q = -1;
    result = true;
  } else if (local_flux == -Nc) {
    *my_B = 1;
    *my_Q = 0;
    result = true;
  }
  // if the gauss law not satisfied kill all the variable associated
  // with the site
  else {

    /* this comments marks the change of code*/
    // Link[x][y][z][0]=0;
    // Link[Prev(x)][y][z][0]=0;
    // Link[x][y][z][1]=0;
    // Link[x][Prev(y)][z][1]=0;
    // Link[x][y][z][2]=0;
    // Link[x][y][Prev(z)][2]=0;
    // my_Q=0;
    // my_B=0;
    result = false;

    /*end of change*/
  }
  return result;
}
// It is best to start off with the unfilled links and vertices ( quarks and
// Y-vertices ) as it is consistent with the Gauss law constraint at all
// points on the lattice
void coldstart(void) {
  for (int x = 0; x < Size; x++) {     // run along the x-axis
    for (int y = 0; y < Size; y++) {   // the same along y-axis
      for (int z = 0; z < Size; z++) { // and lastly along z-axis
        for (int d = 0; d < 3; d++) {
          Link[x][y][z][d] = 0; // no flux along any link
        }
        Q[x][y][z] = 0;
        B[x][y][z] = 0;
      }
    }
  }
  /*
  // We can place a quark antiquark pair ( Meson ) in the lattice
  //
  Q[0][Size/2][Size/2]=1;
  //
  Q[Size/2][Size/2][Size/2]=-1;
  //
  for(int x=0; x < Size/2; x++){
  //
  Link[x][Size/2][Size/2][0]=1;
  //
  }
  */
}
// The updates are performed by picking the values of the three Link variables
// which are attached with the vertex and then place quarks and Y-vertices
// in accordance with the Gauss law:
bool metropolis(double E, double Enext, double beta) {
  double prob = my_rand();
  double min;
  if (E >= Enext)
    min = 1.0;
  else
    min = exp(beta * (E - Enext));
  return (min > prob);
}

void local_update(int x, int y, int z) {
  int dLink[3]; // New Link Values
  int dFlux = 0;
  double Eflux = 0, Enew, Eold;
  // update the links first
  for (int k = 0; k < 3; k++) {
    dLink[k] = fill();
    dFlux += dLink[k];
    Eflux += dLink[k] * dLink[k];
  }
  int new_local_flux = dFlux - Link[Prev(x)][y][z][0] - Link[x][Prev(y)][z][1] -
                       Link[x][y][Prev(z)][2];
  int new_Q, new_B;
  bool consistency_check = false;
  consistency_check = place_vertex(new_local_flux, x, y, z, &new_Q,
                                   &new_B); //(new line)check whether the new
                                            //flux line is consistence with the
                                            // flux conservation
  if (consistency_check) { //(new code)only update if it's consistence
    Eold = Site_Energy(x, y, z);
    // Next calculate the energy shift due to this move and
    // acccept it in accordance with metropolis algorithm
    Enew = Eflux + m * new_Q * new_Q + v * new_B * new_B;
    // next check the viability of the move if accordance
    //  with the metropolis algorithm
    if (metropolis(Eold, Enew, beta)) {
      // update the config
      for (int j = 0; j < 3; j++) {
        Link[x][y][z][j] = dLink[j];
      }
      Q[x][y][z] = new_Q;
      B[x][y][z] = new_B;
    }
  }
}
// More Observables
void measure_vertices(double *qnum, double *qbarnum, double *bnum,
                      double *bbarnum) {
  *qnum = 0.0;
  *qbarnum = 0.0;
  *bnum = 0.0;
  *bbarnum = 0.0;
  for (int x = 0; x < Size; x++) {
    for (int y = 0; y < Size; y++) {
      for (int z = 0; z < Size; z++) {
        if (Q[x][y][z] == 1) {
          *qnum += 1.0;
        }
        if (Q[x][y][z] == -1) {
          *qbarnum += 1.0;
        }
        if (B[x][y][z] == 1) {
          *bnum += 1.0;
        }
        if (B[x][y][z] == -1) {
          *bbarnum += 1.0;
        }
      }
    }
  }
  // cout<<"baryon number= "<<*bnum<<endl;
}
// ofstream run("run10.data");
void sweep() {
  for (int x = 0; x < Size; x++) {
    for (int y = 0; y < Size; y++) {
      for (int z = 0; z < Size; z++) {
        local_update(x, y, z);
        // Check whether it is working
        // run << "x ="<< x << "
        // y=" << y << "
        // z=" << z << endl;
        // run <<"Link values ="<< Link[x][y][z][0] << " " << Link[x][y][z][1]
        // <<
        // " " << Link[x][y][z][2] <<endl;
        // run<< "Quark="<< Q[x][y][z]<< endl;
        // run << " Y vertex=" << B[x][y][z] << endl;
        // run << "Site Flux =" << site_flux(x,y,z)<<endl;
      }
    }
  }
}

int main(void) {
  double devEsq[SampleSize]; // Array for measuring deviation in Energy
  double quark[SampleSize], antiquark[SampleSize];
  double baryon[SampleSize], antibaryon[SampleSize];
  // Files for writing data
  ofstream energy("phasem_10new30.data");
  ofstream cv("CVm_10new30.data");
  ofstream Qfile("m_10quark30.data");
  ofstream Bfile("m_10baryon30.data");
  // Quark Mass
  m = 10;
  // Y-Vertex Energy
  v = 0.0;
  // Number of lattice sites
  int V = Size * Size * Size;
  double qnum, qbarnum, bnum, bbarnum;
  for (beta = 0.0; beta <= 10.0; beta = beta + 0.2) {
    // beta=2.0;
    qnum = 0.0;
    qbarnum = 0.0;
    bnum = 0.0;
    bbarnum = 0.0;
    coldstart();
    // First thermalize
    for (int j = 0; j < Size * Size * Size; j++) {
      // double output = total_energy();
      // cout << "This is original energy "<< output << endl;
      //
      counter = 0;
      sweep();
    }
    cout << "Thermalization done" << endl;
    for (int k = 0; k < SampleSize; k++) {
      sweep();
      // run << total_energy() << endl;
      //  measure the physical quantities
      E[k] = total_energy();
      measure_vertices(&qnum, &qbarnum, &bnum, &bbarnum);
      quark[k] = qnum;
      antiquark[k] = qbarnum;
      baryon[k] = bnum;
      antibaryon[k] = bbarnum;
      // cout<<k<<"	"<<baryon[k]<<"	"<<bnum<<endl;
    }
    double Eav, varE, cache;
    multijack_average(E, &Eav, &varE);
    cout << "beta=" << beta << ",Average Energy measured=" << Eav << endl;
    for (int j = 0; j < SampleSize; j++) {
      cache = beta * (E[j] - Eav);
      devEsq[j] = cache * cache / V;
    }
    double CV, CV_err;
    multijack_average(devEsq, &CV, &CV_err);
    double qav, qav_err, qbarav, qbarav_err;
    double bav, bav_err, bbarav, bbarav_err;
    multijack_average(quark, &qav, &qav_err);
    multijack_average(antiquark, &qbarav, &qbarav_err);
    multijack_average(baryon, &bav, &bav_err);
    multijack_average(antibaryon, &bbarav, &bbarav_err);
    energy << beta << " " << Eav / V << " " << varE / V << endl;
    cv << beta << " " << CV << " " << CV_err << endl;
    Qfile << beta << " " << qav << " " << qav_err << " " << qbarav << " "
          << qbarav_err << endl;
    Bfile << beta << " " << bav << " " << bav_err << " " << bbarav << " "
          << bbarav_err << endl;
  }
  return (0);
}
