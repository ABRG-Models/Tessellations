#include <morph/tools.h>
#include <morph/HdfData.h>
#include <morph/Random.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#include <hdf5.h>
#include <unistd.h>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <cctype>
// #include <boost/math/special_functions/bessel.hpp>
#define PI 3.1415926535897932
#define NUMPOINTS 41 //just the A-E rows.
//#define NUMPOINTS 79 //just the A-E rows.

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HdfData;
using morph::Tools;
using namespace std;

/*!
 * For mixing up bits of three args; used to generate a good random
 * seed using time() getpid() and clock().
 */
unsigned int
mix (unsigned int a, unsigned int b, unsigned int c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}



int main (int argc, char **argv)
{
    vector <pair <float, float>> centres; //seed points for regions
    centres.resize(NUMPOINTS);
    std::string num;
    unsigned int off;
    if (argc > 2) {
    	off = atoi(argv[2]);
    }
    else {
        num = "";
	off = 1;
    }

    double maxX = stod(argv[1]);
    double minX = -maxX;
    double maxY = maxX;
	double minY = -maxY;

    ofstream afile ( "./centres.h");
    ofstream bfile ( "./centres.inp");


   //   unsigned int seed = off;
   unsigned int seed = mix(clock(), time(NULL), getpid());
   cout << "seed " << seed << endl;

    cout << "numpoints " << NUMPOINTS << " maxX " << maxX << " minX " << minX << " maxY " << maxY << " minY " << minY << endl;
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);
    cout << "Random float number is " << ruf.get() << endl;
    // You can find the min and max:
    cout << "That float RNG has min and max: " << ruf.min() << "/" << ruf.max() << endl;

    int count = 0;
    while (count < NUMPOINTS) {
       double choice = ruf.get();
       if ((0 < choice) && (choice <= 0.25)) {
           double x = ruf.get() * maxX;
           double y  = ruf.get() * maxY;
           cout << " radius " << x*x + y*y << endl;
           if (x*x + y*y < maxX*maxX) {
               centres[count].first = x;
               centres[count].second = y  ;
               cout << " if 1 " << endl;
               count++;
	       }
       }
       if ((0.25 < choice) && (choice <= 0.5)) {
           double x = -ruf.get() * maxX;
           double y  = ruf.get() * maxY;
           cout << " radius " << x*x + y*y << endl;
           if (x*x + y*y < maxX*maxX) {
               centres[count].first = x;
               centres[count].second = y;
               cout << " if 2 " << endl;
               count++;
	       }
       }
	   if ((0.5 < choice) && (choice  <= 0.75)) {
           double x = -ruf.get() * maxX;
           double y  = -ruf.get() * maxY;
           cout << " radius " << x*x + y*y << endl;
           if (x*x + y*y < maxX*maxX) {
               centres[count].first = x;
               centres[count].second = y  ;
               cout << " if 3 " << endl;
               count++;
	       }
       }
	   if ((0.75 < choice) && (choice  <= 1.0)) {
           double x = ruf.get() * maxX;
           double y  = -ruf.get() * maxY;

           cout << " radius " << x*x + y*y << endl;
           if (x*x + y*y < maxX*maxX){
               centres[count].first = x;
               centres[count].second = y;
               cout << " if 4 " << endl;
               count++;
	       }
       }
    } //end of setting of random value

    std::sort(centres.begin(),centres.end());

    for (int i=0; i < NUMPOINTS; i++) {
	  string sindex = to_string(i);
	  string centres1 = " centres[" + sindex + "].first = ";
	  string  centres2 = " centres[" + sindex + "].second = ";
      afile << centres1 << centres[i].first << " ; " << centres2 << centres[i].second << " ; " << endl;
      bfile << centres[i].first << " " << centres[i].second << endl;
    }

    return 0;
};
