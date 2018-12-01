#include <chrono>
#include <iostream>
#include <array>
#include <vector>

#include "kdtree.h"

class MyPoint : public std::array<double, 2>
{
public:
	static const int DIM = 2;

	MyPoint() {}
	MyPoint(double x, double y)
	{ 
		(*this)[0] = x;
		(*this)[1] = y;
	}
};


int main() {

  using namespace std::chrono;
  
  // set parameters
  const int width = 500;
  const int height = 500;
  const int resolution = 1000000;
  const int nquerys = 10000;
  
  // iterate over different sizes for a tree
  for (int npoints = 10; npoints <= 1000000; npoints *= 10)
    {
      // generate points
      std::vector<MyPoint> points(npoints);
      for (int i = 0; i < npoints; i++)
	{
	  const double x = ((double) (rand() % (resolution*width)))/resolution;
	  const double y = ((double) (rand() % (resolution*height)))/resolution;
	  points[i] = MyPoint(x, y);
	}

      // set query points
      std::vector<MyPoint> querys(nquerys);
      for (int i = 0; i < nquerys; i++)
	{
	  const double x = ((double) (rand() % (resolution*width)))/resolution;
	  const double y = ((double) (rand() % (resolution*height)))/resolution;
	  querys[i] = MyPoint(x, y);
	}
      const MyPoint query(0.5 * width, 0.5 * height);

      // place to return values to force evaluation
      std::vector<int> responses1(nquerys);
      std::vector<int> responses2(nquerys);

      // time some basic operations
      steady_clock::time_point t1 = steady_clock::now();
      
      kdt::KDTree<MyPoint> kdtree(points);
      
      steady_clock::time_point t2 = steady_clock::now();
      
      for (int i = 0; i < nquerys; i++)
      	{
	  responses1[i]=kdtree.nnSearch(querys[i]);
	}
      
      steady_clock::time_point t3 = steady_clock::now();
      
      kdtree.clear();
      
      steady_clock::time_point t4 = steady_clock::now();
      
      //an implementation of linear search for comparison/testing
      int nearestPoint = 0;
      double distToNearest = std::numeric_limits<double>::max();
      double nextDist;
      for (int i = 0; i < npoints; i++)
	{
	  nextDist=0;
	  for (int j = 0; j < MyPoint::DIM;j++)
	    {
	      nextDist+=(querys[0][j]-points[i][j])*(querys[0][j]-points[i][j]);
	    }
	  if (nextDist<distToNearest)
	    {
	      distToNearest=nextDist;
	      nearestPoint=i;
	    }
	}

      steady_clock::time_point t5 = steady_clock::now();

      // build by making insertions
      for (int i = 0; i < npoints; i++)
	{
	  kdtree.insert(points[i]);
	}

      steady_clock::time_point t6 = steady_clock::now();

      for (int i = 0; i < nquerys; i++)
	{
	  responses2[i]=kdtree.nnSearch(querys[i]);
	}

      steady_clock::time_point t7 = steady_clock::now();

      if (kdtree.validate())
	{
	  std::cout << "Insertion built tree was valid" << std::endl;
	}
      else
	{
	  std::cout << "Error: Insertion built tree was invalid" << std::endl;
	}
      kdtree.clear();
      
      // Force compiler to compute the above
      int force=0;
      for (int i = 0; i < nquerys; i++)
	{
	  force+=responses1[i];
	  force+=responses2[i];
	  force+=nearestPoint;
	}
      if (force==17)
	{
	  std::cout << "An Astounding Thing Happened!" << std::endl;
	}
      
      // Display the times to cout
      duration<double> time_span;
      double avgtime;
      
      std::cout << "Times for " << npoints << " points:" << std::endl;
      
      time_span = duration_cast<duration<double>>(t2 - t1);
      std::cout << "Building:" << time_span.count() << " seconds." << std::endl;
      
      time_span = duration_cast<duration<double>>(t3 - t2);
      avgtime=time_span.count()/nquerys;
      std::cout << "nnSearch(average of " << nquerys << " runs):" << avgtime << " seconds." << std::endl;
      
      time_span = duration_cast<duration<double>>(t4 - t3);
      std::cout << "clear:" << time_span.count() << " seconds." << std::endl;
      
      time_span = duration_cast<duration<double>>(t5 - t4);
      std::cout << "linear nnSearch:" << time_span.count() << " seconds." << std::endl;

      std::cout << "Built by insertion case:" << std::endl;

      time_span = duration_cast<duration<double>>(t6 - t5);
      std::cout << "Building:" << time_span.count() << " seconds." << std::endl;

      time_span = duration_cast<duration<double>>(t7 - t6);
      avgtime=time_span.count()/nquerys;
      std::cout << "nnSearch(average of " << nquerys << " runs):" << avgtime << " seconds." << std::endl;
      
      // Test nnSearch against the linear search result
      if (nearestPoint!=responses1[0])
	{
	  std::cout << "The built in one tree returned an incorrect nearest neighbour." << std::endl;
	}

      if (nearestPoint!=responses2[0])
	{
	  std::cout << "The built by insertion tree returned an incorrect nearest neighbour." << std::endl;
	}

      std::cout << std::endl;
    }
}
