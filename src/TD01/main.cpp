///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <iomanip> // set_precision
///////////////////////////////////////////////////////////////////////////////

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

//shape and digitizer
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/GaussDigitizer.h"

//tracking grid curve
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/GridCurve.h"

#include "DGtal/io/boards/Board2D.h"
///////////////////////////////////////////////////////////////////////////////

using namespace DGtal;
using namespace Z2i;

template <typename T, typename F>
T perimeter(const MelkmanConvexHull<Z2i::Point, F>& cvx, const T h) {

  T p = 0.f;
  assert(cvx.size() >= 2 && "convex Hull should contain at least two points");

  for (size_t i = 0; i < cvx.size()-1; i++) {
    p += (cvx[i]-cvx[i+1]).norm();
  }
  p += (cvx[cvx.size()-1]-cvx[0]).norm();

  return p*h;
}

template <typename T, typename F>
T area(const MelkmanConvexHull<Z2i::Point, F>& cvx, const T h) {

  T area = 0.f;
  assert(cvx.size() >= 3 && "convex Hull should contain at three least two points");

  for (size_t i = 0; i < cvx.size()-1; ++i) {
    area += cvx[i][0]*cvx[i+1][1] - cvx[i][1]*cvx[i+1][0];
  }
  area += cvx[cvx.size()][0]*cvx[0][1] - cvx[cvx.size()][1]*cvx[0][0];

  return abs(area*h*h)/2.0f;
}

template <typename T>
std::string toString(const T& number, const int precision) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << number;
  return ss.str();
}

int main() {
  // define an Euclidean shape (disk)
  typedef ImplicitBall<Z2i::Space> Disk;

  Disk disk(Z2i::Point(0,0), 10);

  // Gauss discretization
  double h = 0.1f;
  std::cout << "discretization scale : " << h << std::endl;
  
  GaussDigitizer<Z2i::Space,Disk> dig;
  dig.attach(disk);
  dig.init(disk.getLowerBound()+Z2i::Vector(-1,-1),disk.getUpperBound()+Z2i::Vector(1,1), h);
  
  // Domain domain(dig.getLowerBound(), dig.getUpperBound());
  Domain domain = dig.getDomain();
  
  Z2i::KSpace ks; // make a Kovalevsky-Khalimsky space
  ks.init( dig.getLowerBound(), dig.getUpperBound(), true);
  
  SurfelAdjacency<2> sAdj(true); // set an adjacency (4-connectivity)

  // search for one boundary element
  Z2i::SCell bel = Surfaces<Z2i::KSpace>::findABel(ks, dig, 1000);
  std::vector<Z2i::Point> boundaryPoints; // boundary tracking
  Surfaces<Z2i::KSpace>::track2DBoundaryPoints( boundaryPoints, ks, sAdj, dig, bel );

  Z2i::Curve c;
  c.initFromVector(boundaryPoints);  
  {
    // draw a boundary curve and make a pdf file
    Board2D aBoard;
    aBoard << c; 
    std::string name = "../img/boundaryCurve"+toString(h, 2)+".png";
    aBoard.saveCairo(name.c_str(), Board2D::CairoPNG);
  }

  Z2i::DigitalSet set(domain);

  Shapes<Domain>::digitalShaper(set, dig);

  std::cout << "using digital shaper" << std::endl;

  float perimeterDS = static_cast<double>(boundaryPoints.size()-1)*h;
  std::cout << "perimeter:" << perimeterDS << std::endl;
  
  float airDS = static_cast<double>(set.size()) * h * h;
  std::cout << "area:" << airDS << std::endl;

  // Make ConvexHull
  typedef InHalfPlaneBySimple3x3Matrix<Z2i::Point, DGtal::int64_t> Functor;
  Functor functor;
  MelkmanConvexHull<Z2i::Point, Functor> cvx(functor);
  for (const auto &p : boundaryPoints)
    cvx.add(p);
  
  {
    // scan the CVX points & draw the edges
    Board2D aBoard;
    aBoard.setPenColor(Color::Red);

    auto p = *cvx.begin();
    for(auto it = cvx.begin()+1; it != cvx.end(); ++it) {
      const auto q = *it;
      aBoard.drawArrow(p[0]-0.5, p[1]-0.5, q[0]-0.5, q[1]-0.5);
      p=q;
    }
    aBoard.drawArrow(p[0]-0.5, p[1]-0.5, (*cvx.begin())[0]-0.5, (*cvx.begin())[1]-0.5);

    std::string name = "../img/convexHull"+toString(h, 2)+".png";
    aBoard.saveCairo(name.c_str(), Board2D::CairoPNG);
  }

  std::cout << "using convex hull" << std::endl;
  std::cout << "perimeter:" << perimeter(cvx, h) << std::endl;
  std::cout << "area:" << area(cvx, h) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////