#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/Color.h"

#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>

#include "cmath" 

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef ImageSelector<Domain, unsigned char >::Type ImageType; // type of image
typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSetType; // Digital set type
typedef Object<DT4_8, DigitalSetType> ObjectType48;
typedef Object<DT8_4, DigitalSetType> ObjectType84;

typedef FreemanChain<int> Contour4; 
typedef ArithmeticalDSSComputer<Contour4::ConstIterator,int,4> DSS4;
typedef GreedySegmentation<DSS4> Decomposition4;

template <class T>
vector<Z2i::Point> getBoundariesPoints(T &object, const bool is4_8) {
    KSpace kSpace; // make a Kovalevsky-Khalimsky space
    // we need to add a margine to prevent situations such that an object touch the bourder of the domain
    kSpace.init(object.domain().lowerBound() - Point(1, 1), object.domain().upperBound() + Point(1, 1), true);

    SurfelAdjacency<2> sAdj(is4_8); // set an adjacency (4-connectivity)

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(kSpace, object.pointSet(), 100000);

    // boundary points
    vector<Z2i::Point> boundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(boundaryPoints, kSpace, sAdj, object.pointSet(), bel);

    return boundaryPoints;
}

template <class T>
Curve boundaryCurve(T &object, const bool is4_8) {
    // obtain a curve from points
    Z2i::Curve boundaryCurve;
    boundaryCurve.initFromVector(getBoundariesPoints(object, is4_8));

    return boundaryCurve;
}

bool CurveIsInBound(const Z2i::Curve& curve, const Domain& domain) {

    const int limitX = domain.upperBound()[0] * 2;
    const int limitY = domain.upperBound()[1] * 2;

    bool in = true;
    for (auto &p : curve) {
        PointVector<2, Integer> point = p.preCell().coordinates;
        if (point[0] <= 0 || point[0] >= limitX || point[1] <= 0 || point[1] >= limitY) {
            in = false;
            break;
        }
    }
    return in;
}



template <class T>
void writeCsv(std::string filename, std::vector<std::pair<std::string, std::vector<T>>> dataset) {
    std::ofstream csv(filename);

    struct comma_separator : std::numpunct<char> {
        virtual char do_decimal_point() const override { return ','; }
    };
    const std::string separator(";");
    csv.imbue(std::locale(std::cout.getloc(), new comma_separator));

    for(std::size_t c = 0; c < dataset.size(); c++) {
        csv << dataset[c].first;
        if(c != dataset.size() - 1) { csv << separator; }
    }
    csv << std::endl;

    for(std::size_t i = 0; i < dataset[0].second.size(); i++) {
        for(std::size_t c = 0; c < dataset.size(); c++) {
            csv << dataset[c].second[i];
            if(c != dataset.size() - 1) { csv << separator; }
        }
        csv << std::endl;
    }
    csv.close();
    std::cout << "file " << filename << " saved." << std::endl; 
}

int main(int argc, char** argv) {
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings

    if( argc < 2 ){
        std::cout << "please give fileName" << std::endl;
        return 0;
    }

    const std::string fStart = "./RiceGrains/Rice_";
    const std::string fileName(argv[1]);
    const std::string fEnd = "_seg_bin.pgm";
    std::cout << "filePath: " << fStart << fileName << fEnd << std::endl;

    ImageType image = PGMReader<ImageType>::importPGM(fStart + fileName + fEnd);

    // 1) Create a digital set of proper size
    Z2i::DigitalSet set2d(image.domain());
    // 2) Use SetFromImage::append() to populate a digital set from the input image
    SetFromImage<Z2i::DigitalSet>::append<ImageType>(set2d, image, 1, 255);

    // 3) Create a digital object from the digital set
    // (4,8) adjacency
    vector<ObjectType48> objects48; // All conected components are going to be stored in it
    back_insert_iterator<vector<ObjectType48>> inserter48(objects48); // Iterator used to populated "objects".

    // 4) Use method writeComponents to obtain connected components
    // (4,8) adjacency
    ObjectType48 connectedComponents48(dt4_8, set2d);
    connectedComponents48.writeComponents(inserter48);

    std::vector<double> perimeters;
    std::vector<double> areas;
    std::vector<double> circularities;

    const bool is4_8 = true;
    for (auto &obj : objects48) {

        double perimeter = 0;
        double area = 0;

        vector<Z2i::Point> boundaryPoints = getBoundariesPoints(obj, is4_8);
        Z2i::Curve curve;
        curve.initFromVector(boundaryPoints);

        if(CurveIsInBound(curve, image.domain())) {
            // Construct the Freeman chain
            Contour4 boundaryBorder(boundaryPoints);
            
            // Segmentation
            Decomposition4 boundaryDecomposition(boundaryBorder.begin(),boundaryBorder.end(), DSS4());
            
            // for each segment
            Z2i::Point firstPoint = boundaryDecomposition.begin().get().front();
            Z2i::Point lastPoint;
            for (auto it = boundaryDecomposition.begin(); it != boundaryDecomposition.end(); ++it) {
                const Z2i::Point a = it.get().front();
                const Z2i::Point b = it.get().back();
                perimeter += (b-a).norm();
                area += a[0] * b[1] - a[1] * b[0];
                lastPoint = b;
            }

            perimeter += (firstPoint-lastPoint).norm();
            area += lastPoint[0] * firstPoint[1] - lastPoint[1] * firstPoint[0];
            area = abs(area) / 2.0;

            const double circularity = (4 * M_PI * area) / (perimeter * perimeter);

            perimeters.push_back(perimeter);
            areas.push_back(area);
            circularities.push_back(circularity);
        }
    }

    // export to csv
    
    std::vector<std::pair<std::string, std::vector<double>>> dataset = {{"perimeter", perimeters}, {"area", areas}, {"circularity", circularities}};
    writeCsv("./data_"+ fileName +".csv", dataset);

    /*
    ofstream data("./data_"+ fileName +".csv");
    
    data << "perimeter, ";
    // for(const double x : perimeters) { data << x << ", "; }
    std::copy(perimeters.begin(), perimeters.end(), std::ostream_iterator<double>(data, ", "));
    data << std::endl;

    data << "area, ";
    std::copy(areas.begin(), areas.end(), std::ostream_iterator<double>(data, ", "));
    data << std::endl;

    data << "circularity, ";
    std::copy(circularities.begin(), circularities.end(), std::ostream_iterator<double>(data, ", "));
    data << std::endl;

    data.close();
    */

    return 0;
}
