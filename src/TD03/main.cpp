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

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef ImageSelector<Domain, unsigned char >::Type ImageType; // type of image
typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSetType; // Digital set type
typedef Object<DT4_8, DigitalSetType> ObjectType48;
typedef Object<DT8_4, DigitalSetType> ObjectType84;

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

template <class T>
void exportCurve(const Domain& domain, const T& obj, const bool is4_8, const std::string& fileName) {
    int count = 0;
    Board2D aBoard;

    const Z2i::Curve curve = boundaryCurve(obj, is4_8);
    aBoard << curve;

    const std::string filePath = "./img/" + fileName + "_" + (is4_8 ? "4_8" : "8_4")  + ".pdf";
    aBoard.saveCairo(filePath.c_str(), Board2D::CairoPDF);
    std::cout << "boundaries saved as path : " << filePath << std::endl;
}

const ObjectType48& keepLargerComponent(const vector<ObjectType48>& objs) {
    size_t maxId = 0;
    for(size_t i = 0; i < objs.size(); i++) {
        if(objs[i].size() > objs[maxId].size()) {
            maxId = i;
        }
    }
    return objs[maxId];
}

int main(int argc, char** argv) {
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings

    if( argc < 2 ) {
        std::cout << "please give fileName" << std::endl;
        return 0;
    }

    const std::string fStart = "./binary/IMG_";
    const std::string fileName(argv[1]);
    const std::string fEnd = ".pgm";
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

    // keep only larger component
    ObjectType48 max = keepLargerComponent(objects48);

    exportCurve(image.domain(), max, true, fileName);
    return 0;
}
