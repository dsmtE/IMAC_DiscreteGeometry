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
int exportCurves(const Domain& domain, const std::vector<T>& objects, const bool is4_8, const std::string& fileName, const std::string& suffix, const bool boundElimination = false) {
    int count = 0;
    Board2D aBoard;
    for (auto &obj : objects) {
        const Z2i::Curve curve = boundaryCurve(obj, is4_8);

        if(!boundElimination || CurveIsInBound(curve, domain)) {

            aBoard << curve;
            /*
            for (const Z2i::SCell &point : curve) {
                auto tmp = point.preCell().coordinates;
                aBoard << tmp;
            }
            */
            ++count;
        }
    }
    const std::string filePath = "./img/boundaryCurve_" + fileName + "_" + (is4_8 ? "4_8" : "8_4")  + "_" + suffix + ".pdf";
    aBoard.saveCairo(filePath.c_str(), Board2D::CairoPDF);
    std::cout << "boundaries saved as path : " << filePath << std::endl;

    return count;
}

template <class T>
void exportCurvesAndDss(const Domain& domain, const std::vector<T>& objects, const bool is4_8, const std::string& fileName) {

    Board2D aBoard;
    for (auto &obj : objects) {
        vector<Z2i::Point> boundaryPoints = getBoundariesPoints(obj, is4_8);
        Z2i::Curve curve;
        curve.initFromVector(boundaryPoints);

        if(CurveIsInBound(curve, domain)) {
            aBoard << curve;

            // Construct the Freeman chain
            Contour4 boundaryBorder(boundaryPoints);
            
            // Segmentation
            Decomposition4 boundaryDecomposition(boundaryBorder.begin(),boundaryBorder.end(), DSS4());

            // Draw each segment
            for (auto it = boundaryDecomposition.begin(), itEnd = boundaryDecomposition.end(); it != itEnd; ++it) {
                aBoard << SetMode("ArithmeticalDSS", "BoundingBox");
                aBoard << CustomStyle("ArithmeticalDSS/BoundingBox", new CustomPenColor(Color::Blue));
                aBoard << it->primitive();
            }
        }
    }
    const std::string filePath = "./img/CurvesAndDss_" + fileName + "_" + (is4_8 ? "4_8" : "8_4") + ".pdf";
    aBoard.saveCairo(filePath.c_str(), Board2D::CairoPDF);
    std::cout << "boundaries saved as path : " << filePath << std::endl;
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
    // (8,4) adjacency
    vector<ObjectType84> objects84;
    back_insert_iterator<vector<ObjectType84>> inserter84(objects84);

    // 4) Use method writeComponents to obtain connected components
    // (4,8) adjacency
    ObjectType48 connectedComponents48(dt4_8, set2d);
    connectedComponents48.writeComponents(inserter48);
    // (8,4) adjacency
    ObjectType84 connectedComponents84(dt8_4, set2d);
    connectedComponents84.writeComponents(inserter84);

    const bool boundariesElimination = true;
    int count48;
    int count84;

    if(boundariesElimination) {
        auto limit = image.domain().upperBound()*2;
        count48 = exportCurves(image.domain(), objects48, true, fileName, "boundariesElimination", true);
        count84 = exportCurves(image.domain(), objects84, false, fileName, "boundariesElimination", true);
    }else {
        count48 = exportCurves(image.domain(), objects48, true, fileName, "NotWellFormed");
        count84 = exportCurves(image.domain(), objects84, false, fileName, "NotWellFormed");
    }

    std::cout << "Nombre de grains de riz " << fileName <<  " 4_8: " << count48 << std::endl;
    std::cout << "Nombre de grains de riz " << fileName <<  " 8_4: " << count84 << std::endl;

    exportCurvesAndDss(image.domain(), objects48, true, fileName);
    exportCurvesAndDss(image.domain(), objects84, false, fileName);

    return 0;
}
