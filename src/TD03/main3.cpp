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

const ObjectType48& keepLargerComponent(const vector<ObjectType48>& objs) {
    size_t maxId = 0;
    for(size_t i = 0; i < objs.size(); i++) {
        if(objs[i].size() > objs[maxId].size()) {
            maxId = i;
        }
    }
    return objs[maxId];
}

Z2i::Point getMassCenter(const ObjectType48& obj) {
    Z2i::Point c;
    for(const Z2i::Point& p : obj) {
        c += p;
    }
    c[0] /= (double)obj.size();
    c[1] /= (double)obj.size();

    return c;
}

int main(int argc, char** argv) {
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings

    const std::string fStart = "./binary/IMG_";
    const std::string fEnd = ".pgm";

    std::vector<Z2i::Point> massCenters;

    for(int i = 2763; i < 2767; ++i) {
        
        const std::string filePath = fStart + std::to_string(i) + fEnd;
        std::cout << "filePath: " << filePath << std::endl;

        ImageType image = PGMReader<ImageType>::importPGM(filePath);

        Z2i::DigitalSet set2d(image.domain());
        SetFromImage<Z2i::DigitalSet>::append<ImageType>(set2d, image, 1, 255);

        // Create a digital object from the digital set(4,8) adjacency
        vector<ObjectType48> objects48;
        back_insert_iterator<vector<ObjectType48>> inserter48(objects48);

        ObjectType48 connectedComponents48(dt4_8, set2d);
        connectedComponents48.writeComponents(inserter48);

        const ObjectType48& larger = keepLargerComponent(objects48); // keep only larger component

        massCenters.push_back(getMassCenter(larger));
    }

    /*
    for(const Z2i::Point& p: massCenters) {
        std::cout << p << ", ";
    }
    std::cout << std::endl;
    */

    std::cout << "sucessives translations : ";
    for(size_t i = 0; i < massCenters.size()-1; ++i) {
        const Z2i::Point& c1 = massCenters[i];
        const Z2i::Point& c2 = massCenters[i+1];
        std::cout << c2-c1 << ", ";
    }
    std::cout << std::endl;
    
    return 0;
}
