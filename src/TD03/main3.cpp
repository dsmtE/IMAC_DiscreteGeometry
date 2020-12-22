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

#include  <numeric>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "DGtal/images/RigidTransformation2D.h"
#include "DGtal/images/ConstImageAdapter.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;
using namespace functors;

typedef ImageSelector<Domain, unsigned char >::Type ImageType; // type of image
typedef DigitalSetSelector<Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSetType; // Digital set type
typedef Object<DT4_8, DigitalSetType> ObjectType48;
typedef Object<DT8_4, DigitalSetType> ObjectType84;

// Rigid transformation
typedef ForwardRigidTransformation2D < Space > ForwardTrans;
typedef BackwardRigidTransformation2D < Space > BackwardTrans;
typedef ConstImageAdapter<ImageType, Domain, BackwardTrans, ImageType::Value, Identity > MyImageBackwardAdapter;
typedef DomainRigidTransformation2D < Domain, ForwardTrans > MyDomainTransformer;
typedef MyDomainTransformer::Bounds Bounds;


const ObjectType48& keepLargerComponent(const vector<ObjectType48>& objs) {
    size_t maxId = 0;
    for(size_t i = 0; i < objs.size(); i++) {
        if(objs[i].size() > objs[maxId].size()) {
            maxId = i;
        }
    }
    return objs[maxId];
}

// source exemple : https://projet.liris.cnrs.fr/dgtal/doc/0.9.1/exampleRigidtransformation2d_8cpp_source.html
void rigidTransformAndExport(const ImageType& img, const Z2i::Point& origin , const double& angle, const Z2i::Point& translation, const std::string& fileName) {
    // origin, angle in radians, and translation
    ForwardTrans forwardTrans(origin, angle, translation);
    BackwardTrans backwardTrans(origin, angle, translation);
    MyDomainTransformer domainTransformer(forwardTrans);
    Identity idD;

    // Bounds bounds = domainTransformer(img.domain());
    // Domain transformedDomain(bounds.first, bounds.second);
    // MyImageBackwardAdapter backwardImageAdapter(img, transformedDomain, backwardTrans, idD);

    // keep the same domain for comparaison
    MyImageBackwardAdapter backwardImageAdapter(img, img.domain(), backwardTrans, idD);
    backwardImageAdapter >> (fileName + ".pgm");
}

std::vector<Z2i::Point> getPts(const ObjectType48& obj) {
    std::vector<Z2i::Point> vec;

    for(const Z2i::Point& p : obj) {
        vec.push_back(p);
    }
    return vec;
}

Eigen::MatrixXd toEigenMat(const std::vector<Z2i::Point>& pts) {
    Eigen::MatrixXd mat(pts.size(), 2);
    for(std::size_t i = 0; i < pts.size(); i++) {
        mat.row(i) << pts[i][0], pts[i][1];
    }

    return mat;
}

Z2i::Point getMassCenter(const ObjectType48& obj) {
    const std::vector<Z2i::Point> pts = getPts(obj);
    static const Z2i::Point zero(0.0, 0.0);
    Z2i::Point c = std::accumulate(pts.begin(), pts.end(), zero);

    c[0] /= pts.size();
    c[1] /= pts.size();
    return c;
}

Eigen::Vector2d getPrincipalEigenValue(const ObjectType48& obj) {
    Eigen::MatrixXd pts = toEigenMat(getPts(obj));

    // Mean centering pts
    Eigen::MatrixXd centered = pts.rowwise() - pts.colwise().mean();
    
    // Compute the covariance matrix
    // cov(X,Y) = E[(X-E[X])(Y-E[Y])]
    // for centered values it's just cov(X,Y) = E[XY]
    // here X & Y are vectors
    Eigen::MatrixXd cov = centered.transpose() * centered;
    cov = cov / (centered.rows() - 1);

    static const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    
    // std::cout << "cov:" << std::endl << cov << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    // std::cout << "eigenvectors:" << std::endl << eig.eigenvectors().format(CleanFmt) << std::endl;

    // return last eigenVec (max eigenVal)
    return eig.eigenvectors().col(1);
}

int main(int argc, char** argv) {
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings

    const std::string fStart = "./binary/IMG_";
    const std::string fEnd = ".pgm";

    std::vector<Z2i::Point> massCenters;
    std::vector<Eigen::Vector2d> maxEigenVecs;

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
        maxEigenVecs.push_back(getPrincipalEigenValue(larger));
    }

    std::vector<Z2i::Point> translations;
    std::cout << "sucessives translations : ";
    for(size_t i = 0; i < massCenters.size()-1; ++i) {
        const Z2i::Point& c1 = massCenters[i];
        const Z2i::Point& c2 = massCenters[i+1];
        translations.push_back(c2-c1);
        std::cout << translations.back() << ",";
    }
    std::cout << std::endl;

    std::vector<double> rotations;
    std::cout << "sucessives Rotation Angle : ";
    for(size_t i = 0; i < maxEigenVecs.size()-1; ++i) {
        const Eigen::Vector2d& v1 = maxEigenVecs[i];
        const Eigen::Vector2d& v2 = maxEigenVecs[i+1];
        rotations.push_back(atan2(v2[1], v2[0]) - atan2(v1[1], v1[0]));
        std::cout << rotations.back() << ",";
    }
    std::cout << std::endl;

    static const Z2i::Point zero(0.0, 0.0);

    ImageType image = PGMReader<ImageType>::importPGM(fStart + "2763" + fEnd);
    double currentAngle = 0.0;
    Z2i::Point currentTranslation(0.0, 0.0);
    for(int i = 0; i < rotations.size(); ++i) {
        ImageType image = PGMReader<ImageType>::importPGM(fStart + std::to_string(2764 + i) + fEnd);

        currentAngle -= rotations[i];
        currentTranslation -= translations[i];
        std::cout << "rot: " << currentAngle << " trans: " << currentTranslation << std::endl;
        // rotate around the center of mass of the corresponding image (i+1)
        rigidTransformAndExport(image, massCenters[i+1], currentAngle, currentTranslation, "./out/IMG_" + std::to_string(2764 + i) + "_motionTo2763");
    }

    return 0;
}
