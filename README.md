# IMAC_DiscreteGeometry
Repository that contains the TDs for our IMAC3 "Discrete Geometry" class

## How to install and build the project :
- install vcpkg by cloning this repository were you want: git clone https://github.com/microsoft/vcpkg and start install using ".\vcpkg\bootstrap-vcpkg.bat"
- install DGtal dependencies through vcpkg : vcpkg install zlib boost --triplet x64-windows
> informations about vcpkg can be found here : https://github.com/microsoft/vcpkg
- change your path to vcpkg in the cmakeInstall.bat file
- Launch the bat file to start the installation
- open the sln (visual studio solution) file in the build build directory
- build our target in visual studio

tester only using MSVC with visual studio
