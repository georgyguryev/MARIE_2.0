OLDDIR="$(pwd)"
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
cd "$SCRIPTPATH"
mkdir -p ./Vc/build-1.4
cd ./Vc/build-1.4
cmake -DCMAKE_INSTALL_PREFIX=../../../../ -DBUILD_TESTING=OFF ../
make -j$((($(nproc)/2)))
make install
cd "$OLDDIR"
