./regen.sh
./configure --prefix=$PWD/local
export PKG_CONFIG_PATH=$PWD/local/lib/pkgconfig:$PKG_CONFIG_PATH
