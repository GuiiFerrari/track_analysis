#!/bin/sh

cd build
rm -r *
cmake ../
make
mv pATTPCexe ../
#mv libTeBATlib.so ../../ext/lib
#cd ../include
#cp *.h ../../ext/include
cd ../
