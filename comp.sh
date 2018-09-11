#LHAPDFinc="/usr/local/Cellar/lhapdf/6.1.5/include"
#LHAPDFlib="/usr/local/Cellar/lhapdf/6.1.5/lib -lLHAPDF"
LHAPDFinc="/usr/local/Cellar/lhapdf/6.2.1/include"
LHAPDFlib="/usr/local/Cellar/lhapdf/6.2.1/lib -lLHAPDF"

ROOTlib=$(root-config --libs) 
ROOTinc=$(root-config --incdir) 

target="evgen"
g++ -std=c++11 -fPIC -o $target -I/Users/kazuki/Packages/herwig7/include -I$LHAPDFinc -I$ROOTinc -L/Users/kazuki/Packages/herwig7/lib/ThePEG -L$LHAPDFlib $ROOTlib -lThePEG $target.cc

target="muon_decay"
g++ -std=c++11 -fPIC -o $target -I/Users/kazuki/Packages/herwig7/include -I$LHAPDFinc -I$ROOTinc -L/Users/kazuki/Packages/herwig7/lib/ThePEG -L$LHAPDFlib $ROOTlib -lThePEG $target.cc
#g++ -std=c++11 -fPIC -o xsec -I/Users/kazuki/Packages/herwig7/include -I$LHAPDFinc -I$ROOTinc -L/Users/kazuki/Packages/herwig7/lib/ThePEG -L$LHAPDFlib $ROOTlib -lThePEG full.cc
#g++ -std=c++11 -fPIC -o rsgen -I$LHAPDFinc -I$ROOTinc -L$LHAPDFlib $ROOTlib rsgen.cc
