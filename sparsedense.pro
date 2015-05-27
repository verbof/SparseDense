TEMPLATE	= app
LANGUAGE	= C++

#LIBS	        = -L/users/drosos/pardiso -lmetis41-P_pardiso -lmetis41_pardiso -lpardiso -lpils_pardiso -lmetis41-P_pardiso -lmetis41_pardiso -lpardiso -lpils_pardiso

LIBS	        = -L${HOME}/Libraries/linuxAMD64 -lpardiso500-GNU481-X86-64
LIBS           += -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64             \
                  -lmkl_lapack95_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core \ 
                  -lgfortran -lgomp /opt/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64/libirc.a

CCX              = /opt/cray/craype/2.2.1/bin/CC
#CCX              = /opt/cray/craype/2.2.1/bin/CC
QMAKE_CXX       = $$CCX
QMAKE_CC        = $$CCX
QMAKE_LINK      = $$CCX
QMAKE_LFLAGS    = -dynamic

DEFINES         = 

INCLUDEPATH	+= src 


CONFIG-=qt
CONFIG+=release
SOURCES+=src/CSRdouble.cpp   \
	 src/CSRcomplex.cpp  \
	 src/RealMath.cpp    \
	 src/ParDiSO.cpp     \
	 src/tools.cpp       \
	 src/schur.cpp       \
	 src/IO.cpp          \
	 src/readinput.cpp   \
	 src/readdist.cpp    \
	 src/main.cpp        \

DESTDIR     = bin
OBJECTS_DIR = obj


unix {
  OBJECTS_DIR = obj
}



