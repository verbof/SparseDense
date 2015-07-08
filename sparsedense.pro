TEMPLATE	= app
LANGUAGE	= C++

LIBS	        = -L${HOME}/Libraries/linuxAMD64 -lpardiso500-GNU481-X86-64
LIBS           += -L${MKLROOT}/lib/intel64/                         \
                  -lmkl_lapack95_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core                     \ 
                  -lgfortran -lgomp ${INTEL_PATH}.${INTEL_MINOR_VERSION}/compiler/lib/intel64/libirc.a

CCX             = /opt/cray/craype/2.2.1/bin/CC
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



