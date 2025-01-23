LDFLAGS= -L${shell root-config --libdir}  -L${ARA_DEPS_INSTALL_DIR}/lib -L${ARA_UTIL_INSTALL_DIR}/lib -L${ARA_UTIL_INSTALL_DIR}/../../RHEL_7_x86_64/source/AraSim 
CXXFLAGS= -I${shell root-config --cflags} -I${ARA_DEPS_INSTALL_DIR}/include -I${ARA_UTIL_INSTALL_DIR}/include -I${ARA_UTIL_INSTALL_DIR}/../../RHEL_7_x86_64/source/AraSim -I../tools
LDLIBS += ${shell root-config --libs} -lAraEvent -lRootFftwWrapper -lAra -lAraCorrelator -lAraVertex  -lMathMore

