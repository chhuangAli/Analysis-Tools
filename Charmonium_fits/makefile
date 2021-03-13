ANA_HOME=${HOME}/cernbox/AnaSourceCode
FIT_INCLUDE=${ANA_HOME}/SelfHeadFiles/FitFunction
FIT_LIB=${ANA_HOME}/SelfLibrary/FitFunction
READCOUNT_INCLUDE=${ANA_HOME}/SelfHeadFiles/MathTechnicalFunc
READCOUNT_LIB=${ANA_HOME}/SelfLibrary/MathTechnicalFunc
PTCENT=${ANA_HOME}/SelfHeadFiles/PtCentRange
SRC=Fit_jpsi_psi2S.cpp
FUNC=read_access_histos.cpp
FUNC2=read_access_subtracted_histos.cpp
FUNC3=read_access_histos_y_dependence.cpp

OBJS = $(SRC:.c=.o)
OBJS += $(FUNC:.c=.o)
OBJS += $(FUNC2:.c=.o)
OBJS += $(FUNC3:.c=.o)

GOrun: $(OBJS)
	g++ -c -o read_access_histos.o $(FUNC) -I${PTCENT} `root-config --cflags --libs`
	g++ -c -o read_access_subtracted_histos.o $(FUNC2) `root-config --cflags --libs`        
	g++ -c -o read_access_histos_y_dependence.o $(FUNC3) `root-config --cflags --libs`
	g++ -o GOrun $(OBJS) -I${READCOUNT_INCLUDE} -L${READCOUNT_LIB} -lMathTechnicalFunc `root-config --cflags --libs` -I${ALICE_ROOT}/include -L${ALICE_ROOT}/lib -lSTEERBase -I${PTCENT} -I${FIT_INCLUDE} -L${FIT_LIB} -lFitFunction



clean:
	rm GOrun *.o 
