ANA_HOME=${HOME}/cernbox/AnaSourceCode
FIT_INCLUDE=${ANA_HOME}/SelfHeadFiles/FitFunction
FIT_LIB=${ANA_HOME}/SelfLibrary/FitFunction
MATH_TECHNICAL_INCLUDE=${ANA_HOME}/SelfHeadFiles/MathTechnicalFunc
MATH_TECHNICAL_LIB=${ANA_HOME}/SelfLibrary/MathTechnicalFunc
PTCENT=${ANA_HOME}/SelfHeadFiles/PtCentRange
#SRC=backup_AcceptanceEffi_PbPb_for_tuning.cpp
SRC=AcceptanceEffi_PbPb_for_tuning.cpp

OBJS = $(SRC:.c=.o)

GOrun: $(OBJS)
	g++ -o GOrun $(SRC) `root-config --cflags --libs` -I${ALICE_ROOT}/include -L${ALICE_ROOT}/lib -lSTEERBase -I${PTCENT} -I${FIT_INCLUDE} -L${FIT_LIB} -lFitFunction -I${MATH_TECHNICAL_INCLUDE} -L${MATH_TECHNICAL_LIB} -lMathTechnicalFunc

clean:
	rm GOrun
