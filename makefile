#Paramétrage du compilateur et des options de compilations des modules
COMP	= gcc
OPTIONS	= -c -ansi -Wall


#Paramétrage de l'exécutable et de la liste des objets qui le composent
EXECUTABLE	= mainProgram
OBJETS		= main.o kineticSolver.o solverRusanov.o testCase.o initialConditions.o usefull.o timeSetting.o topography.o boundaryCondition.o -lm

#Description des cibles
${EXECUTABLE}: ${OBJETS}
	${COMP} -o ${EXECUTABLE} ${OBJETS}

main.o: main.c kineticSolver.h initialConditions.h usefull.h timeSetting.h topography.h testCase.h boundaryCondition.h constant.h
	${COMP} ${OPTIONS} main.c

testCase.o: testCase.c testCase.h usefull.h constant.h
	${COMP} ${OPTIONS} testCase.c
	
kineticSolver.o: kineticSolver.c kineticSolver.h usefull.h constant.h
	${COMP} ${OPTIONS} kineticSolver.c
	
solverRusanov.o: solverRusanov.c solverRusanov.h usefull.h constant.h
	${COMP} ${OPTIONS} solverRusanov.c

initialConditions.o: initialConditions.c initialConditions.h constant.h
	${COMP} ${OPTIONS} initialConditions.c
	
boundaryCondition.o: boundaryCondition.c boundaryCondition.h constant.h
	${COMP} ${OPTIONS} boundaryCondition.c

usefull.o: usefull.c usefull.h constant.h
	${COMP} ${OPTIONS} usefull.c
	
topography.o: topography.c topography.h constant.h
	${COMP} ${OPTIONS} topography.c
	
timeSetting.o: timeSetting.c timeSetting.h usefull.h constant.h
	${COMP} ${OPTIONS} timeSetting.c
