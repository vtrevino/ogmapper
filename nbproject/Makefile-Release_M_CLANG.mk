#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=clang
CCC=clang++
CXX=clang++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=CLang-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release_M_CLANG
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/BGZF.o \
	${OBJECTDIR}/CircularArray.o \
	${OBJECTDIR}/SimpleStatisticsCalculator.o \
	${OBJECTDIR}/binnedBigCounters.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/ogAligner.o \
	${OBJECTDIR}/ogBitwiseAT1GC0Encoding.o \
	${OBJECTDIR}/ogCandidatePosManager.o \
	${OBJECTDIR}/ogCigarOperations.o \
	${OBJECTDIR}/ogCountingFunctions.o \
	${OBJECTDIR}/ogFastAQGZreader.o \
	${OBJECTDIR}/ogGappedBitwiseAT1GC0Encoding.o \
	${OBJECTDIR}/ogGenome.o \
	${OBJECTDIR}/ogGenomePositions.o \
	${OBJECTDIR}/ogGuider.o \
	${OBJECTDIR}/ogHomoPolymerCompressedEnc.o \
	${OBJECTDIR}/ogIndex-original.o \
	${OBJECTDIR}/ogIndex.o \
	${OBJECTDIR}/ogKeyEncoding.o \
	${OBJECTDIR}/ogKeys.o \
	${OBJECTDIR}/ogMappingFunctions.o \
	${OBJECTDIR}/ogPlainEncoding.o \
	${OBJECTDIR}/ogReadKeyMapping.o \
	${OBJECTDIR}/ogReadsMapper.o \
	${OBJECTDIR}/ogReusableBinaryTree.o \
	${OBJECTDIR}/ogSamWriter.o \
	${OBJECTDIR}/ogSingleRead.o \
	${OBJECTDIR}/ogStateMachineGuider.o \
	${OBJECTDIR}/ogSwapBitwiseAT1GC0Encoding.o \
	${OBJECTDIR}/ogTupleGuider.o \
	${OBJECTDIR}/sequenceRead.o


# C Compiler Flags
CFLAGS=-std=c++11 -target arm64-apple-darwin

# CC Compiler Flags
CCFLAGS=-std=c++11 -target arm64-apple-darwin
CXXFLAGS=-std=c++11 -target arm64-apple-darwin

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=lib/WFA2/MAC-M/libwfacpp.a lib/LIBOMP/MAC-M/libomp.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ogmapper

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ogmapper: lib/WFA2/MAC-M/libwfacpp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ogmapper: lib/LIBOMP/MAC-M/libomp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ogmapper: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ogmapper ${OBJECTFILES} ${LDLIBSOPTIONS} -lz -lomp

${OBJECTDIR}/BGZF.o: BGZF.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/BGZF.o BGZF.cpp

${OBJECTDIR}/CircularArray.o: CircularArray.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CircularArray.o CircularArray.cpp

${OBJECTDIR}/SimpleStatisticsCalculator.o: SimpleStatisticsCalculator.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SimpleStatisticsCalculator.o SimpleStatisticsCalculator.cpp

${OBJECTDIR}/binnedBigCounters.o: binnedBigCounters.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/binnedBigCounters.o binnedBigCounters.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/ogAligner.o: ogAligner.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogAligner.o ogAligner.cpp

${OBJECTDIR}/ogBitwiseAT1GC0Encoding.o: ogBitwiseAT1GC0Encoding.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogBitwiseAT1GC0Encoding.o ogBitwiseAT1GC0Encoding.cpp

${OBJECTDIR}/ogCandidatePosManager.o: ogCandidatePosManager.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogCandidatePosManager.o ogCandidatePosManager.cpp

${OBJECTDIR}/ogCigarOperations.o: ogCigarOperations.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogCigarOperations.o ogCigarOperations.cpp

${OBJECTDIR}/ogCountingFunctions.o: ogCountingFunctions.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogCountingFunctions.o ogCountingFunctions.cpp

${OBJECTDIR}/ogFastAQGZreader.o: ogFastAQGZreader.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogFastAQGZreader.o ogFastAQGZreader.cpp

${OBJECTDIR}/ogGappedBitwiseAT1GC0Encoding.o: ogGappedBitwiseAT1GC0Encoding.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogGappedBitwiseAT1GC0Encoding.o ogGappedBitwiseAT1GC0Encoding.cpp

${OBJECTDIR}/ogGenome.o: ogGenome.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogGenome.o ogGenome.cpp

${OBJECTDIR}/ogGenomePositions.o: ogGenomePositions.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogGenomePositions.o ogGenomePositions.cpp

${OBJECTDIR}/ogGuider.o: ogGuider.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogGuider.o ogGuider.cpp

${OBJECTDIR}/ogHomoPolymerCompressedEnc.o: ogHomoPolymerCompressedEnc.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogHomoPolymerCompressedEnc.o ogHomoPolymerCompressedEnc.cpp

${OBJECTDIR}/ogIndex-original.o: ogIndex-original.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogIndex-original.o ogIndex-original.cpp

${OBJECTDIR}/ogIndex.o: ogIndex.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogIndex.o ogIndex.cpp

${OBJECTDIR}/ogKeyEncoding.o: ogKeyEncoding.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogKeyEncoding.o ogKeyEncoding.cpp

${OBJECTDIR}/ogKeys.o: ogKeys.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogKeys.o ogKeys.cpp

${OBJECTDIR}/ogMappingFunctions.o: ogMappingFunctions.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogMappingFunctions.o ogMappingFunctions.cpp

${OBJECTDIR}/ogPlainEncoding.o: ogPlainEncoding.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogPlainEncoding.o ogPlainEncoding.cpp

${OBJECTDIR}/ogReadKeyMapping.o: ogReadKeyMapping.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogReadKeyMapping.o ogReadKeyMapping.cpp

${OBJECTDIR}/ogReadsMapper.o: ogReadsMapper.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogReadsMapper.o ogReadsMapper.cpp

${OBJECTDIR}/ogReusableBinaryTree.o: ogReusableBinaryTree.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogReusableBinaryTree.o ogReusableBinaryTree.cpp

${OBJECTDIR}/ogSamWriter.o: ogSamWriter.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogSamWriter.o ogSamWriter.cpp

${OBJECTDIR}/ogSingleRead.o: ogSingleRead.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogSingleRead.o ogSingleRead.cpp

${OBJECTDIR}/ogStateMachineGuider.o: ogStateMachineGuider.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogStateMachineGuider.o ogStateMachineGuider.cpp

${OBJECTDIR}/ogSwapBitwiseAT1GC0Encoding.o: ogSwapBitwiseAT1GC0Encoding.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogSwapBitwiseAT1GC0Encoding.o ogSwapBitwiseAT1GC0Encoding.cpp

${OBJECTDIR}/ogTupleGuider.o: ogTupleGuider.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ogTupleGuider.o ogTupleGuider.cpp

${OBJECTDIR}/sequenceRead.o: sequenceRead.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sequenceRead.o sequenceRead.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
