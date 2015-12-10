#!/bin/bash
#****************************************************************************
#****************************************************************************
#
# File: makefile
# Project: unitTests
# Descrition: This script builds the source code for unitTests.  The source
#	is located in WORKINGDIR.  The shared source code and object files
#	are stored in SRC.  Old object files are replaced if their source file
#	is newer.
#
# Special Thanks to Tim Mclaren at Virginia Tech for helping develop
#	this makefile.
# Date: October 7, 2015
#
#****************************************************************************
#****************************************************************************

#Define Directories**********************************************************
WORKINGDIR=$(pwd)
INCLUDE="../include/"
SRC="../src/"
projectName="euler2D"

#Compile Shared Source Code**************************************************
cd $SRC

# Compile all .cpp files
for sourceFile in *.cpp
do
	objectFile="${sourceFile%%.*}.o"
	headerNew=false

	for headerFile in $INCLUDE*.h
	do
		if [[ -e $objectFile && $headerFile -nt $objectFile ]]; then
			headerNew=true
		fi
	done

	if [[ $objectFile -nt $sourceFile && "$headerNew" = false ]]; then
		echo "$sourceFile: Compiled."
	else
		echo "$sourceFile: Compiling..."
		mpic++ -Wall -std=c++0x -fPIC -I$INCLUDE -O2 -c -fopenmp $sourceFile
	fi

done

#Create the shared object file***********************************************
mpic++ -Wall -std=c++0x -shared -O2 -fopenmp -o libmeshBlock.so *.o

#Copy shared object file into WORKINGDIR*************************************
cp libmeshBlock.so $WORKINGDIR

#Compile Project source files************************************************
cd $WORKINGDIR

for sourceFile in *.cpp
do
	objectFile="${sourceFile%%.*}.o"
	headerNew=false

	for headerFile in *.h
	do
		if [[ -e $objectFile && $headerFile -nt $objectFile ]]; then
			headerNew=true
		fi
	done

	for headerFile in $INCLUDE*.h
	do
		if [[ -e $objectFile && $headerFile -nt $objectFile ]]; then
			headerNew=true
		fi
	done

	if [[ $objectFile -nt $sourceFile && "$headerNew" = false ]]; then
		echo "$sourceFile: Compiled."
	else
		echo "$sourceFile: Compiling..."
		mpic++ -Wall -std=c++0x -I$INCLUDE -O2 -c -fopenmp $sourceFile
	fi

done

#Create Executable***********************************************************
mpic++ -Wall -std=c++0x -I$INCLUDE -O2 -fopenmp -o $projectName *.o -L$WORKINGDIR -lmeshBlock
