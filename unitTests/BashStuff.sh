#!/bin/bash

for i in {0..100}
do 
	number=$(printf %.6d $i)
	fileName=output@-$number.vtk
	#echo $i
	mv $fileName EulerOutputSample/
#	rm *@-$i.vtk
done

#for fileToMove in *\).vtk
#do 
#	stringLength=${#fileToMove}

#	#Find the new file name.
#	i=$stringLength
#	i=$[$i-1]
#	while [ "${fileToMove:$i:1}" != "(" ]
#	do
#		i=$[$i-1]
#	done
#	
#	i=$[$i-1]
#	newStringLength=$i

#	newFile=${fileToMove:0:$newStringLength}.vtk

#	mv "$fileToMove" "$newFile"
#done
