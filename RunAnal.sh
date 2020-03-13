#!/bin/bash

# Run as: ./RunAnal.sh Part1

dir="Outputs/"$1

if [ -z "$dir" ]; then
	echo "You need to declare target directory."
	echo "Example: ./RunAnal.sh Part1"
	exit
fi

make -C build/

if [ $? -eq 0 ] ; then
  echo "Make was successfully done..."
else
  echo "Make has failed..."
  exit
fi



if [ -d $dir ]; then
	echo "Target file exists"
	echo "Results will be updated"
	if [ ! -d $dir/Figs ]; then
		mkdir $dir/trackQuality
		echo "Creating folder: " $dir/trackQuality
	fi
else
	mkdir $dir
	echo "Creating folder: " $dir
	mkdir $dir/Figs
	echo "Creating folder: " $dir/Figs
fi

build/Analysis $1

echo "Done! Analysis has been finished"