#!/bin/bash

mkdir -p Animate

for i in `find P* -name "*pathview.png"`
	do      
		#echo "path: $i"
		base="$(basename "$i" .png)"
		#echo "basename: $base"
		dir="$(dirname "$i")"
		#echo "Directory: $dir"
		new=$base\_$dir.png
		#echo "New name: $new"
		new_dir="${base%%.*}"
		#echo $new_dir
		mkdir -p Animate/$new_dir
		cp $i Animate/$new_dir/$new	
	done

cd Animate

for i in ssc*/
	do
		#echo $i
		input_sig=$(ls $i/*sig.png)		
		#echo $input_sig
		convert -delay 20 -loop 0 -morph 5 $input_sig $i/Timeline_sig.gif 
		input_nonsig=$(ls $i/*P?.png)
		convert -delay 20 -loop 0 -morph 5 $input_nonsig $i/Timeline_all.gif
	done
