mkdir -p GC_output

find . -name "per_sequence_gc_content.png" -exec sh -c '
	for file do
		parent=$(dirname "$file")
		echo $parent
		grandparent=$(dirname "$parent")
		echo $grandparent
		echo $file
		dest_name=${grandparent#*/Raw/*}
		echo $dest_name
		cp $file GC_output/$dest_name\_per_sequence_gc_content.png
	done' sh {} +
