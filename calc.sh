
for f in *chr_depths*
do 
	x=$(awk -F ' ' '{sum+=$2;}END{print sum;}' $f)
	echo $f $x
done
