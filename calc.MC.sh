
for f in *depth.stats.txt
do 
	min=$(cat $f | awk '{print $3}' | grep -v 'MEDIAN' | sort -h | tail -n 1)
	max=$(cat $f | awk '{print $3}' | grep -v 'MEDIAN' | sort -h | head -n 1)
	echo $f $min $max
done
