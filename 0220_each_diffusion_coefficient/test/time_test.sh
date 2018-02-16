#!/bin/bash

count=10

isptime="-ptime"

mkdir -p result

touch result/time.txt
echo '' > ./result/time.txt

for file in ./dat/*.dat
do
	echo "java -jar my_sim.jar -pfile: $file -batchRun $isptime"
	SECONDS=0
	for i in `seq 1 $count`
	do
    java -jar src/my_sim.jar -pfile: $file -batchRun $isptime
    echo $SECONDS
  done
  echo "$file,$SECONDS,$count" >> ./result/time.txt
	mv *.txt ./result
done

python time_test.py
