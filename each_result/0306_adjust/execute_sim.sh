#!/bin/bash

function setStarttime() {
        start_time=`date +%s`
}

function getEndtime() {
        end_time=`date +%s`

        SS=`expr ${end_time} - ${start_time}`

        HH=`expr ${SS} / 3600`
        SS=`expr ${SS} % 3600`
        MM=`expr ${SS} / 60`
        SS=`expr ${SS} % 60`

        echo "${HH}:${MM}:${SS}"
}

read -p "Count(default=1): " input
expr $input + 1 > /dev/null 2>&1
RET=$?
if [ $RET -lt 2 ]; then
    count=$input
else
    count=1
fi

read -p "Wait (y/N): " yn
if [[ $yn = 'y' ]] ; then
  isptime="-ptime"
else
  isptime=""
fi

num=`ls -U1 ./dat/*.dat | wc -l`
file_num=`expr $num`
file_index=0

mkdir -p result

touch result/time.txt
echo '' > ./result/time.txt

for file in ./dat/*.dat
do
	echo "java -jar my_sim.jar -pfile: $file -batchRun $isptime"
	file_index=$(( file_index + 1 ))
  SECONDS=0
	for i in `seq 1 $count`
	do
    setStarttime
		echo "Count: $i / $count --- $file_index of $file_num"
		java -jar src/my_sim.jar -pfile: $file -batchRun $isptime
    getEndtime
	done
  echo "$file: $SECONDS / $count"
  mv *.txt ./result
done
