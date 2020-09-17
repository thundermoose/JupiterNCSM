#!/bin/bash


for programs in $(find ./test_no_logging -regex ".*\\.x")
do
	echo "Running testsin $programs"
	$programs --run-all-tests
	if [ $? -ne 0 ]
	then
		exit 1
	fi
done
exit 0
