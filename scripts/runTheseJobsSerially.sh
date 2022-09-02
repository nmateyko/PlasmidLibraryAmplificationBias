#!/bin/bash
while read p; do
	$1 $p
done < $2
