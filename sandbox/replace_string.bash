#!/bin/bash

if [ "$#" != "2" ]; then
    echo
    echo "Usage: $0 string1 string2"
    echo "   Replaces string 1 with string 2 in most"
    echo "   text file types in subdirectories."
    echo
    exit 1
fi

find . -name "*.cpp" | \
    while read file; do
	sed -e "s/${1}/${2}/g" ${file} > $$.tmp
	cat $$.tmp > ${file}
    done

find . -name "*.h" | \
    while read file; do
	sed -e "s/${1}/${2}/g" ${file} > $$.tmp
	cat $$.tmp > ${file}
    done

find . -name "*.r" | \
    while read file; do
	sed -e "s/${1}/${2}/g" ${file} > $$.tmp
	cat $$.tmp > ${file}
    done

find . -name "*.in" | \
    while read file; do
	sed -e "s/${1}/${2}/g" ${file} > $$.tmp
	cat $$.tmp > ${file}
    done

find . -name "*.bash" | \
    while read file; do
	sed -e "s/${1}/${2}/g" ${file} > $$.tmp
	cat $$.tmp > ${file}
    done
