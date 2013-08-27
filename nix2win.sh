#!/bin/bash 
# A simple script to convert a unix file to Windows format
f=$1
f_nix="${f}_nix"
#echo "mv $f $f_nix"
mv $f $f_nix
#echo "awk 'sub("$", "\r")' $f_nix > $f"
awk 'sub("$", "\r")' $f_nix > $f
rm $f_nix