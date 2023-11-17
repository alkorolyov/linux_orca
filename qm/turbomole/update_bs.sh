#!/bin/bash

# create new folder
_TMP_DIR="$TMP_DIR"0
mkdir $_TMP_DIR
# copy coords and hessian
cp $TMP_DIR/hessapprox $_TMP_DIR/hessapprox
cp $TMP_DIR/coord $_TMP_DIR/coord
cd $_TMP_DIR
export TMP_DIR=$_TMP_DIR

# define with new basis set
define <<EOF


a coord
*
no
b all $1
*
eht



*
EOF

match='$grad    file=gradient'
insert='$hessapprox   file=hessapprox'
file='control'
sed -i "s/$match/$match\n$insert/" $file

#uff
jobex