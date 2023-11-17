#!/bin/bash

xyz_file="$1"
basis_set="$2"

# Check if the XYZ file exists
if [ ! -e "$xyz_file" ]; then
    echo "Error: XYZ file not found."
    exit 1
fi

# Check basis set supplied
if [ -z "$2" ]; then
    echo "No basis set supplied"
    exit 1
fi

# create tmp dir or delete existing
if [ ! -e "$TMP_DIR" ]
then
    mkdir $TMP_DIR
else
    rm -rf $TMP_DIR
fi


# copy coords
x2t "$xyz_file" > $TMP_DIR/coord
cd $TMP_DIR

define <<EOF


a coord
*
no
b all $basis_set
*
eht



dft
on
func b3-lyp
*
*
EOF

# single point
#dscf coord

# geometry optimization
uff
jobex