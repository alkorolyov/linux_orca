#!/bin/bash

# Run ORCA calculation with provided options
#
# syntax ./run.sh structure.xyz "orca options"(optional) solvent(optional)
# example usage:
# ./run.sh methanol.xyz
# ./run.sh methanol.xyz 'B3LYP DEF2-TZVP OPT'
# ./run.sh methanol.xyz 'MP2 DEF2-TZVP OPT' THF

tmp_dir=.tmp
xyz_file="$1"

# Check if the second argument is provided
if [ -z "$2" ]; then
    options="B3LYP DEF2-SVP OPT" # default
else
    options="$2"
fi

# Check if the third argument is provided
if [ -n "$3" ]; then
#  options="$options CPCM($3)"
  solvent="%cpcm
  smd true
  SMDsolvent \"$3\"
  end"
fi
inp_file="input.inp"


# Check if the XYZ file exists
if [ ! -e "$xyz_file" ]; then
    echo "Error: XYZ file not found."
    exit 1
fi

# Create the input file
cat <<EOF >"$inp_file"
!$options
%pal nprocs 8 end
%geom
   MaxIter 100
   end
$solvent
* xyz 0 1
EOF
tail -n +3 "$xyz_file" >> "$inp_file"
cat <<EOF >>"$inp_file"
*

EOF

echo "Input file '$inp_file' generated successfully."

# Create tmp dir and copy inputs
rm -rf $tmp_dir ; mkdir $tmp_dir
mv $inp_file $tmp_dir
#cp $xyz_file $tmp_dir
cd $tmp_dir

echo "Running orca ... "
echo "Check .tmp/output for progress and results"
/opt/orca-5.0.3/orca/orca input.inp --use-hwthread-cpus > output
#tail -f output | grep --line-buffered -E "GEOMETRY OPTIMIZATION CYCLE|TOTAL RUN TIME"
echo "Finished"