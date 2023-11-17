#!/bin/bash

tmp_dir=.tmp
xyz_file="$1"
inp_file="input.inp"

# Check if the XYZ file exists
if [ ! -e "$xyz_file" ]; then
    echo "Error: XYZ file not found."
    exit 1
fi

# Create the input file
cat <<EOF >"$inp_file"
!B3LYP DEF2-SVP OPT
%output
  PrintLevel Mini
  end
%pal nprocs 24 end
%geom
   MaxIter 100
   end
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

echo "Running orca. Check /orca/.tmp/output for progess ..."
/opt/orca-5.0.3/orca/orca input.inp --use-hwthread-cpus > output

echo "Finished"