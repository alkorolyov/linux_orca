if [ ! -e "$1" ] ; then
  echo "No input file provided"
  exit 1
fi



java -jar molden2molden.jar -i $1 -o output.PURE -fromorca3bf -orca3signs > M2M.log
java -jar janpa.jar -i output.PURE  > output.JANPA