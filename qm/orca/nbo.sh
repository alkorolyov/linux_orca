name="input"
janpa_dir="~/linux_qm/qm/janpa"
orca_2mkl $name -molden
java -jar $janpa_dir/molden2molden.jar -i $name.molden.input -o $name.PURE -fromorca3bf -orca3signs  > M2M.log
java -jar $janpa_dir/janpa.jar -i $name.PURE  > $name.JANPA