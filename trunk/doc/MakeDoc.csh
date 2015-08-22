#!/bin/csh -f
# Prints the documentation for the input blocks of the SEM2D sources

echo "      ============================================================="
echo "      = Self-documentation for the INPUT BLOCKS of the SEM2D code ="
echo "      ============================================================="
echo " "
echo " "
#foreach file ($0:h/../SRC/*.f90) 
foreach file (../SRC/*.f90) 
  # In the source code the input blocks are comment lines
  # bounded by "BEGIN INPUT BLOCK" and "END INPUT BLOCK"
  sed -n '/BEGIN INPUT BLOCK/,/END INPUT BLOCK/p' $file \
   | sed -e '/END INPUT BLOCK/d' \
   -e 's/BEGIN INPUT BLOCK/----------------------------------------------------------------------------/g' \
   -e 's/^\!//g' -e 's/ARG://g'
end
