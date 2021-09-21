#!/bin/csh

set Path = $0:h/
set Suffix = ".emc2_bd"

echo '\n'
echo '  HREF   : H-REFINEMENT OF AN EMC2 MESH'
echo '  WARNING: all file names will be automatically given'
echo '           the .emc2_bd suffix'
echo ' '

echo -n '  Enter the original file name : '
set file1 = $<
set file1 = ${file1}$Suffix
if !( -e $file1 ) then
  echo '  FATAL ERROR: file does not exist'
  echo '  HREF: abort'
  echo '\n'
  exit
endif

echo -n '  Enter the new mesh file name : '
set file2 = $<
set file2 = ${file2}$Suffix

echo -n '  Refinement ratio (integer>1) : '
set ratio = $<

awk -f ${Path}href.awk ratio=$ratio $file1 > $file2

echo ' '
echo '  You can restore now the new mesh in EMC2'
echo '  (in PREP_MESH mode) and save it in FTQ format'
echo '  HREF: end'
echo '\n'
