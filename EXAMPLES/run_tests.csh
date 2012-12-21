#!/bin/csh -f
# runs a test in each directory containing a matlab script named analyze_test.m

echo
echo "------------------------------"
echo "-  Test suite for SEM2DPACK  -"
echo "------------------------------"
echo

@ ntotal = 0
@ nok = 0
set start=`pwd`
foreach name (`ls`)
  if (! -d $name) continue
  if (! -e $name/analyze_test.m) continue

  @ ntotal ++
  cd $name 
  echo -n $name "............"
    
 # run test 
  rm -f *_sem2d.* test.out
  sem2dsolve > info
  matlab -nosplash -nojvm < analyze_test.m > test.out
  grep -sq "Test = 1" test.out  
 # exit code from last command $status==0 means succesful

  if ($status == 0) then
    @ nok ++
    echo "......... [OK]"
  else
    echo "......... [FAILED]"
  endif

#  grep -q "FLEXlm error" test.out  
#  if ($status == 0) then
#  else
#    echo "status = " $status
#    echo "Matlab license error !"
#    cat test.out
#    exit 1
#  endif

  cd $start
end

echo
echo "Summary: "
if ($nok == $ntotal) then
  echo "[OK] Passed each of " $ntotal " tests"
else
  echo "[FAILED] Only passed " $nok " tests out of " $ntotal 
endif

echo
