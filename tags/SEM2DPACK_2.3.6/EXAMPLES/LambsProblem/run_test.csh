#!/bin/csh -f
# run_test.csh exit code 0 = succesful

# run SEM2DPACK
rm -f *_sem2d.*
sem2dsolve > info

# run comparison to reference solution (matlab script)
rm -f test.out
matlab -nosplash -nojvm < analyze_test.m > test.out

# test: exit code 0 = succesful
grep -sq "Test = 1" test.out  

# takes exit code from last command
# value in caller script is in $status
exit   
