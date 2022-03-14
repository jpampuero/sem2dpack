#!/bin/csh
set comps = (x y z)

foreach comp ($comps)
  set name = U${comp}_sem2d
  if (-e $name.dat) then
    echo 'Writing ascii file $name.tab'
    post_seis.exe << END > tmp
2
$comp
4
${name}.tab
1
END

  endif
end

rm -f tmp

echo 'NOTE: first column = time'
