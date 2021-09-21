#!/bin/csh
set comps = (x y z)

foreach comp ($comps)
  set name = U${comp}_sem2d
  if (-e $name.dat) then
    echo $name
    post_seis.exe << END > tmp
2
$comp
5
$name
3
2.0
1
END

    gv $name.ps &
  endif
end

rm -f tmp
