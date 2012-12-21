#!/bin/csh
post_seis.exe << END > tmp
2
6
Ux_sem2d.ps
3
2.0
1
END

post_seis.exe << END > tmp
3
6
Uz_sem2d.ps
3
2.0
1
END

gv ux_sem2d.ps &
gv uz_sem2d.ps &

rm -f tmp
