#!/bin/csh
post_seis.exe << END > tmp
2
5
Ux_sem2d.tab
1
END


post_seis.exe << END > tmp
3
5
Uz_sem2d.tab
1
END

rm -f tmp

echo 'Wrote ascii files Ux_sem2d.tab and Uz_sem2d.tab'
echo 'NOTE: the first column is the time axis'
