#!/bin/csh
# Generates a movie from Snapshot*.ps (output of SEM2D)

if ( (! -X convert) | (! -X animate) ) then
  echo "  This script requires the ImageMagick package (convert and animate)"
  exit 1
endif
if (-X gifsicle) then
  set MERGE_GIF = gifsicle
  set MERGE_GIF_OPT =  "-m -d 10 -l"
else if (-X gifmerge) then 
  set MERGE_GIF = gifmerge
  set MERGE_GIF_OPT =  -l0
else
  echo "  This script requires the GIFSICLE or GIFMERGE packages"
  exit 1
endif

echo "  Generating a movie"
#@ lxo = 562
@ lxo = 533
@ lzo = 588
@ pxo = 37
@ pzo = 158
@ marge = 20

@ i = 1
foreach ps_file ( `/bin/ls Snapshot*.ps` )

  echo "Converting $ps_file..."
  set Num = $ps_file:r
  # crop out the text
  convert -crop ${lxo}x${lzo}+${pxo}+${pzo} $ps_file tmp.gif

  # get effective image size for tight crop
  if ($i == 1) then

    convert -trim tmp.gif tmp3.gif
    set nx = `identify  -format "%w" tmp3.gif`
    set nz = `identify  -format "%h" tmp3.gif `
    @ nx = $nx + $marge + $marge
    @ nz = $nz + $marge + $marge
    @ px = $lxo - $nx
    @ pz = $lzo - $nz
    @ i = 0
    rm -f tmp3.gif

  endif

  convert -crop ${nx}x${nz}+${px}+${pz} -geometry 100% -rotate 90 tmp.gif $Num.gif

end
rm -f tmp.gif

#gifmerge -l0 Snapshot*.gif > movie_sem2d.gif
#gifsicle -m -d 10 -l Snapshot*.gif > movie_sem2d.gif
$MERGE_GIF $MERGE_GIF_OPT Snapshot*.gif > movie_sem2d.gif
rm -f Snapshot*.gif
echo "Your movie is in movie_sem2d.gif"
echo "To watch it again: animate movie_sem2d.gif"
animate movie_sem2d.gif &
#xanim movie_sem2d.gif &
