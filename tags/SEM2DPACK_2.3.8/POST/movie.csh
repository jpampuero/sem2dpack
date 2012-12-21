#!/bin/csh
# Generates a movie from *_*_sem2d.ps (output of SEM2D)

#-- check arguments
if ( $#argv == 0 | $#argv > 2 ) then
  echo "movie.csh [-crop] xx"
  echo "where xx = dx, vx or ax, etc"
  echo "makes a movie_xx_sem2d.gif from xx_*_sem2d.ps" 
  echo "-crop crops tightly the image, leaving out text"
  exit
endif
  @ do_crop = ($1 == '-crop')
if ( $#argv > 1 ) then
  set xx = $2
else
  set xx = $1
endif

#-- check gif packages
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

#-- make movie
echo "  Generating a movie"
#@ lxo = 562
@ lxo = 533
@ lzo = 588
@ pxo = 37
@ pzo = 158
@ marge = 20

@ i = 1
foreach ps_file ( `/bin/ls ${xx}_*.ps` )

  echo "Converting $ps_file..."
  set Num = $ps_file:r

  if ($do_crop) then
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

  else
    convert -rotate 90 -resize 640x480 $ps_file $Num.gif
    #convert -rotate 90 -resize 320x240 $ps_file $Num.gif
  endif

end
rm -f tmp.gif

#gifmerge -l0 Snapshot*.gif > movie_sem2d.gif
#gifsicle -m -d 10 -l Snapshot*.gif > movie_sem2d.gif
$MERGE_GIF $MERGE_GIF_OPT ${xx}_*.gif > movie_${xx}_sem2d.gif
rm -f ${xx}_*.gif
echo "Created movie_${xx}_sem2d.gif"
echo "To watch it again: animate movie_${xx}_sem2d.gif"
animate movie_${xx}_sem2d.gif &
#xanim movie_sem2d.gif &
