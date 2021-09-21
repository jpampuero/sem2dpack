#!/bin/csh
# Breaks a PostScript article into pages
# and creates a .tex to include in your .tex

set doc = $1
set NPAGES = `grep %%Pages: $doc | tail -1 | awk '{print $2}'`
set tex = $doc:t.tex

echo "% include this file in your .tex :" > $tex 
echo "%   \input{$tex}" >> $tex 
@ page = 0
while ( $page < $NPAGES ) 
  @ page ++ 
  set pagename = $doc:t:r_$page.ps
  psselect -p$page $doc  $pagename
  #psselect -p$page $doc tmp.ps
  #set pagename = $doc:t:r_$page.epsi
  #ps2epsi tmp.ps $pagename
  # \ImgC is my own LaTeX macro to insert figures
  echo "\ImgC{$pagename}{1}" >> $tex
end
