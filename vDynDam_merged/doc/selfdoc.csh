#!/bin/csh -f
echo 'Running selfdoc.csh ...'
./MakeDoc.csh > doc.txt
echo "\\begin{verbatim}" > tmp1
echo "\\end{verbatim}" > tmp2
cat tmp1 doc.txt tmp2 > USERS_GUIDE/selfdoc.tex

rm -f tmp1 tmp2
