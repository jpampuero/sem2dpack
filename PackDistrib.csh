#!/bin/tcsh -f
# Wraps the SEM2DPACK package for public distribution

set Sem2DPath =  $0:h
set SrcDir    = $Sem2DPath/SRC

cd $Sem2DPath
echo Preparing a distribution package for:
cat CurrentVersion

set Version   = `head -1 CurrentVersion | awk '{print $3}'`
set NewPackDir  = ../SEM2DPACK_$Version

echo Last packaged version was:
ls -lth ../repository/*gz | head -1
echo

echo Checklist:
cat PackDistrib.checklist

# create a dummy dir, I'll put there all the files to be distributed
/bin/rm -rf $NewPackDir
mkdir $NewPackDir

# make a general file header
cat CurrentVersion Copyright AuthorContact LicenseNotice > header.tmp

# process the source files: add headers
echo -n Processing the source files ...
mkdir $NewPackDir/SRC
sed 's/^/! /g' header.tmp > f90_header.tmp
foreach f90file ($SrcDir/*.f90)
  cat f90_header.tmp $f90file >  $NewPackDir/SRC/$f90file:t
end
sed 's/^/# /g' header.tmp > make_header.tmp
foreach makefile ($SrcDir/Makefile*)
  cat make_header.tmp $makefile > $NewPackDir/SRC/$makefile:t
end
echo ... OK
 
# make documentation
make -C doc/USERS_GUIDE/

# process general files
echo -n Processing general files ...
cp header.tmp $NewPackDir/ReadMeLicense
cp ChangeLog CurrentVersion ReadMeInstall GeneralPublicLicense ToDo Copyright AuthorContact $NewPackDir
cp doc/USERS_GUIDE/users_guide_sem2dpack.pdf $NewPackDir
cp -R EXAMPLES $NewPackDir
# exclude output files from the sample directories
/bin/rm -f $NewPackDir/EXAMPLES/*/*_sem2d*
/bin/rm -f $NewPackDir/EXAMPLES/*/*.gif
/bin/rm -rf $NewPackDir/EXAMPLES/*/data_*
cp -R POST $NewPackDir
cp -R PRE $NewPackDir
echo ... OK

echo -n Packaging and compressing sem2dpack_$Version.tar.gz ...
find $NewPackDir/ -name "*~" -exec /bin/rm -f {} \;
tar cfhz ../repository/sem2dpack_$Version.tar.gz -C .. SEM2DPACK_$Version
echo ... OK

/bin/rm -f *.tmp

echo Checklist:
cat PackDistrib.checklist

