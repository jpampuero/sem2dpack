# TASK: h-refinement of an emc2 mesh
# USAGE: gawk -f href.awk ratio=(integer >1) original_file.emc2_bd > refined_file.emc2_bd
# INPUT: refinement ratio (integer >1) supplied at command line
#        + mesh file *.emc2_bd
# OUTPUT: modified mesh file *.emc2_bd
#       The output file must be opened in emc2 and saved to *.ftq format
# 
BEGIN { NFpre = 0
}
{ if ( $1 == "'SEGMENT'" || $1 == "'SPLINE'" ) 
   { if (NF >= 9) 
      {for (i=1;i<=8;i++) printf("%s ",$i) 
       printf("%s ", ratio*($9-1) +1)
       for (i=10;i<=NF;i++) printf("%s ",$i)
       print " "
      }
    else
      {NFpre = NF;print$0; next}
   }
  else 
   { if ( NFpre == 0 )
      print $0
    else
      {pos = 9-NFpre
       for (i=1;i<=pos-1;i++) printf("%s ",$i)
       printf("%s ", ($pos-1)*ratio +1 )
       for (i=pos+1;i<=NF;i++) printf("%s ",$i)
       print " "
       NFpre = 0 
      }
   }
}
END{
}
