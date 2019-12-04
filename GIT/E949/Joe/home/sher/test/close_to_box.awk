#
#    File:  close_to_box.awk
# Purpose:  count events that are close to the box
#   Usage:  > awk -f close_to_box.awk close_to_box.{skim}.{year}.dat
#
# History:  2001 Sep  PB  created
#

BEGIN {
# final box cut
   nfpv = 1.0;
   nfkp2 = 0.00147;
   nftd = 1.0;
   nfkm2t = 0.1346800E-01;
   nfkm2b = 0.2964120E-01;
   nfbm1 = 1.0;
   nfbm2 = 1.0;
   factor = 5.0;
}

{
   if (($3<factor*nfpv)&&
       ($4<factor*nfkp2)&&
       ($5<factor*nftd)&&
       ($6<factor*nfkm2t)&&
       ($7<factor*nfkm2b)&&
       ($8<factor*nfbm1)&&
       ($9<factor*nfbm2)) {
      print;
   }
}

END {
}
