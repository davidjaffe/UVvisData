#
#    File:  bgtotal.awk
# Purpose:  locate events in cells of the total background function
#   Usage:  > awk -f bgtotal.awk bgtotal.dat
#
# History:  2001 Sep  PB  created
#

BEGIN {
#  1998 candidate event (run 37443, event 176348)
#   nfpv = 0.1953363;
#   nfkp2 = 0.0002985278;
#   nftd = 0.0;
#   nfkm2t = 0.0;
#   nfkm2b = 0.007959974
#   nfbm1 = 0.3383143;
#   nfbm2 = 0.0;

#  1998 signal-like event (run 39180, event 56327)
#   nfpv = 0.1343
#   nfkp2 = 0.0002985278;
#   nftd = 3.1010;
#   nfkm2t = 0.001700;
#   nfkm2b = 0.007020;
#   nfbm1 = 0.0;
#   nfbm2 = 0.0;

#  1998 L/R event (run 37761, event 142498)
   nfpv = 1.3551;
   nfkp2 = 0.2810E-02;
   nftd = 0.0000;
   nfkm2t = 0.9199E-03;
   nfkm2b =  0.000;
   nfbm1 = 0.0411;
   nfbm2 = 0.9401
}

{
   if (($6>nfpv)&&
       ($7>nfkp2)&&
       ($8>nftd)&&
       ($9>nfkm2t)&&
       ($10>nfkm2b)&&
       ($11>nfbm1)&&
       ($12>nfbm2)) {
      print NR-1" "$0;
   }
}

END {
}
