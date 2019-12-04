      real function probbeam(ityp,iflag,bkg,acc)
*
*    File:  probbeam.f
* Purpose:  calculate single-beam and double-beam function values
*   Usage:  status = probbeam(ityp,iflag,bkg,acc,delc_val,trs,b4strob,
*                             bw1hrs,bw2hrs,bwtrs_val,b4trs_val)
*              where ityp = 1 for single-beam bkg and acc
*                    ityp = 2 for double-beam bkg and acc
*                   iflag = 0 returns bkg and acc given
*                             delc_val,trs-b4strob,bwhrs,bwtrs_val,b4tcut
*                   iflag = 1 returns acc given bkg
*                   iflag = 2 returns bkg given acc
*                  status = 1 if function is called correctly,
*                           0 otherwise.
*
      implicit none

      include 'pass3nt.inc'

      integer i,ityp,iflag
      real bkg,acc
      integer nrowbm1,nrowbm2,ncolbm1,ncolbm2
      parameter (nrowbm1 = 23, nrowbm2 = 16, ncolbm1 = 8, ncolbm2 = 7)
      real bm1(ncolbm1,nrowbm1),bm2(ncolbm2,nrowbm2)
      real delc_cut(nrowbm1),trsb4_cut(nrowbm1),
     &     bm1k(nrowbm1),bm1p(nrowbm1),cex(nrowbm1),bm1tot(nrowbm1),
     &     nbm1(nrowbm1),abm1(nrowbm1)
      real bwh_cut(nrowbm2),bwb4_cut(nrowbm2),
     &     bm2k(nrowbm2),bm2p(nrowbm2),bm2tot(nrowbm2),
     &     nbm2(nrowbm2),abm2(nrowbm2)
      real bkg1,bkg2,acc1,acc2,linint
      real bwhrs,offset,b4tcut,bwb4tcut
      logical first
      data first /.true./

*  Single-beam function
*
*      cut_position     bkg_k  bkg_pi bkg_cex  bkg_tot     N      A
*      DELC  trs-b4str  (E-2)  (E-2)  (E-2)    (E-2)
      data bm1/
     & 0.185, -99.00, 136.358, 21.647, 25.092, 183.097, 131.630, 1.104,
     & 0.265, -99.00,  67.793, 10.762, 18.928,  97.483,  70.081, 1.102,
     & 0.290, -99.00,  67.793, 10.762, 17.442,  95.997,  69.013, 1.102,
     & 0.330, -99.00,  33.649,  5.342, 15.038,  54.029,  38.842, 1.100,
     & 0.375, -99.00,  16.829,  2.672, 12.415,  31.916,  22.945, 1.098,
     & 0.415, -99.00,  16.829,  2.672, 10.448,  29.949,  21.531, 1.096,
     & 0.450, -99.00,   9.586,  1.522,  9.267,  20.375,  14.678, 1.093,
     & 0.515, -99.00,   4.478,  0.711,  6.951,  12.140,   8.748, 1.087,
     & 0.605, -99.00,   1.873,  0.297,  5.158,   7.328,   5.268, 1.076,
     & 0.640, -99.00,   1.152,  0.183,  4.197,   5.532,   3.977, 1.071,
     & 0.700, -99.00,   0.575,  0.091,  3.410,   4.076,   2.930, 1.061,
     & 1.000, -99.00,   0.144,  0.023,  1.224,   1.391,   1.000, 1.000,
     & 1.000,  -1.80,   0.144,  0.023,  1.224,   1.391,   1.000, 1.000,
     & 1.050,  -0.50,   0.139,  0.021,  1.093,   1.253,   0.901, 0.988,
     & 1.100,   0.40,   0.114,  0.015,  0.962,   1.091,   0.784, 0.976,
     & 1.200,   0.80,   0.084,  0.011,  0.831,   0.926,   0.666, 0.952,
     & 1.240,   1.25,   0.064,  0.005,  0.656,   0.725,   0.521, 0.943,
     & 1.300,   1.40,   0.045,  0.004,  0.525,   0.574,   0.412, 0.929,
     & 1.400,   1.60,   0.040,  0.002,  0.437,   0.479,   0.344, 0.906,
     & 1.450,   1.70,   0.025,  0.002,  0.350,   0.377,   0.271, 0.894,
     & 1.500,   2.00,   0.015,  0.001,  0.262,   0.278,   0.200, 0.883,
     & 1.800,   2.10,   0.010,  0.001,  0.131,   0.142,   0.102, 0.815,
     & 1.900,   2.20,   0.005,  0.001,  0.044,   0.050,   0.036, 0.793/

*  Double-beam function
*
*      cut_position  bkg_k   bkg_pi  bkg_tot      N      A
*     BWHRS BW/B4TRS (E-2)   (E-2)   (E-2)
      data bm2/
     & 0.00, 0.020, 11.222,   2.451,  13.673, 1243.000, 1.022,
     & 0.00, 0.050,  8.648,   1.945,  10.593,  963.000, 1.021,
     & 0.00, 0.120,  5.039,   1.158,   6.197,  563.364, 1.020,
     & 0.00, 0.190,  2.919,   0.595,   3.514,  319.455, 1.019,
     & 0.00, 0.330,  0.957,   0.148,   1.105,  100.455, 1.016,
     & 0.00, 0.380,  0.634,   0.086,   0.720,   65.455, 1.015,
     & 0.00, 0.500,  0.262,   0.019,   0.281,   25.546, 1.013,
     & 0.00, 0.600,  0.115,   0.008,   0.123,   11.182, 1.011,
     & 0.00, 0.800,  0.018,   0.003,   0.021,    1.909, 1.006,
     & 0.00, 1.000,  0.00918, 0.00227, 0.01145,  1.000, 1.000,
     & 1.85, 1.000,  0.00918, 0.00206, 0.01124,  0.892, 0.978,
     & 2.15, 1.000,  0.00734, 0.00186, 0.00920,  0.803, 0.975,
     & 2.50, 1.000,  0.00734, 0.00165, 0.00899,  0.785, 0.971,
     & 4.00, 1.000,  0.00550, 0.00165, 0.00715,  0.624, 0.950,
     & 6.00, 1.000,  0.00367, 0.00145, 0.00512,  0.447, 0.921,
     & 9.50, 1.000,  0.00367, 0.00124, 0.00491,  0.429, 0.872/

*
*  Initialize
*
      probbeam = 1.

      if (first) then
         first = .false.
         do i = 1,nrowbm1
            delc_cut(i) = bm1(1,i)
            trsb4_cut(i) = bm1(2,i)
            bm1k(i) = bm1(3,i)*0.01
            bm1p(i) = bm1(4,i)*0.01
            cex(i) = bm1(5,i)*0.01
            bm1tot(i) = bm1(6,i)*0.01
            nbm1(i) = bm1(7,i)
            abm1(i) = bm1(8,i)
         enddo
         do i = 1,nrowbm2
            bwh_cut(i) = bm2(1,i)
            bwb4_cut(i) = bm2(2,i)
            bm2k(i) = bm2(3,i)*0.01
            bm2p(i) = bm2(4,i)*0.01
            bm2tot(i) = bm2(5,i)*0.01
            nbm2(i) = bm2(6,i)
            abm2(i) = bm2(7,i)
         enddo
      endif

*
*  Calculate needed variables
*
      bwhrs = min(abs(bw1hrs),abs(bw2hrs))

      if (run.lt.29000) then       ! 95
         offset = -0.5
      else if (run.lt.35550) then  ! 96 - 97
         offset = 0.0
      else                         ! 98
         offset = 1.0
      endif
      b4tcut = abs((b4trs_val+offset)/3.)
      bwb4tcut = min(b4tcut,bwtrs_val)

*
*  Single-beam evalution
*
      if (ityp.eq.1) then
         if (iflag.eq.0) then
            if ((delc_val.lt.delc_cut(1)).or.
     &          ((trs-b4strob).lt.trsb4_cut(1))) then
               bkg = 9999.
               acc = 9999.
            else if ((delc_val.ge.delc_cut(nrowbm1)).and.
     &               ((trs-b4strob).ge.trsb4_cut(nrowbm1))) then
               bkg = 0.
               acc = 0.
            else
               bkg1 = linint(nbm1,delc_cut,nrowbm1,delc_val)
               bkg2 = linint(nbm1,trsb4_cut,nrowbm1,trs-b4strob)
               bkg = max(bkg1,bkg2)
               acc1 = linint(abm1,delc_cut,nrowbm1,delc_val)
               acc2 = linint(abm1,trsb4_cut,nrowbm1,trs-b4strob)
               acc = max(acc1,acc2)
            end if
         else if (iflag.eq.1)then
            if (bkg.gt.nbm1(1)) then
               acc = 9999.
            else if (bkg.le.nbm1(nrowbm1)) then
               acc = 0.
            else
               acc = linint(abm1,nbm1,nrowbm1,bkg)
            endif
         else if (iflag.eq.2)then
            if (acc.gt.abm1(1)) then
               bkg = 9999.
            else if (acc.le.abm1(nrowbm1)) then
               bkg = 0.
            else
               bkg = linint(nbm1,abm1,nrowbm1,acc) 
            endif
         else
            probbeam = 0.
            write (6,*) 'Wrong input iflag!'
            stop
         end if

*
*  Double-beam evalution
*
      else if (ityp.eq.2) then
         if (iflag.eq.0) then
            if ((bwhrs.lt.bwh_cut(1)).or.
     &          (bwb4tcut.lt.bwb4_cut(1))) then
               bkg = 9999.
               acc = 9999.
            else if ((bwhrs.ge.bwh_cut(nrowbm2)).and.
     &               (bwb4tcut.ge.bwb4_cut(nrowbm2))) then
               bkg = 0.
               acc = 0.
            else
               bkg1 = linint(nbm2,bwh_cut,nrowbm2,bwhrs)
               bkg2 = linint(nbm2,bwb4_cut,nrowbm2,bwb4tcut)
               bkg = max(bkg1,bkg2)
               acc1 = linint(abm2,bwh_cut,nrowbm2,bwhrs)
               acc2 = linint(abm2,bwb4_cut,nrowbm2,bwb4tcut)
               acc = max(acc1,acc2)
            end if
         else if (iflag.eq.1)then
            if (bkg.gt.nbm2(1)) then
               acc = 9999.
            else if (bkg.le.nbm2(nrowbm2)) then
               acc = 0.
            else
               acc = linint(abm2,nbm2,nrowbm2,bkg)
            endif
         else if (iflag.eq.2)then
            if (acc.gt.abm2(1)) then
               bkg = 9999.
            else if (acc.le.abm2(nrowbm2)) then
               bkg = 0.
            else
               bkg = linint(nbm2,abm2,nrowbm2,acc) 
            endif
         else
            probbeam = 0.
            write (6,*) 'Wrong input iflag!'
            stop
         end if

      else
         probbeam = 0.
         write (6,*) 'The input ityp value does not match 1 or 2'
         stop         
      end if

      return
      end

      include 'linint.f'

