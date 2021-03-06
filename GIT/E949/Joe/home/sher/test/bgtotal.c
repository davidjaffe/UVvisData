/*
 *    File:  bgtotal.c
 * Purpose:  define the inside-the-box (and outside-the-box)
 *           global background function
 *   Usage:  > bgtotal
 *
 * History:  2001 Jul 25  PB  created
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bgtotal.h"

/*
 *  Filenames for the ascii data files which contain N and A values
 *  for the 7 background functions.  Note that the functions are designed
 *  such that:
 *  - the pv and kp2 functions are linear in Kp2 background;
 *  - the td function is linear in both Km2 range-tail and
 *    muon-band background;
 *  - the km2t function is linear in Km2 range-tail background;
 *  - the km2b function is linear in muon-band background;
 *  - the bm1 function is linear in the sum of single-beam kaon- and
 *    pion-entering background plus charge-exchange background;
 *  - the bm2 function is linear in the sum of double-beam kaon- and
 *    pion-entering background.
 *  This is why there are 7 functions but 8 background types.
 *  The code below is written specifically with these particular
 *  functions in mind.  The code is also written on the assumption
 *  that N and A values are tabulated in decreasing order.
 */
char FilePv[15] = "pvfunc.dat";
char FileKp2[15] = "kp2func.dat";
char FileTd[15] = "tdfunc.dat";
char FileKm2t[15] = "km2tfunc.dat";
char FileKm2b[15] = "km2bfunc.dat";
char FileBm1[15] = "bm1func.dat";
char FileBm2[15] = "bm2func.dat";

/*
 *  Expected background levels at the final box position for the
 *  8 background types.
 */
float BoxBgkp2 = 0.012;
float BoxBgkm2t = 0.0245;
float BoxBgkm2b = 0.0092;
float BoxBgbm1 = 0.0196;
/*  bm1 = bm1k + bm1p + cex
float BoxBgbm1k = 0.0034;
float BoxBgbm1p = 0.0005;
float BoxBgcex = 0.0157;
*/
float BoxBgbm2 = 0.0004;
/*  bmm2 = bm2k + bm2p
float BoxBgbm2k = 0.00032;
float BoxBgbm2p = 0.00008;
*/

/*
 *  Function (N and A) values at the final box position for the 7 functions
 *  (bm1 and bm2 N values are also broken down into kaon, pion, and cex
 *  values).  These must be EXACTLY the same as the correponding values in
 *  the data files above.
 */
float BoxNpv = 1.0;
float BoxNkp2 = 0.00147;
float BoxNtd = 1.0;
float BoxNkm2t = 0.1346800E-01;
float BoxNkm2b = 0.2964120E-01;
float BoxNbm1 = 1.0;
/*  if Nbm1k, Nbm1p, and Ncex are known for each Nbm1,
    these might be useful.
float BoxNbm1k = 0.00144;
float BoxNbm1p = 0.00023;
float BoxNcex = 0.01224;
*/
float BoxNbm2 = 1.0;
/*  if Nbm2k and Nbm2p are known for each Nbm2,
    these might be useful.
float BoxNbm2k = 0.0000918;
float BoxNbm2p = 0.0000227;
*/

float BoxApv = 0.90020;
float BoxAkp2 = 1.0;
float BoxAtd = 1.0;
float BoxAkm2t = 1.0;
float BoxAkm2b = 1.0;
float BoxAbm1 = 1.0;
float BoxAbm2 = 1.0;

/*
 *  Start of main program.
 */
int main(void)
{
   FILE *ifp,*ofp,*ofpin,*ofpout;
   DKFUNC bgpv,bgkp2,bgtd,bgkm2t,bgkm2b;
   BMFUNC bgbm1,bgbm2;
   TOTFUNC bgtotal;
   PTRDKVALUE p,p1,p2,p3,p4,p5;
   PTRBMVALUE b,b1,b2;
   PTRTOTVALUE t,t1,t2a,t2b,t2c;
   float SNmax,kinacc;
   float dBgkp2,dBgkm2t,dBgkm2b,dBgbm1,dBgbm2;
   float dNtotal,dAtotal;
   float Nlo,Alo,Nin,Ain,Ninlo,Ainlo,Nout,Aout,Noutlo,Aoutlo;
   int i,first,ninside,noutside,firstin,firstout;
   float d;

   /*
    *  Initialize.
    */
   bgpv.nvalues = 0;
   bgpv.head = bgpv.box = bgpv.tail = NULL;
   bgkp2.nvalues = 0;
   bgkp2.head = bgkp2.box = bgkp2.tail = NULL;
   bgtd.nvalues = 0;
   bgtd.head = bgtd.box = bgtd.tail = NULL;
   bgkm2t.nvalues = 0;
   bgkm2t.head = bgkm2t.box = bgkm2t.tail = NULL;
   bgkm2b.nvalues = 0;
   bgkm2b.head = bgkm2b.box = bgkm2b.tail = NULL;
   bgbm1.nvalues = 0;
   bgbm1.head = bgbm1.box = bgbm1.tail = NULL;
   bgbm2.nvalues = 0;
   bgbm2.head = bgbm2.box = bgbm2.tail = NULL;
   bgtotal.nvalues = 0;
   bgtotal.head = bgtotal.tail = NULL;

   /*
    *  Read the PV function into a structure.
    */
   if ((ifp = fopen(FilePv,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using PV function, because",
             FilePv,"does not exist.");
   } else {
      bgpv.head = malloc(sizeof(DKVALUE));
      bgpv.head -> prev = NULL;
      bgpv.tail = bgpv.head;
      p = bgpv.head;
      while (fscanf(ifp,"%f %f",&(p->n),&(p->a)) == 2) {
         (p -> n) = 46.918/(p -> n);
         bgpv.nvalues += 1;
         if ((p -> n) == BoxNpv) {
            bgpv.box = p;
         }
         bgpv.tail = p;
         bgpv.tail -> next = malloc(sizeof(DKVALUE));
         bgpv.tail -> next -> prev = bgpv.tail;
         p = bgpv.tail -> next;
      }
      p -> n = 0;
      p -> a = 0;
      bgpv.nvalues += 1;
      bgpv.tail = p;
      bgpv.tail -> next = NULL;
   }

   /*
    *  Read the Kp2 kinematic function into a structure.
    */
   if ((ifp = fopen(FileKp2,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using Kp2 kinematic function, because",
             FileKp2,"does not exist.");
   } else {
      bgkp2.head = malloc(sizeof(DKVALUE));
      bgkp2.head -> prev = NULL;
      bgkp2.tail = bgkp2.head;
      p = bgkp2.head;
      while (fscanf(ifp,"%f %f %f %f %f",&d,&d,&d,&(p->a),&(p->n)) == 5) {
         bgkp2.nvalues += 1;
         if ((p -> n) == BoxNkp2) {
            bgkp2.box = p;
         }
         bgkp2.tail = p;
         bgkp2.tail -> next = malloc(sizeof(DKVALUE));
         bgkp2.tail -> next -> prev = bgkp2.tail;
         p = bgkp2.tail -> next;
      }
      p -> n = 0;
      p -> a = 0;
      bgkp2.nvalues += 1;
      bgkp2.tail = p;
      bgkp2.tail -> next = NULL;
   }

   /*
    *  Read the TD function into a structure.
    */
   if ((ifp = fopen(FileTd,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using TD function, because",
             FileTd,"does not exist.");
   } else {
      bgtd.head = malloc(sizeof(DKVALUE));
      bgtd.head -> prev = NULL;
      bgtd.tail = bgtd.head;
      p = bgtd.head;
      while (fscanf(ifp,"%f %f %f",&(p->n),&(p->a),&d) == 3) {
         (p -> n) /= 39;
         bgtd.nvalues += 1;
         if ((p -> n) == BoxNtd) {
            bgtd.box = p;
         }
         bgtd.tail = p;
         bgtd.tail -> next = malloc(sizeof(DKVALUE));
         bgtd.tail -> next -> prev = bgtd.tail;
         p = bgtd.tail -> next;
      }
      p -> n = 0;
      p -> a = 0;
      bgtd.nvalues += 1;
      bgtd.tail = p;
      bgtd.tail -> next = NULL;
   }

   /*
    *  Read the Km2 range-tail kinematic function into a structure.
    */
   if ((ifp = fopen(FileKm2t,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using Km2 range-tail kinematic function, because",
             FileKm2t,"does not exist.");
   } else {
      bgkm2t.head = malloc(sizeof(DKVALUE));
      bgkm2t.head -> prev = NULL;
      bgkm2t.tail = bgkm2t.head;
      p = bgkm2t.head;
      while (fscanf(ifp,"%f %f",&(p->n),&(p->a)) == 2) {
         bgkm2t.nvalues += 1;
         if ((p -> n) == BoxNkm2t) {
            bgkm2t.box = p;
         }
         bgkm2t.tail = p;
         bgkm2t.tail -> next = malloc(sizeof(DKVALUE));
         bgkm2t.tail -> next -> prev = bgkm2t.tail;
         p = bgkm2t.tail -> next;
      }
      p -> n = 0;
      p -> a = 0;
      bgkm2t.nvalues += 1;
      bgkm2t.tail = p;
      bgkm2t.tail -> next = NULL;
   }

   /*
    *  Read the muon-band kinematic function into a structure.
    */
   if ((ifp = fopen(FileKm2b,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using muon-band kinematic function, because",
             FileKm2b,"does not exist.");
   } else {
      bgkm2b.head = malloc(sizeof(DKVALUE));
      bgkm2b.head -> prev = NULL;
      bgkm2b.tail = bgkm2b.head;
      p = bgkm2b.head;
      while (fscanf(ifp,"%f %f",&(p->n),&(p->a)) == 2) {
         bgkm2b.nvalues += 1;
         if ((p -> n) == BoxNkm2b) {
            bgkm2b.box = p;
         }
         bgkm2b.tail = p;
         bgkm2b.tail -> next = malloc(sizeof(DKVALUE));
         bgkm2b.tail -> next -> prev = bgkm2b.tail;
         p = bgkm2b.tail -> next;
      }
      p -> n = 0;
      p -> a = 0;
      bgkm2b.nvalues += 1;
      bgkm2b.tail = p;
      bgkm2b.tail -> next = NULL;
   }

   /*
    *  Read the single-beam function into a structure.
    */
   if ((ifp = fopen(FileBm1,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using single-beam function, because",
             FileBm1,"does not exist.");
   } else {
      bgbm1.head = malloc(sizeof(BMVALUE));
      bgbm1.head -> prev = NULL;
      bgbm1.tail = bgbm1.head;
      b = bgbm1.head;
      while (fscanf(ifp,"%f %f %f %f %f %f %f",&d, \
             &(b->kaon),&(b->pion),&(b->cex),&d,&(b->n),&(b->a)) == 7) {
         b -> kaon *= 0.01;
         b -> pion *= 0.01;
         b -> cex *= 0.01;
         bgbm1.nvalues += 1;
         if ((b -> n) == BoxNbm1) {
            bgbm1.box = b;
         }
         bgbm1.tail = b;
         bgbm1.tail -> next = malloc(sizeof(BMVALUE));
         bgbm1.tail -> next -> prev = bgbm1.tail;
         b = bgbm1.tail -> next;
      }
      b -> kaon = 0;
      b -> pion = 0;
      b -> cex = 0;
      b -> n = 0;
      b -> a = 0;
      bgbm1.nvalues += 1;
      bgbm1.tail = b;
      bgbm1.tail -> next = NULL;
   }

   /*
    *  Read the double-beam function into a structure.
    */
   if ((ifp = fopen(FileBm2,"r")) == NULL) {
      printf("%s \n %s %s \n",
             "Not using double-beam function, because",
             FileBm2,"does not exist.");
   } else {
      bgbm2.head = malloc(sizeof(BMVALUE));
      bgbm2.head -> prev = NULL;
      bgbm2.tail = bgbm2.head;
      b = bgbm2.head;
      while (fscanf(ifp,"%f %f %f %f %f %f",&d, \
             &(b->kaon),&(b->pion),&d,&(b->n),&(b->a)) == 6) {
         b -> cex = 0.;
         b -> kaon *= 0.01;
         b -> pion *= 0.01;
         bgbm2.nvalues += 1;
         if ((b -> n) == BoxNbm2) {
            bgbm2.box = b;
         }
         bgbm2.tail = b;
         bgbm2.tail -> next = malloc(sizeof(BMVALUE));
         bgbm2.tail -> next -> prev = bgbm2.tail;
         b = bgbm2.tail -> next;
      }
      b -> kaon = 0;
      b -> pion = 0;
      b -> cex = 0;
      b -> n = 0;
      b -> a = 0;
      bgbm2.nvalues += 1;
      bgbm2.tail = b;
      bgbm2.tail -> next = NULL;
   }

   /*
    *  Check that the PV function is
    *  constructed properly.
    */
   for (p = bgpv.head; p -> next != NULL; p = p -> next) {
      p -> dn = ((p -> n) - (p -> next -> n))/BoxNpv;
      p -> da = ((p -> a) - (p -> next -> a))/BoxApv;
      if ((p -> dn <= 0)||(p -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the PV function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                p -> n, "Exiting...");
         exit(1);
      }
      p -> SN = (p -> da)/(p -> dn);
   }
   bgpv.tail -> dn = 0;
   bgpv.tail -> da = 0;
   bgpv.tail -> SN = 0;
   SNmax = 0;
   for (p = bgpv.head; p -> next != NULL; p = p -> next) {
      if (p -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the PV function are",
                "not increasing as the function is",
                "tightened, at N =",
                p -> n, "Exiting...");
         exit(2);
      }
      SNmax = p -> SN;
   }

   /*
    *  Check that the Kp2 kinematic function is
    *  constructed properly.
    */
   for (p = bgkp2.head; p -> next != NULL; p = p -> next) {
      p -> dn = ((p -> n) - (p -> next -> n))/BoxNkp2;
      p -> da = ((p -> a) - (p -> next -> a))/BoxAkp2;
      if ((p -> dn <= 0)||(p -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the Kp2 kinematic function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                p -> n, "Exiting...");
         exit(1);
      }
      p -> SN = (p -> da)/(p -> dn);
   }
   bgkp2.tail -> dn = 0;
   bgkp2.tail -> da = 0;
   bgkp2.tail -> SN = 0;
   SNmax = 0;
   for (p = bgkp2.head; p -> next != NULL; p = p -> next) {
      if (p -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the Kp2 kinematic function are",
                "not increasing as the function is",
                "tightened, at N =",
                p -> n, "Exiting...");
         exit(2);
      }
      SNmax = p -> SN;
   }

   /*
    *  Check that the TD function is
    *  constructed properly.
    */
   for (p = bgtd.head; p -> next != NULL; p = p -> next) {
      p -> dn = ((p -> n) - (p -> next -> n))/BoxNtd;
      p -> da = ((p -> a) - (p -> next -> a))/BoxAtd;
      if ((p -> dn <= 0)||(p -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the TD function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                p -> n, "Exiting...");
         exit(1);
      }
      p -> SN = (p -> da)/(p -> dn);
   }
   bgtd.tail -> dn = 0;
   bgtd.tail -> da = 0;
   bgtd.tail -> SN = 0;
   SNmax = 0;
   for (p = bgtd.head; p -> next != NULL; p = p -> next) {
      if (p -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the TD function are",
                "not increasing as the function is",
                "tightened, at N =",
                p -> n, "Exiting...");
         exit(2);
      }
      SNmax = p -> SN;
   }

   /*
    *  Check that the Km2 range-tail kinematic function is
    *  constructed properly.
    */
   for (p = bgkm2t.head; p -> next != NULL; p = p -> next) {
      p -> dn = ((p -> n) - (p -> next -> n))/BoxNkm2t;
      p -> da = ((p -> a) - (p -> next -> a))/BoxAkm2t;
      if ((p -> dn <= 0)||(p -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the Km2 range-tail kinematic function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                p -> n, "Exiting...");
         exit(1);
      }
      p -> SN = (p -> da)/(p -> dn);
   }
   bgkm2t.tail -> dn = 0;
   bgkm2t.tail -> da = 0;
   bgkm2t.tail -> SN = 0;
   SNmax = 0;
   for (p = bgkm2t.head; p -> next != NULL; p = p -> next) {
      if (p -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the km2 range-tail kinematic function are",
                "not increasing as the function is",
                "tightened, at N =",
                p -> n, "Exiting...");
         exit(2);
      }
      SNmax = p -> SN;
   }

   /*
    *  Check that the muon-band kinematic function is
    *  constructed properly.
    */
   for (p = bgkm2b.head; p -> next != NULL; p = p -> next) {
      p -> dn = ((p -> n) - (p -> next -> n))/BoxNkm2b;
      p -> da = ((p -> a) - (p -> next -> a))/BoxAkm2b;
      if ((p -> dn <= 0)||(p -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the muon-band kinematic function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                p -> n, "Exiting...");
         exit(1);
      }
      p -> SN = (p -> da)/(p -> dn);
   }
   bgkm2b.tail -> dn = 0;
   bgkm2b.tail -> da = 0;
   bgkm2b.tail -> SN = 0;
   SNmax = 0;
   for (p = bgkm2b.head; p -> next != NULL; p = p -> next) {
      if (p -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the muon-band kinematic function are",
                "not increasing as the function is",
                "tightened, at N =",
                p -> n, "Exiting...");
         exit(2);
      }
      SNmax = p -> SN;
   }

   /*
    *  Check that the single-beam function is
    *  constructed properly.
    */
   for (b = bgbm1.head; b -> next != NULL; b = b -> next) {
      b -> dn = ((b -> n) - (b -> next -> n))/BoxNbm1;
      b -> da = ((b -> a) - (b -> next -> a))/BoxAbm1;
      if ((b -> dn <= 0)||(b -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the single-beam function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                b -> n, "Exiting...");
         exit(1);
      }
      b -> SN = (b -> da)/(b -> dn);
   }
   bgbm1.tail -> dn = 0;
   bgbm1.tail -> da = 0;
   bgbm1.tail -> SN = 0;
   SNmax = 0;
   for (b = bgbm1.head; b -> next != NULL; b = b -> next) {
      if (b -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the single-beam function are",
                "not increasing as the function is",
                "tightened, at N =",
                b -> n, "Exiting...");
         exit(2);
      }
      SNmax = b -> SN;
   }

   /*
    *  Check that the double-beam function is
    *  constructed properly.
    */
   for (b = bgbm2.head; b -> next != NULL; b = b -> next) {
      b -> dn = ((b -> n) - (b -> next -> n))/BoxNbm2;
      b -> da = ((b -> a) - (b -> next -> a))/BoxAbm2;
      if ((b -> dn <= 0)||(b -> da <= 0)) {
         printf(" %s %s \n %s %f %s \n",
                "N,A values for the double-beam function are",
                "non-optimal, non-unique, and/or not",
                "in order of decreasing N,A, at N =",
                b -> n, "Exiting...");
         exit(1);
      }
      b -> SN = (b -> da)/(b -> dn);
   }
   bgbm2.tail -> dn = 0;
   bgbm2.tail -> da = 0;
   bgbm2.tail -> SN = 0;
   SNmax = 0;
   for (b = bgbm2.head; b -> next != NULL; b = b -> next) {
      if (b -> SN <= SNmax) {
         printf(" %s %s \n %s %f %s \n",
                "S:N ratios for the double-beam function are",
                "not increasing as the function is",
                "tightened, at N =",
                b -> n, "Exiting...");
         exit(2);
      }
      SNmax = b -> SN;
   }


   /*
    *  Calculate the signal:noise ratio in each 7-dimensional cell
    *  as defined by differential changes in the 7 background functions.
    */
   printf("Calculating signal:noise ratio in each 7-dimensional cell... \n");

   first = 1;

   /*
    *  We loop over all possible combinations of PV, TD, BM1, and BM2
    *  function values, because we assume that the rejections and
    *  acceptances of these 4 cuts are all completely independent
    *  (which is not strictly true, but is hopefully adequate).
    */
   for (p1 = bgpv.head; p1 -> next != NULL; p1 = p1 -> next) {
      for (p3 = bgtd.head; p3 -> next != NULL; p3 = p3 -> next) {
         for (b1 = bgbm1.head; b1 -> next != NULL; b1 = b1 -> next) {
            for (b2 = bgbm2.head; b2 -> next != NULL; b2 = b2 -> next) {

	       /*
		*  We don't loop over all possible combinations of
		*  Kp2, Km2 range-tail, and muon-band kinematic
		*  function values, because these 3 kinematic functions
		*  all cut on the same variables (one or more of range,
		*  energy, and momentum).  Therefore, their rejections
		*  and acceptances are NOT independent.  Here we assume
		*  that
		*  (1) In the cell defined by the tightest Kp2, Km2
		*  range-tail, and muon-band kinematic cuts, the
		*  total kinematic acceptance is additive.  That is,
		*  these 3 cuts remove mutually exclusive regions of
		*  the K -> pi nu nubar phase space.
		*  (2) The Kp2 background which passes the tightest
		*  Kp2 kinematic cut position is distributed like
		*  signal in ptot and chirm.  That is, it is ~flat in
		*  ptot, and falling in chirm as chirm increases away
		*  from 0.  So the Kp2 background which gets into the
		*  Km2 range tail and muon band is distributed like dA
		*  in these 2 kinematic functions.
		*  (3) The Km2 range-tail background which passes the
		*  tightest Km2 range-tail kinematic cut position is
		*  distributed like signal in ptot,etot,rtot and chirm.
                *  That is, it is ~flat in ptot,etot,rtot, and falling
                *  in chirm as chirm increases away from 0.  So the Km2
                *  range-tail background which gets into the Kp2 and
                *  muon-band regions is distributed like dA in these 2
                *  kinematic functions.
		*  (4) The muon-band background which passes the
		*  tightest muon-band kinematic cut position is
		*  distributed like signal in ptot, etot, and rtot.
		*  That is, it is ~flat in ptot, etot, and rtot.
		*  So the muon-band background which gets into the Kp2
                *  and range-tail regions is distributed like dA in
                *  these 2 kinematic functions.
		*
		*  First compute dN and dA for the tightest Kp2,
		*  Km2 range-tail, and muon-band kinematic region.
		*/
               p2 = bgkp2.tail -> prev;
               p4 = bgkm2t.tail -> prev;
               p5 = bgkm2b.tail -> prev;

               if (first) {
                  first = 0;
                  t = malloc(sizeof(TOTVALUE));
                  t -> prev = NULL;
                  bgtotal.head = t;
               } else {
                  t -> next = malloc(sizeof(TOTVALUE));
                  t -> next -> prev = t;
                  t = t -> next;
               }

               t1 = t;

               dBgkp2 = (p1 -> dn)*(p2 -> dn)*
                        (p3 -> da)*((p4 -> da) + (p5 -> da) - 1.)*
                        (b1 -> da)*(b2 -> da)*
                        BoxBgkp2;
               dBgkm2t = (p3 -> dn)*(p4 -> dn)*
                         (p1 -> da)*((p2 -> da) + (p5 -> da) - 1.)*
                         (b1 -> da)*(b2 -> da)*
                         BoxBgkm2t;
               dBgkm2b = (p3 -> dn)*(p5 -> dn)*
                         (p1 -> da)*((p2 -> da) + (p4 -> da) - 1.)*
                         (b1 -> da)*(b2 -> da)*
                         BoxBgkm2b;
               kinacc = ((p2 -> da) + (p4 -> da) + (p5 -> da) - 2.);
               if (kinacc < 0.) kinacc = 0.;
               dBgbm1 = (b1 -> dn)*
                        (p1 -> da)*(p3 -> da)*
                        kinacc*(b2 -> da)*
                        BoxBgbm1;
               dBgbm2 = (b2 -> dn)*
                        (p1 -> da)*(p3 -> da)*
                        kinacc*(b1 -> da)*
                        BoxBgbm2;

               dNtotal = (dBgkp2 + dBgkm2t + dBgkm2b +
                          dBgbm1 + dBgbm2)/
                         (BoxBgkp2 + BoxBgkm2t + BoxBgkm2b +
                          BoxBgbm1 + BoxBgbm2);
               dAtotal = (p1 -> da)*(p3 -> da)*kinacc*
                         (b1 -> da)*(b2 -> da);

               t -> dn = dNtotal;
               t -> da = dAtotal;
               t -> SN = dAtotal/dNtotal;
               t -> xord = dAtotal/dNtotal;
               t -> dkp2 = dBgkp2;
               t -> dkm2t = dBgkm2t;
               t -> dkm2b = dBgkm2b;
               t -> dbm1 = dBgbm1;
               t -> dbm2 = dBgbm2;
               t -> ptrpv = p1;
               t -> ptrkp2 = p2;
               t -> ptrtd = p3;
               t -> ptrkm2t = p4;
               t -> ptrkm2b = p5;
               t -> ptrbm1 = b1;
               t -> ptrbm2 = b2;

               t -> next = NULL;
               bgtotal.tail = t;
               bgtotal.nvalues += 1;

	       /*
		*  Now compute dN and dA for the remaining Kp2 regions.
		*  Km2 range-tail and muon-band background are present
                *  proportional to dA in the Kp2 kinematic function.
		*  The kinematic acceptance is just the Kp2 kinematic
		*  acceptance, because these regions are sub-regions
		*  of the tightest Km2 range-tail and muon-band region.
		*/
               p4 = bgkm2t.tail -> prev;
               p5 = bgkm2b.tail -> prev;
               for (p2 = bgkp2.head; p2 -> next -> next != NULL;
                    p2 = p2 -> next) {
                  t -> next = malloc(sizeof(TOTVALUE));
                  t -> next -> prev = t;
                  t = t -> next;

                  if (p2 == bgkp2.tail -> prev -> prev) {
                     t2a = t;
		  }

                  dBgkp2 = (p1 -> dn)*(p2 -> dn)*
                           (p3 -> da)*
                           (b1 -> da)*(b2 -> da)*
                           BoxBgkp2;
                  dBgkm2t = (p3 -> dn)*(p4 -> dn)*
                            (p1 -> da)*(p2 -> da)*
                            (b1 -> da)*(b2 -> da)*
                            BoxBgkm2t;
                  dBgkm2b = (p3 -> dn)*(p5 -> dn)*
                            (p1 -> da)*(p2 -> da)*
                            (b1 -> da)*(b2 -> da)*
                            BoxBgkm2b;
                  kinacc = (p2 -> da);
                  dBgbm1 = (b1 -> dn)*
                           (p1 -> da)*(p3 -> da)*
                           kinacc*(b2 -> da)*
                           BoxBgbm1;
                  dBgbm2 = (b2 -> dn)*
                           (p1 -> da)*(p3 -> da)*
                           kinacc*(b1 -> da)*
                           BoxBgbm2;

                  dNtotal = (dBgkp2 + dBgkm2t + dBgkm2b +
                             dBgbm1 + dBgbm2)/
                            (BoxBgkp2 + BoxBgkm2t + BoxBgkm2b +
                             BoxBgbm1 + BoxBgbm2);
                  dAtotal = (p1 -> da)*(p3 -> da)*kinacc*
                            (b1 -> da)*(b2 -> da);

                  t -> dn = dNtotal;
                  t -> da = dAtotal;
                  t -> SN = dAtotal/dNtotal;
                  t -> xord = dAtotal/dNtotal;
                  t -> dkp2 = dBgkp2;
                  t -> dkm2t = dBgkm2t;
                  t -> dkm2b = dBgkm2b;
                  t -> dbm1 = dBgbm1;
                  t -> dbm2 = dBgbm2;
                  t -> ptrpv = p1;
                  t -> ptrkp2 = p2;
                  t -> ptrtd = p3;
                  t -> ptrkm2t = p4;
                  t -> ptrkm2b = p5;
                  t -> ptrbm1 = b1;
                  t -> ptrbm2 = b2;

                  t -> next = NULL;
                  bgtotal.tail = t;
                  bgtotal.nvalues += 1;
	       }

	       /*
		*  Now compute dN and dA for the remaining Km2
		*  range-tail regions.
		*  Kp2 and muon-band background are present
                *  proportional to dA in the Km2 range-tail
                *  kinematic function.
		*  The kinematic acceptance is just the Km2
		*  range-tail kinematic acceptance, because these
		*  regions are sub-regions of the tightest Kp2 
		*  and muon-band region.
		*/
               p2 = bgkp2.tail -> prev;
               p5 = bgkm2b.tail -> prev;
               for (p4 = bgkm2t.head; p4 -> next -> next != NULL;
                    p4 = p4 -> next) {
                  t -> next = malloc(sizeof(TOTVALUE));
                  t -> next -> prev = t;
                  t = t -> next;

                  if (p4 == bgkm2t.tail -> prev -> prev) {
                     t2b = t;
		  }

                  dBgkp2 = (p1 -> dn)*(p2 -> dn)*
                           (p3 -> da)*(p4 -> da)*
                           (b1 -> da)*(b2 -> da)*
                           BoxBgkp2;
                  dBgkm2t = (p3 -> dn)*(p4 -> dn)*
                            (p1 -> da)*
                            (b1 -> da)*(b2 -> da)*
                            BoxBgkm2t;
                  dBgkm2b = (p3 -> dn)*(p5 -> dn)*
                            (p1 -> da)*(p4 -> da)*
                            (b1 -> da)*(b2 -> da)*
                            BoxBgkm2b;
                  kinacc = (p4 -> da);
                  dBgbm1 = (b1 -> dn)*
                           (p1 -> da)*(p3 -> da)*
                           kinacc*(b2 -> da)*
                           BoxBgbm1;
                  dBgbm2 = (b2 -> dn)*
                           (p1 -> da)*(p3 -> da)*
                           kinacc*(b1 -> da)*
                           BoxBgbm2;

                  dNtotal = (dBgkp2 + dBgkm2t + dBgkm2b +
                             dBgbm1 + dBgbm2)/
                            (BoxBgkp2 + BoxBgkm2t + BoxBgkm2b +
                             BoxBgbm1 + BoxBgbm2);
                  dAtotal = (p1 -> da)*(p3 -> da)*kinacc*
                            (b1 -> da)*(b2 -> da);

                  t -> dn = dNtotal;
                  t -> da = dAtotal;
                  t -> SN = dAtotal/dNtotal;
                  t -> xord = dAtotal/dNtotal;
                  t -> dkp2 = dBgkp2;
                  t -> dkm2t = dBgkm2t;
                  t -> dkm2b = dBgkm2b;
                  t -> dbm1 = dBgbm1;
                  t -> dbm2 = dBgbm2;
                  t -> ptrpv = p1;
                  t -> ptrkp2 = p2;
                  t -> ptrtd = p3;
                  t -> ptrkm2t = p4;
                  t -> ptrkm2b = p5;
                  t -> ptrbm1 = b1;
                  t -> ptrbm2 = b2;

                  t -> next = NULL;
                  bgtotal.tail = t;
                  bgtotal.nvalues += 1;
	       }

	       /*
		*  Now compute dN and dA for the remaining
		*  muon-band regions.
		*  Kp2 and Km2 range-tail background are present
		*  proportional to dA in the muon-band kinematic
		*  function.
		*  The kinematic acceptance is just the
		*  muon-band kinematic acceptance, because these
		*  regions are sub-regions of the tightest Kp2 
		*  and Km2 range-tail region.
		*/
               p2 = bgkp2.tail -> prev;
               p4 = bgkm2t.tail -> prev;
               for (p5 = bgkm2b.head; p5 -> next -> next != NULL;
                    p5 = p5 -> next) {
                  t -> next = malloc(sizeof(TOTVALUE));
                  t -> next -> prev = t;
                  t = t -> next;

                  if (p5 == bgkm2b.tail -> prev -> prev) {
                     t2c = t;
		  }

                  dBgkp2 = (p1 -> dn)*(p2 -> dn)*
                           (p3 -> da)*(p5 -> da)*
                           (b1 -> da)*(b2 -> da)*
                           BoxBgkp2;
                  dBgkm2t = (p3 -> dn)*(p4 -> dn)*
                            (p1 -> da)*(p5 -> da)*
                            (b1 -> da)*(b2 -> da)*
                            BoxBgkm2t;
                  dBgkm2b = (p3 -> dn)*(p5 -> dn)*
                            (p1 -> da)*
                            (b1 -> da)*(b2 -> da)*
                            BoxBgkm2b;
                  kinacc = (p5 -> da);
                  dBgbm1 = (b1 -> dn)*
                           (p1 -> da)*(p3 -> da)*
                           kinacc*(b2 -> da)*
                           BoxBgbm1;
                  dBgbm2 = (b2 -> dn)*
                           (p1 -> da)*(p3 -> da)*
                           kinacc*(b1 -> da)*
                           BoxBgbm2;

                  dNtotal = (dBgkp2 + dBgkm2t + dBgkm2b +
                             dBgbm1 + dBgbm2)/
                            (BoxBgkp2 + BoxBgkm2t + BoxBgkm2b +
                             BoxBgbm1 + BoxBgbm2);
                  dAtotal = (p1 -> da)*(p3 -> da)*kinacc*
                            (b1 -> da)*(b2 -> da);

                  t -> dn = dNtotal;
                  t -> da = dAtotal;
                  t -> SN = dAtotal/dNtotal;
                  t -> xord = dAtotal/dNtotal;
                  t -> dkp2 = dBgkp2;
                  t -> dkm2t = dBgkm2t;
                  t -> dkm2b = dBgkm2b;
                  t -> dbm1 = dBgbm1;
                  t -> dbm2 = dBgbm2;
                  t -> ptrpv = p1;
                  t -> ptrkp2 = p2;
                  t -> ptrtd = p3;
                  t -> ptrkm2t = p4;
                  t -> ptrkm2b = p5;
                  t -> ptrbm1 = b1;
                  t -> ptrbm2 = b2;

                  t -> next = NULL;
                  bgtotal.tail = t;
                  bgtotal.nvalues += 1;
	       }

	       /*
		*  Due to the simplicity of the model above, it's possible
		*  for the tightest kinematic region to have smaller
		*  signal:noise than the 2nd tightest region (typically by
                *  by a small amount, ~1%).  If this happens, force the
                *  2nd tightest region to have marginally smaller
                *  signal:noise than the tighest region.  The ordering of
                *  signal:noise regions is done below in "quicksort" using
                *  the parameter "xord", so if one desires, one can leave
                *  the signal:noise ratios (SN) as is, and simply change
                *  the values of xord.
		*/
               if ((t1 -> SN < t2a -> SN)||
                   (t1 -> SN < t2b -> SN)||
                   (t1 -> SN < t2c -> SN)) {
                  printf ("Wrong order! %13.10f %13.10f %13.10f %13.10f %f %f %f %f \n",
                          t1 -> SN, t2a -> SN, t2b -> SN, t2c -> SN,
                          t1 -> ptrpv -> n,
                          t1 -> ptrtd -> n,
                          t1 -> ptrbm1 -> n,
                          t1 -> ptrbm2 -> n);
		  /*
                  exit(3);
		  */
                  if (t1 -> SN < t2a -> SN)
                     t2a -> xord = t1 -> SN - 0.0000001;
                     t2a -> SN = t1 -> SN - 0.0000001;
                  if (t1 -> SN < t2b -> SN)
                     t2b -> xord = t1 -> SN - 0.0000001;
                     t2b -> SN = t1 -> SN - 0.0000001;
                  if (t1 -> SN < t2c -> SN)
                     t2c -> xord = t1 -> SN - 0.0000001;
                     t2c -> SN = t1 -> SN - 0.0000001;
	       }
	    }
	 }
      }
   }

   /*
    *  Do a "quick sort" on the order parameter "xord" such that
    *  the list is in decreasing order of signal:noise ratios.
    */
   printf("Ordering the signal:noise ratios (quick sort of %d values) \n",
          bgtotal.nvalues);
   quicksort(bgtotal.head, bgtotal.tail);

   /*
    *  Calculate the integrated dn and da values (i.e., N and A values)
    *  for each cut position, and print to a file.
    */
   printf("Calculating total N and A values... \n");
   ofp = fopen("bgtotal.dat","w");
   fprintf(ofp,"%d \n",bgtotal.nvalues);

   Nlo = 0.;
   Alo = 0.;
   for (t = bgtotal.head; t != NULL; t = t -> next) {
      if (t == bgtotal.head) {
         t -> n = Nlo + t -> dn;
         t -> a = Alo + t -> da;
      } else {
         t -> n = (t -> prev -> n) + (t -> dn);
         t -> a = (t -> prev -> a) + (t -> da);
      }
      fprintf(ofp,"%f %f %14.11f %14.11f %f %f %f %f %f %f %f %f \n",
              t -> n, t -> a, t -> dn, t -> da, t -> SN,
              t -> ptrpv -> n,
              t -> ptrkp2 -> n,
              t -> ptrtd -> n,
              t -> ptrkm2t -> n,
              t -> ptrkm2b -> n,
              t -> ptrbm1 -> n,
              t -> ptrbm2 -> n);
   }
   fclose(ofp);

   /*
    *  Calculate the integrated dn and da values (i.e., N and A values)
    *  for each cut position and print to a file, this time for
    *  inside-the-box and outside-the-box cells separately.
    */
   ofpin = fopen("bgtotal.inside.dat","w");
   ofpout = fopen("bgtotal.outside.dat","w");

   ninside = 0;
   noutside = 0;
   for (t = bgtotal.head; t != NULL; t = t -> next) {
      if ((t -> ptrpv -> n <= bgpv.box -> n)&&
          (t -> ptrkp2 -> n <= bgkp2.box -> n)&&
          (t -> ptrtd -> n <= bgtd.box -> n)&&
          (t -> ptrkm2t -> n <= bgkm2t.box -> n)&&
          (t -> ptrkm2b -> n <= bgkm2b.box -> n)&&
          (t -> ptrbm1 -> n <= bgbm1.box -> n)&&
          (t -> ptrbm2 -> n <= bgbm2.box -> n)) {
         ninside += 1;
      } else {
         noutside += 1;
      }
   }

   fprintf(ofpin,"%d \n",ninside);
   fprintf(ofpout,"%d \n",noutside);

   Ninlo = 0.;
   Ainlo = 0.;
   Noutlo = 1.;
   Aoutlo = 1.;
   firstin = 1;
   firstout = 1;
   for (t = bgtotal.head; t != NULL; t = t -> next) {
      if ((t -> ptrpv -> n <= bgpv.box -> n)&&
          (t -> ptrkp2 -> n <= bgkp2.box -> n)&&
          (t -> ptrtd -> n <= bgtd.box -> n)&&
          (t -> ptrkm2t -> n <= bgkm2t.box -> n)&&
          (t -> ptrkm2b -> n <= bgkm2b.box -> n)&&
          (t -> ptrbm1 -> n <= bgbm1.box -> n)&&
          (t -> ptrbm2 -> n <= bgbm2.box -> n)) {
         if (firstin) {
            firstin = 0;
            Nin = Ninlo + t -> dn;
            Ain = Ainlo + t -> da;
         } else {
            Nin += (t -> dn);
            Ain += (t -> da);
         }
         fprintf(ofpin,"%f %f %f %f %f %f %f %f %f %f %f %f \n",
                 Nin, Ain, t -> dn, t -> da, t -> SN,
                 t -> ptrpv -> n,
                 t -> ptrkp2 -> n,
                 t -> ptrtd -> n,
                 t -> ptrkm2t -> n,
                 t -> ptrkm2b -> n,
                 t -> ptrbm1 -> n,
                 t -> ptrbm2 -> n);
      } else {
         if (firstout) {
            firstout = 0;
            Nout = Noutlo + t -> dn;
            Aout = Aoutlo + t -> da;
         } else {
            Nout += (t -> dn);
            Aout += (t -> da);
         }
         fprintf(ofpout,"%f %f %f %f %f %f %f %f %f %f %f %f \n",
                 Nout, Aout, t -> dn, t -> da, t -> SN,
                 t -> ptrpv -> n,
                 t -> ptrkp2 -> n,
                 t -> ptrtd -> n,
                 t -> ptrkm2t -> n,
                 t -> ptrkm2b -> n,
                 t -> ptrbm1 -> n,
                 t -> ptrbm2 -> n);
      }
   }
   fclose(ofpin);
   fclose(ofpout);

   printf("Done. \n");

   /*
    *  Free the used memory.
    */
   for (p = bgpv.head; p != NULL; p = p -> next) free(p);
   for (p = bgkp2.head; p != NULL; p = p -> next) free(p);
   for (p = bgtd.head; p != NULL; p = p -> next) free(p);
   for (p = bgkm2t.head; p != NULL; p = p -> next) free(p);
   for (p = bgkm2b.head; p != NULL; p = p -> next) free(p);
   for (b = bgbm1.head; b != NULL; b = b -> next) free(b);
   for (b = bgbm2.head; b != NULL; b = b -> next) free(b);

   exit(0);
}
