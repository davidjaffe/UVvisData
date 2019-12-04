/*
 *    File:  bgtotal.h
 * Purpose:  include file for bgtotal.c
 *   Usage:  #include "bgtotal.h"
 *
 * History:  2001 Jul 25  PB  created
 */

/*
 *  Template structure for the PV, TD, Kp2 and Km2 kinematic
 *  background functions.
 */
struct dkvalue {
   float n;
   float a;
   float dn;
   float da;
   float SN;
   struct dkvalue *next;
   struct dkvalue *prev;
};
typedef struct dkvalue DKVALUE;
typedef DKVALUE *PTRDKVALUE;

struct dkfunc {
   int nvalues;
   PTRDKVALUE head,box,tail;
};
typedef struct dkfunc DKFUNC;

/*
 *  Template structure for the beam background functions.
 */
struct bmvalue {
   float n;
   float a;
   float dn;
   float da;
   float SN;
   float kaon;
   float pion;
   float cex;
   struct bmvalue *next;
   struct bmvalue *prev;
};
typedef struct bmvalue BMVALUE;
typedef BMVALUE *PTRBMVALUE;

struct bmfunc {
   int nvalues;
   PTRBMVALUE head,box,tail;
};
typedef struct bmfunc BMFUNC;

/*
 *  Template structure for the total background function.
 */
struct totvalue {
   float n;
   float a;
   float dn;
   float da;
   float SN;
   float xord;
   float dkp2;
   float dkm2t;
   float dkm2b;
   float dbm1;
   float dbm2;
   PTRDKVALUE ptrpv,ptrkp2,ptrtd,ptrkm2t,ptrkm2b;
   PTRBMVALUE ptrbm1,ptrbm2;
   struct totvalue *next;
   struct totvalue *prev;
};
typedef struct totvalue TOTVALUE;
typedef TOTVALUE *PTRTOTVALUE;

struct totfunc {
   int nvalues;
   PTRTOTVALUE head,tail;
};
typedef struct totfunc TOTFUNC;

