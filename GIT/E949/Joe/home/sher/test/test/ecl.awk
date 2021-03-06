#
#    File:  ecl.awk
# Purpose:  make ecl*.vect files from ecl*.dat files, for use by ecl.kumac
#           and eclresults.c
#   Usage:  > awk -f ecl.awk ecl{label}.dat >& ecl{label}.vect
#
# History:  2001 Sep  PB  created
#

BEGIN {
    i = 0;

# 1998 single-event-sensitivity in units of E-10
    sens = 1.892;
}

/signal scaling:/ {
    scale[i] = $3;
    br[i] = $3*sens;
    if (br[i] < 0.001) {
	br[i] = 0.001;
    }
}

/Xobs:/ {
    xobs[i] = $2;
}

/CL\(s\):/ {
    cls[i] = $2;

    if (cls[i] < 0.5) {
       cls2[i] = 1 - cls[i];
    } else {
       cls2[i] = cls[i];
    }

    if (i > 0) {
       dy = -(cls[i] - cls[i-1]);
       dx = br[i] - br[i-1];
       dydx[i-1] = dy/dx;
    }
}

/CL1\(b\):/ {
    clb[i] = $2;

    i += 1;
}

END {
# normalize the area under the dCLs/dBR curve to 1.
    nval = i;
    dydysum = 0;
    for (j = 0; j < nval; j++) {
       dydxsum += dydx[j];
    }
    for (j = 0; j < nval; j++) {
       dydx[j] /= dydxsum;
    }

    for (j = 0; j < nval; j++) {
       printf("%10.6f %10.6f %10.6f %10.6f %10.6f \n",
              xobs[j],br[j],cls[j],dydx[j],clb[j]);
    }
}
