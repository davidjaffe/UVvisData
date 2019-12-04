c common storage for combined probability distributions.
c  -- necessary so that the average confidence levels can be
c  computed after the CL's without needing a recomputation of
c  the pdf for the estimator.

	integer naccsb,naccb
	real*8 sbaccp(maxacc),sbaccest(maxacc)
        real*8 baccp(maxacc),baccest(maxacc)
	common /eclplist/ naccsb,naccb,
     >                    sbaccp,sbaccest,baccp,baccest
