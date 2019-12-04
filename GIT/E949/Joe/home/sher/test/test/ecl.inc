c maximum number of separate probabilities to compute for a single channel
c (i.e., max number of events in any one channel for which the probability
c of occurrence is calculated)

	integer maxlist
	parameter (maxlist=100)

c number of probability bins in the linear binning region above 1%.
c a value of 990.0D0 means the bins are spaced by 99%/990 = 0.1%.

	real*8 scaleplin
	parameter(scaleplin=990.0D0)

c maximum number of separate probabilities to carry around for the
c combined channels -- see the subroutine e_probbin for computation
c of the probability bins.  
c It's a function of the fineness of the linear scan region.

        integer maxacc
	parameter (maxacc=200+scaleplin)
c	parameter (maxacc=500+scaleplin)

c don't bother assigning individual probabilities in lists
c if the cumulative probability sum exceeds pmaxcut.

	real*8 pmaxcut
	parameter (pmaxcut=0.99999)

c write some debug output if the inputs do not make sense
	logical check_inputs
	parameter (check_inputs=.true.)

c for evaluation of systematic uncertainty -- integral over Gaussian
c distributions of s and b in number of sigma and the step
c for example, systwid of 3.0, dsystwid of 0.6 steps the s and b
c from -3 sigma to +3 sigma in steps of 0.6 sigma (11 steps)

	real*8 systwid,dsystwid
	parameter (systwid=3.0D0)
	parameter (dsystwid=0.6D0)

c a big estimator for those times when the background has zero
c probability.  For those machines that do not have an inf.
c keep it within a single precision limit so sortzv doesn't have trouble

	real*8 estbig
	parameter(estbig=1.0D25)

c need to declare some internal storage for systematic error fluctuations
c on s and b  -- have a problem with size, upgrade this!
c Also need this to allocate some space for a sort order.

	integer maxbin
	parameter (maxbin=8000)

c prameter for e_simadd -- the width chosen for adding search channel
c results -- discriminator here is log(1+s/b) -- bigger number combines
c more together by plain addition; setting it to zero asks for no
c combination of search results by addition.  (also one does not need
c to call e_simadd)

	real*8 addwid
	parameter (addwid=0.1D0)

c parameter for computation of CL with the test statistic for the data --
c use this many MC evaluations of the data test statistic with different
c values of s and b within their systematic uncertainties

	integer nexpsys
	parameter (nexpsys=200)
