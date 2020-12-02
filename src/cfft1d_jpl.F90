      subroutine cfft1d_jpl(n,c,dir)

      implicit none

      integer*4  n
      integer dir
      complex*8  c,cout
      dimension c(n)
      dimension cout(n)


#if defined(FFTW) || defined(HAVE_FFTW)
!     NOTE: if above condition changed, also need to update fftw3stub.c

#include <fftw3.f>

      integer*8 plan

      if(dir.eq.-1)then
          call sfftw_plan_dft_1d(plan,n,c,cout,FFTW_FORWARD,FFTW_ESTIMATE)
      end if
      if(dir.eq. 1)then
          call sfftw_plan_dft_1d(plan,n,c,cout,FFTW_BACKWARD,FFTW_ESTIMATE)
      end if
      call sfftw_execute_dft(plan,c,cout)
      call sfftw_destroy_plan(plan)

      c=cout

#else
!     NO FFT routine has been specified
!     force compilation to fail, with below "ABORT" syntax error
!     rather than old behavior, of having this routine
!     silently do nothing
!     ABORT NO FFT ROUTINE DEFINED
      stop  "NO FFT ROUTINE DEFINED"

#endif

      return
      end
