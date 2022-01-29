c-----------------------------------------------------------
c     Application-specific common block for communication 
c     with derivative evaluating routine 'fcn'.
c
c      w:    Frequency eigenvalue (sought eigenvalue)
c      A, B: Parameters defining the potential (arbitrary)
c
c-----------------------------------------------------------

      real*8   w, A, B
      
      common / com_fcn / w, A, B 
