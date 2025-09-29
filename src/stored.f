      subroutine stored(   neq,
     &                  dfact0, dfact,
     &                   idep0, idep1,
     &                      u0,    u1,dhist0,dhist1,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension idep0(ngaus*nelx),idep1(ngaus*nelx)
c
      dimension u0(neq),u1(neq)
      dimension dhist0(20*ngaus*nelx),dhist1(20*ngaus*nelx)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
c **********************************************************************
c
c ***** Update Displacements & etc... **********************************
      u0 = u1 ! u0(:) = u1(:)
c
      dfact0 = dfact
c
c      print *,"before"
c      print *,dhist0
      dhist0 = dhist1 ! dhist0(:) = dhist1(:)
      idep0 = idep1 ! idep1(:) = idep1(:)
c      print *,"after"
c      print *,dhist0
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
