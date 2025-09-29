      subroutine restor(   neq,
     &                      u0,    u1,dhist0,dhist1,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension u0(neq),u1(neq)
      dimension dhist0(20*ngaus*nelx),dhist1(20*ngaus*nelx)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
c **********************************************************************
c
c ***** Restore Displacements & etc... *********************************
      u1 = u0 ! u1(:) = u0(:)
c
      dhist1 = dhist0 ! dhist1(:) = dhist0(:)
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
