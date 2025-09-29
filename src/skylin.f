      subroutine skylin(     a,     b, jdiag,   nes,   neq,   kkk )
c
c ---------------------------------------------------------------------
c (c) copyright 2002 by K. Terada & K. Matsui, Tohoku University
c ---------------------------------------------------------------------
c
c   Purpose:
c       solve the symmetric system of linear equations ( r.l.taylor )
c
c   Arguments:
c       a     = coeffecient matrix
c       b     = load vector at the time of input, and solution at the
c               end of the routine
c       jdiag = diagonal pointers
c       nes   = first equation number to be solved
c       neq   = last equation number to be solved
c       kkk   = 0  ,  if the solver of the linear system is assumed
c
c ----------------------------------------------------------------------
c
      implicit double precision(a-h,o-z)
c
      dimension a(1),b(1),jdiag(1)
c
c     do ne=1,neq
c       write(*,*) b(ne)
c     enddo
c
      if(nes.gt.neq) RETURN
      jr=0
      if(nes.gt.1) jr=jdiag(nes-1)
c
c     ( factor the matrix "a" to "ut*d*u" and reduce the vector "b" )
c
      do 1600 j=nes,neq
        jd = jdiag(j)
        jh = jd -jr
        il = j -jh +2
        if(jh-2) 600,300,100
c ++++++++++++++++++++++++++
  100   continue
          if(kkk.eq.2) go to 500
          ie = j -1
          k = jr +2
          id = jdiag(il -1)
c
c     < reduce all equations except diagonal >
c
          do 200 i=il,ie
            ir = id
            id = jdiag(i)
            ih = min0(id-ir-1,i-il+1)
            if(ih.gt.0) a(k)=a(k)-dot(a(k-ih),a(id-ih),ih)
            k = k +1
  200     continue
c
c     < reduce the diagnal >
c ++++++++++++++++++++++++++
  300     continue
          if(kkk.eq.2) go to 500
          ir = jr +1
          ie = jd -1
          k = j -jd
          do 400 i=ir,ie
            id = jdiag(k+i)
            if(a(id).eq.0.d0) go to 400
            d = a(i)
            a(i) = a(i)/a(id)
            a(jd) = a(jd)-d*a(i)
  400     continue
c
c     < reduce the load vector >
c
  500     if(kkk.ne.1) b(j)=b(j)-dot(a(jr+1),b(il-1),jh-1)
c ++++++++++++++++++++++++++
  600   continue
 1600 jr = jd
      if(kkk.eq.1) RETURN
c
c     ( divided by the diagonal pivots )
c
      do 700 i=1,neq                       
        id = jdiag(i)
        if(a(id).ne.0.d0) b(i)=b(i)/a(id)
  700 continue
c
c     ( back substitution )
c
 1100 j = neq
      jd = jdiag(j)
  800 d = b(j)
      j = j -1
      if(j.le.0) RETURN
      jr = jdiag(j)
      if(jd-jr.le.1) go to 1000
      il = j -jd +jr +2
      k = jr -il +1
      do 900 i=il,j
        b(i)=b(i)-a(i+k)*d
  900 continue
 1000 jd = jr
      go to 800
c
      END
c
c ----------------------------------------------------------------------
c
      function dot(a,b,n)
c
c   purpose:
c       take the dot product of two vectors  "a" and "b" .
c
c ----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
      dimension a(1),b(1)
c
      dot=0.
      do 100 i=1,n
  100 dot = dot +a(i)*b(i)
c
      RETURN
      END
