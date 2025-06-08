      subroutine constr(  NGSK,   neq,  neqm, NGSKo,  neqo,
     &                    mdof,  idof, jdiag,index0,index1,
     &                      df,
     &                     res,  disp,  del,    sk,
     &                  ierror )
c
      implicit double precision (a-h,o-z)
c
      dimension mdof(neq),idof(neq)
      dimension jdiag(neq+1)
      dimension index0(neqo+1)
      dimension index1(NGSKo)
c
      dimension res(neq),disp(neq),del(neq)
      dimension sk(NGSK)
c
      common /basic/ nx,nelx,ndf,node,nsn,lmat,lang,ngaus
c
c **********************************************************************
      NGSKm = NGSK -NGSKo
      del = 0.d0 ! del(:) = 0.d0
c
c ***** Re-order the Total Force Vector ********************************
      do ne=1,neq
        m_ne = mdof(ne)
        del(m_ne) = res(ne)
      enddo
c
c     do ne=1,neq
c       write(*,'( i5, e15.5)') ne,mdof(ne)
c       write(*,'(     e15.5)')             res(ne)
c       write(*,'(     e15.5)')                     del(ne)
c       write(*,'(     e15.5)')                             disp(ne)
c     enddo
c     do ne=1,neq
c       write(*,'(2i5,3e15.5)') ne,mdof(ne),res(ne),del(ne),disp(ne)
c     enddo
c
c ***** Convert Enforced Displacement to Equivalent Force **************
c     do ia=1,neqm
c       fdisp = 0.d0
c       do jb=neqm+1,neq
c         i_jb = idof(jb)
c         kij = jb -ia
c         nd = jdiag(jb) -kij
c         if(nd.gt.jdiag(jb-1)) then
c           fdisp = fdisp +sk(nd)*disp(i_jb)
c           write(*,*) i_jb,disp(i_jb),fdisp
c         endif
c       enddo
c       write(*,*) ia,fdisp
c       del(ia) = del(ia) -fdisp
c     enddo
c ***** Convert Enforced Displacement to Equivalent Force **************
c         ( Only for Non-Zero Terms )
      do neo=neqm+1,neq
        ido = neo -neqm
        inds = index0(ido)
        inde = index0(ido+1) -2
        do ind=inds,inde
          m_ia = index1(ind)
          i_jb = idof(neo)
          nd = NGSKm +ind
c         del(m_ia) = del(m_ia) -sk(nd)*disp(i_jb)
          del(m_ia) = del(m_ia) -sk(nd)*df*disp(i_jb)
c         write(*,'(10i5)') neo, ido, ind, m_ia, i_jb, nd
        enddo
      enddo
c
c ***** Set Enforced Displacement and Potential ************************
      do ne=neqm+1,neq
        id = idof(ne)
c       del(ne) = disp(id)
        del(ne) = df*disp(id)
      enddo
c
c     write(*,*) sk
c     do ne=1,neq
c       write(*,'( i5, e15.5)') ne,mdof(ne)
c       write(*,'(     e15.5)')             res(ne)
c       write(*,'(     e15.5)')                     del(ne)
c       write(*,'(     e15.5)')                             disp(ne)
c     enddo
c     do ne=1,neq
c       write(*,'(2i5,3e15.5)') ne,mdof(ne),foc(ne),del(ne),disp(ne)
c     enddo
c **********************************************************************
c **********************************************************************
      RETURN
      END
