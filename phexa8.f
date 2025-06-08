      subroutine phexa8(   nel,  nelx,   ndf,  node, ngaus,
     &                    idep,MATYPE, prope,    xe,    ue,
     &                     sts,   stn, finte,  vone,  epse,
     &                  dhist0,dhist1,
     &                   e_ene, p_ene,tene_e,dene_p, tempe,
     &                   dtemp,dtempe,
     &                  dndx_g, det_g,ctensg,
     &                  ierror, itr, histi0 )
c
      implicit double precision (a-h,o-z)
c
      dimension dndx_g(ndf,node,ngaus,nelx)
      dimension det_g(ngaus,nelx)
      dimension ctensg(3,3,3,3,ngaus,nelx)
      dimension dhist0(20,ngaus,nelx),dhist1(20,ngaus,nelx)
      dimension histi0(50,ngaus,nelx)
c
      dimension idep(8)
c
      dimension prope(20)
      dimension xe(3,8),dn(3,8),dndx(3,8),ue(3,8)
      dimension wg(2),xg(2)
      dimension dxdl(3,3),dldx(3,3),
     &          sig(3,3),sts(3,3),str(3,3),stn(3,3),
     &          sd(3,3),stry(3,3),oun(3,3),seta(3,3)
      dimension ctens(3,3,3,3)
      dimension finte(24)
      dimension ehist(20)
      dimension histi(50)
c
      common /tvalu/ ctol,stol
c **********************************************************************
      itrmax = 20
c
      wg(1) = 1.d0
      wg(2) = 1.d0
      xg(1) = -1.d0/dsqrt(3.d0)
      xg(2) =  1.d0/dsqrt(3.d0)
c
c ***** Initilization **************************************************
      finte = 0.d0 ! finte(:) = 0.d0
c
      area = 0.d0
      vone = 0.d0
      epse = 0.d0
      p_ene = 0.d0
      e_ene = 0.d0
c
      sts = 0.d0 ! sts(:,:) = 0.d0
      stn = 0.d0 ! stn(:,:) = 0.d0
c
c ***** Set Material Properties ****************************************
c    --- thermal parameters
c     row = prope(3)
c     ccc = prope(4)
c     aaa = prope(5)
c
c ***** Start Numerical Integration ************************************
      ig = 0
      do 100 ix=1,2
        xl = xg(ix)
        wx = wg(ix)
        do 110 iy=1,2
          yl = xg(iy)
          wy = wg(iy)
          do 120 iz=1,2
            zl = xg(iz)
            wz = wg(iz)
            ig = ig +1
c
c         === Initilization of Local Tensors ===
            sig = 0.d0 ! sig(:,:) = 0.d0
            str = 0.d0 ! str(:,:) = 0.d0
c
c         === Derivative of Shape Function r.w.t. Local Coordinate ===
c             ( dn_ij = d(N^j)/d(l_i) )
c           dn(1,1) = -0.125d0*(1.d0 -yl)*(1.d0 -zl)
c           dn(2,1) = -0.125d0*(1.d0 -xl)*(1.d0 -zl)
c           dn(3,1) = -0.125d0*(1.d0 -xl)*(1.d0 -yl)
c           dn(1,2) =  0.125d0*(1.d0 -yl)*(1.d0 -zl)
c           dn(2,2) = -0.125d0*(1.d0 +xl)*(1.d0 -zl)
c           dn(3,2) = -0.125d0*(1.d0 +xl)*(1.d0 -yl)
c           dn(1,3) =  0.125d0*(1.d0 +yl)*(1.d0 -zl)
c           dn(2,3) =  0.125d0*(1.d0 +xl)*(1.d0 -zl)
c           dn(3,3) = -0.125d0*(1.d0 +xl)*(1.d0 +yl)
c           dn(1,4) = -0.125d0*(1.d0 +yl)*(1.d0 -zl)
c           dn(2,4) =  0.125d0*(1.d0 -xl)*(1.d0 -zl)
c           dn(3,4) = -0.125d0*(1.d0 -xl)*(1.d0 +yl)
c           dn(1,5) = -0.125d0*(1.d0 -yl)*(1.d0 +zl)
c           dn(2,5) = -0.125d0*(1.d0 -xl)*(1.d0 +zl)
c           dn(3,5) =  0.125d0*(1.d0 -xl)*(1.d0 -yl)
c           dn(1,6) =  0.125d0*(1.d0 -yl)*(1.d0 +zl)
c           dn(2,6) = -0.125d0*(1.d0 +xl)*(1.d0 +zl)
c           dn(3,6) =  0.125d0*(1.d0 +xl)*(1.d0 -yl)
c           dn(1,7) =  0.125d0*(1.d0 +yl)*(1.d0 +zl)
c           dn(2,7) =  0.125d0*(1.d0 +xl)*(1.d0 +zl)
c           dn(3,7) =  0.125d0*(1.d0 +xl)*(1.d0 +yl)
c           dn(1,8) = -0.125d0*(1.d0 +yl)*(1.d0 +zl)
c           dn(2,8) =  0.125d0*(1.d0 -xl)*(1.d0 +zl)
c           dn(3,8) =  0.125d0*(1.d0 -xl)*(1.d0 +yl)
c
c       === Derivative of Current Position r.w.t. Loacal Coordinate ===
c           ( dxdl_ij = d(x_i)/d(l_j) )
c           do ndia=1,ndf
c             do ndjb=1,ndf
c               dxdlj = 0.d0
c               do ia=1,node
c                 dxdlj = dxdlj +xe(ndia,ia)*dn(ndjb,ia)
c               enddo
c               dxdl(ndia,ndjb) = dxdlj
c             enddo
c           enddo
c
c         === Derivative of Local Coordinate r.w.t. Current Position ===
c             ( dldx_ij = d(li)/d(xj) = ( d(xi)/d(lj) )^{-1} )
c           CALL inv_33(  dxdl,   det,  dldx,ierror )
c           if(det.le.0.d0) CALL zerodt(   nel,    ig,   det)
c
            det = det_g(ig,nel)
c           det_g(ig,nel) = det
            detwxyz = det*wx*wy*wz
            area = area +detwxyz
c
c       === Derivative of Local Coordinate r.w.t. Current Position ===
c           ( dldx_ij = (dxdl_ij)^{-1} )
c           deti = 1.d0/(det*2.d0)
c           do jj=1,3
c             do ii=1,3
c               dldx(ii,jj) = dldx(ii,jj)*deti
c             enddo
c           enddo
c
c       === Derivative of Shape Function r.w.t. Current Position ===
c           ( dndx = d(N^i)/d(x_j) = d(N^i)/d(l_k) d(l_k)/d(x_j) )
            dndx(:,:) = dndx_g(:,:,ig,nel)
c           do ndia=1,ndf
c             do ia=1,node
c               dndxj = 0.d0
c               do ndjb=1,ndf
c                 dndxj = dndxj +dn(ndjb,ia)*dldx(ndjb,ndia)
c               enddo
c               dndx(ndia,ia) = dndxj
c               dndx_g(ndia,ia,ig,nel) = dndxj
c             enddo
c           enddo
c
c       === Compute the Local Strain Tensor ===
            do jj=1,ndf
              do ii=1,ndf
                eij = 0.d0
                do no=1,node
                  eij = eij +0.5d0*(
     &                  ue(ii,no)*dndx(jj,no) +ue(jj,no)*dndx(ii,no) )
                enddo
                str(ii,jj) = eij
              enddo
            enddo
c
c         === Set Deformation Histtories ===
            ehist(:) = dhist0(:,ig,nel)
            histi(:) = histi0(:,ig,nel)
c
c         === Compute Local Stresses ===
c            print'("--------------------------------ig=",2I5)',ig,MATYPE
          if( MATYPE.eq.1) then
            CALL stress(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror, itr, histi )
          elseif(MATYPE.eq.2) then
            CALL st_gtn(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror, itr, histi )
          elseif(MATYPE.eq.3) then
            CALL stress_vm(itrmax, idepg,
     &                   prope,   sig,   str, ehist,
     &                    ctol,  vons, e_dns, p_dns,
     &                   ctens,
     &                  ierror, itr, histi )
          else
            STOP 'Something wrong in phexa8'
          endif
            if(ierror.ne.0) RETURN
c
c ***** Store Deformation Histories (at current state) *****************
            dhist1(:,ig,nel) = ehist(:)
            histi0(:,ig,nel) = histi(:)
            ctensg(:,:,:,:,ig,nel) = ctens(:,:,:,:)
            idep(ig) = idepg
c
            alpeg = ehist(1)
c ***** Compute the Elemental Value (Sum-up for Element ) **************
            vone = vone +vons*detwxyz
            epse = epse +alpeg*detwxyz
            sts(:,:) = sts(:,:) +sig(:,:)*detwxyz
            stn(:,:) = stn(:,:) +str(:,:)*detwxyz
            e_ene = e_ene +e_dns*detwxyz
            p_ene = p_ene +p_dns*detwxyz
c
c ***** Compute the Internal Force Vector ******************************
            do no=1,node
              do ii=1,ndf
                ndof = ndf*(no -1) +ii
                qol = 0.d0
                do jj=1,ndf
                  qol = qol +sig(jj,ii)*dndx(jj,no)
                enddo
                finte(ndof) = finte(ndof) +qol*detwxyz
              enddo
            enddo
c
  120     continue
  110   continue
  100 continue
c
c ***** Compute the Elemental Valume (Volume Average) ******************
      vone = vone/area
      epse = epse/area
c
      tene_e = tene_e +e_ene
      dene_p = dene_p +p_ene
c
      e_ene = e_ene/area
      p_ene = p_ene/area
c
      sts = sts/area ! sts(:,:) = sts(:,:)/area
      stn = stn/area ! sts(:,:) = stn(:,:)/area
c
c ***** Convert the Plastic ene. into Temperature *****
c
c     bbb = aaa*row*ccc
c     if( bbb.eq.0.d0) then
c       dtempe = 0.d0
c     else
c       dtempe = (p_ene*10.d0**9)/(aaa*row*ccc)
c     endif
c     dtemp = dtemp + dtempe*area
c
c **********************************************************************
c **********************************************************************
      RETURN
      END
