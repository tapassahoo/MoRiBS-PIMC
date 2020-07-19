      subroutine caleng(com_1,com_2,E_2H2O,Eulang_1,Eulang_2)
c     =======================================================
c     This subroutine calculates pair-wise interaction potential 
c     between two rigid waters given the coordinates of their 
c     centres of mass and their respective Euler angles
c
c     input: com_1, com_2, Eulang_1, Eulang_2
c     output: E_2H2O 
c     e: com_1, com_2 --> COM of two water molecules 
c        Eulang_1 and Eulang_2 --> Eular angles 
c        rotmat_1, rotmat_2 --> rotational matrix 
c        ROwf, R1wf, R2wf, RMwf --> Geometry of a water molecule
c     s: E2H2O --> estimated potential in Kelvin
c     _______________________________________________________
c        rotmat_1, rotmat_2 computed within the code with 
c        Eulang_1, Eulang_2, respectively.
c     _______________________________________________________
c     A general purpose model for the condensed phases of water: 
c     q-TIP4P/F --> J. Chem. Phys. 131, 024501 (2009).
c     _______________________________________________________
      implicit double precision(a-h,o-z)
      parameter(zero=0.d0)
      dimension ROwf(3),RH1wf(3),RH2wf(3),RMwf(3),
     +          com_1(3),com_2(3),Eulang_1(3),Eulang_2(3),
     +          RO_1_sf(3), RO_2_sf(3),
     +          RH1_1_sf(3), RH1_2_sf(3),
     +          RH2_1_sf(3), RH2_2_sf(3),
     +          RM_1_sf(3), RM_2_sf(3), 
     +          rotmat_2(3,3), rotmat_1(3,3)
      parameter(br2ang=0.52917721092d0,hr2k=3.1577502480407d5,
     +          akcal2k=503.228d0)
c     q-TIP4P/F parameters
      parameter(angHOH=107.4d0,dOH=0.9419d0,agamma=0.73612d0,
     +          epsoo1=0.1852d0,sigoo=3.1589d0,qh=0.5564d0)
      data ROwf/zero,zero,zero/
c     units
c     angHOH in degree
c     dOH in Angstrom
c     dOM in Angstrom
c     epsoo1 in kcal/mole
c     sigoo in Angstrom
c     qh in e

      pi=4.0d0*datan(1.0d0)
c     Geometry of water molecules
      ang1=(angHOH*pi)/180.0d0
      zH=ROwf(3)-dsqrt(0.5*dOH*dOH*(1.0+dcos(ang1))) 
      xH=dsqrt(dOH*dOH-(ROwf(3)-zH)*(ROwf(3)-zH))
c...
      RH1wf(1)=xH
      RH1wf(2)=zero
      RH1wf(3)=zH
c...
      RH2wf(1)=-RH1wf(1)
      RH2wf(2)=RH1wf(2)
      RH2wf(3)=RH1wf(3)
c...
      do i=1,3
         RMwf(i)=agamma*ROwf(i)+0.5*(1.0-agamma)*(RH1wf(i)+RH2wf(i))
      enddo
c...
      epsoo=epsoo1*akcal2k
      qm=-2.0d0*qh
c     
c     prepare rotational matrix for water 1
c     obtain the SFF coordinates for H1, H2, and O of water 1
      call matpre(Eulang_1, rotmat_1)
      do i=1,3
         RO_1_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat_1, 3, ROwf, 1, 1.d0, RO_1_sf, 1 )
      call rottrn(rotmat_1, ROwf, RO_1_sf, com_1)
c
      do i=1,3
         RM_1_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat_1, 3, ROwf, 1, 1.d0, RO_1_sf, 1 )
      call rottrn(rotmat_1, RMwf, RM_1_sf, com_1)

      do i=1,3
         RH1_1_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat_1, 3, R1wf, 1, 1.d0, R1_1_sf, 1 )
      call rottrn(rotmat_1, RH1wf, RH1_1_sf, com_1)
c
      do i=1,3
         RH2_1_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, R2wf, 1, 1.d0, R21sf, 1 )
      call rottrn(rotmat_1, RH2wf, RH2_1_sf, com_1)
c
c     prepare rotational matrix for water 2
c     obtain the SFF coordinates for H1, H2, and O of water 2
      call matpre(Eulang_2, rotmat_2)
      do i=1,3
         RO_2_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat_2, 3, ROwf, 1, 1.d0, RO_2_sf, 1 )
      call rottrn(rotmat_2, ROwf, RO_2_sf, com_2)
c
      do i=1,3
         RM_2_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat_2, 3, ROwf, 1, 1.d0, RO_2_sf, 1 )
      call rottrn(rotmat_2, RMwf, RM_2_sf, com_2)

c     call rottrn(rotmat2,R1wf,R12sf,com2)
      do i=1,3
         RH1_2_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, R1wf, 1, 1.d0, R12sf, 1 )
      call rottrn(rotmat_2, RH1wf, RH1_2_sf, com_2)
c
c     call rottrn(rotmat2,R2wf,R22sf,com2)
      do i=1,3
         RH2_2_sf(i)=0.d0
      enddo
c     call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, R2wf, 1, 1.d0, R22sf, 1 )
      call rottrn(rotmat_2, RH2wf, RH2_2_sf, com_2)
c
c
c ... calculate water dimer energies through q-TIP4P/F formula
      E_2H2O=0.d0
c ... O-O interaction
      roo=0.0d0
      rmm=0.0d0
      do i=1,3
        roo=roo+(RO_1_sf(i)-RO_2_sf(i))*(RO_1_sf(i)-RO_2_sf(i))
        rmm=rmm+(RM_1_sf(i)-RM_2_sf(i))*(RM_1_sf(i)-RM_2_sf(i))
      enddo
      rmm=dsqrt(rmm)
      roo4=roo*roo
      roo6=roo4*roo
      roo12=roo6*roo6
c ... LJ interaction between O-O
      sig2=sigoo*sigoo
      sig6=sig2*sig2*sig2
      sig12=sig6*sig6
      v_o2lj=4.0d0*epsoo*(sig12/roo12-sig6/roo6) !in Kelvin
c ... H-O, H-H and O-O Columbic interaction
      rhm1=zero
      rhm2=zero
      rhm3=zero
      rhm4=zero
      rhh1=zero
      rhh2=zero
      rhh3=zero
      rhh4=zero
      do i=1,3
        rhm1=rhm1+(RM_1_sf(i)-RH1_2_sf(i))*(RM_1_sf(i)-RH1_2_sf(i))
        rhm2=rhm2+(RM_1_sf(i)-RH2_2_sf(i))*(RM_1_sf(i)-RH2_2_sf(i))
        rhm3=rhm3+(RM_2_sf(i)-RH1_1_sf(i))*(RM_2_sf(i)-RH1_1_sf(i))
        rhm4=rhm4+(RM_2_sf(i)-RH2_1_sf(i))*(RM_2_sf(i)-RH2_1_sf(i))
c
        rhh1=rhh1+(RH1_1_sf(i)-RH1_2_sf(i))*(RH1_1_sf(i)-RH1_2_sf(i))
        rhh2=rhh2+(RH1_1_sf(i)-RH2_2_sf(i))*(RH1_1_sf(i)-RH2_2_sf(i))
        rhh3=rhh3+(RH2_1_sf(i)-RH1_2_sf(i))*(RH2_1_sf(i)-RH1_2_sf(i))
        rhh4=rhh4+(RH2_1_sf(i)-RH2_2_sf(i))*(RH2_1_sf(i)-RH2_2_sf(i))
      enddo
      rhm1=dsqrt(rhm1)
      rhm2=dsqrt(rhm2)
      rhm3=dsqrt(rhm3)
      rhm4=dsqrt(rhm4)
      rhh1=dsqrt(rhh1)
      rhh2=dsqrt(rhh2)
      rhh3=dsqrt(rhh3)
      rhh4=dsqrt(rhh4)
c
c     v_mhcolm: coulumbic term between M and H from different H2O in Hartree
      v_mhcolm=qm*qh*(1.d0/rhm1+1.d0/rhm2+1.d0/rhm3+1.d0/rhm4)
c     v_hhcolm:  ... between H and H ...
      v_hhcolm=qh*qh*(1.d0/rhh1+1.d0/rhh2+1.d0/rhh3+1.d0/rhh4)
c     v_mmcolm:  ... between M and M ...
      v_mmcolm=qm*qm*(1.d0/rmm)
c
      E_2H2O=v_o2lj+(v_mhcolm+v_mmcolm+v_hhcolm)*hr2k*br2ang
      return
      end
