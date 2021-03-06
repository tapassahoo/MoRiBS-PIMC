c 20130227 TZ add H^2 estimator to calculate heat capacity of the rotor
      program asymrho_driver
      implicit double precision(a-h,o-z)

      character argum*30

      call getarg(1,argum)
      read(argum,*)temprt
      call getarg(2,argum)
      read(argum,*)nslice
      call getarg(3,argum)
      read(argum,*)kmod
      call getarg(4,argum)
      read(argum,*)ithini
      call getarg(5,argum)
      read(argum,*)ithfnl
      call getarg(6,argum)
      read(argum,*)Bz
      call getarg(7,argum)
      read(argum,*)Bxy
      call getarg(8,argum)
      read(argum,*)maxj

      if(maxj.gt.876) stop 'maxj is larger than the limit of 876'

      call symrho(temprt,nslice,kmod,ithini,ithfnl,Bz,
     +                   Bxy,maxj)

      end
c----------------------------------------------------------------------
      subroutine symrho(temprt,nslice,kmod,ithini,ithfnl,Bz,
     +                   Bxy,maxj)
      implicit double precision(a-h,o-z)

      parameter(pi=3.14159265358979323846d+00,eps=1.d-16)
      parameter(zero=0.d0)
      parameter(boltz=0.6950356d0)
      parameter(maxfac=1754)
      dimension erotpr(0:181*361*361-1),erotsq(0:181*361*361-1),
     +          dlist(0:maxj,-maxj:maxj,0:180),rhopro(0:181*361*361-1),
     +          rhounq(0:180,0:360),eunq(0:180,0:360),
     +          esqunq(0:180,0:360),
     +          coslst(0:maxj,0:360)
      real*16 fact(0:maxfac)
      parameter(increm=10)
      character argum*30

c ... calculate tau
      beta=1.d0/(boltz*temprt)
      tau=beta/dfloat(nslice)
      write(6,'(''tau='',f10.5)')tau


c ... calculate factorial
      call calfac(fact,maxfac)

c ... prepare the list of wigner d function
      do j=0,maxj
        do k=-j,j
          do ith=ithini,ithfnl
            th=dfloat(ith)*pi/180.d0
            dlist(j,k,ith)=wigd(j,k,k,th,maxfac,fact)
          enddo
        enddo
      enddo

c ... prepare the cos list
      do k=0,maxj
        do icp=0,360
          cph=dfloat(icp)*pi/180.d0
          coslst(k,icp)=cos(k*cph)
        enddo
      enddo

c ... calculate partition function at tau and beta and energy and heat capacity at beta
      ztau=0.d0
      zbeta=0.d0
      Ebeta=0.d0
      Esqrt=0.d0
      do j=0,maxj
        do k=0,j
          if(mod(k,kmod).eq.0) then
            kgen=2-kdel(k,0)
            ejk=Bxy*j*(j+1)+(Bz-Bxy)*k*k
            ztau=ztau+kgen*exp(-tau*ejk)*(2*j+1)
            zbeta=zbeta+kgen*exp(-beta*ejk)*(2*j+1)
            Ebeta=Ebeta+kgen*exp(-beta*ejk)*(2*j+1)*ejk
            Esqrt=Esqrt+kgen*exp(-beta*ejk)*(2*j+1)*ejk*ejk
c           write(6,*)j,k,kgen
          endif
        enddo
      enddo
      Ebeta=Ebeta/zbeta
      Esqrt=Esqrt/zbeta
c ... convert to K
      Ebeta=Ebeta/boltz
c ... convert to K^2
      Esqrt=Esqrt/(boltz*boltz)
      Cv=(Esqrt-Ebeta*Ebeta)/(temprt*temprt)
      write(6,*)'ztau=',ztau
      write(6,*)'zbeta=',zbeta
      write(6,*)'Ebeta=',Ebeta,'K'
      write(6,*)'Esqrt=',Esqrt,'K^2'
      write(6,*)'Cv=',Cv,'Kb'

c ... calculate emax
      if(Bz.gt.Bxy) then
        emax=bxy*maxj*(maxj+1)+(Bz-Bxy)*maxj*maxj
      else
        emax=bxy*maxj*(maxj+1)
      endif

c ... percentage of the highest level in ztau
      pmax=(2*maxj+1)*exp(-tau*emax)/ztau
      if(pmax.gt.eps) then
        write(6,*)'pmax too large',pmax,'increase maxj'
        stop
      endif


      do 20 ith=ithini,ithfnl

      do 30 icp=0,360
      cph=dfloat(icp)*pi/180.d0
      rho=0.d0
      erot=0.d0
      esq=0.d0
      do j=0,maxj
        pre=dfloat(2*j+1)/(8.d0*pi*pi)
        do k=0,j
          if(mod(k,kmod).eq.0) then
            kgen=2-kdel(k,0)
            ejk=Bxy*j*(j+1)+(Bz-Bxy)*k*k
            rho=rho+pre*kgen*dlist(j,k,ith)*dcos(k*cph)*exp(-tau*ejk)
            erot=erot+
     +           pre*kgen*dlist(j,k,ith)*dcos(k*cph)*exp(-tau*ejk)*ejk
            esq=esq+
     +         pre*kgen*dlist(j,k,ith)*dcos(k*cph)*exp(-tau*ejk)*ejk*ejk
          endif
        enddo

      enddo
      erot=erot/rho
      esq=esq/rho

      rhounq(ith,icp)=rho
      eunq(ith,icp)=erot
      esqunq(ith,icp)=esq

c ... record
      do iph=0,360
        do ich=0,iph
          if(mod(iph+ich,360).eq.icp) then
            iphp=ich
            ichp=iph
            ind=ith*361*361+iph*361+ich
            ind2=ith*361*361+iphp*361+ichp
            rhopro(ind)=rhounq(ith,icp)
            erotpr(ind)=eunq(ith,icp)
            erotsq(ind)=esqunq(ith,icp)
            rhopro(ind2)=rhounq(ith,icp)
            erotpr(ind2)=eunq(ith,icp)
            erotsq(ind2)=esqunq(ith,icp)
          endif
        enddo
      enddo

   30 continue

      if(ith.lt.10) then
         write(argum,'(a,i1)')'rho.den00',ith
      elseif(ith.ge.10.and.ith.le.99) then
         write(argum,'(a,i2)')'rho.den0',ith
      elseif(ith.gt.99.and.ith.le.180) then
         write(argum,'(a,i3)')'rho.den',ith
      else
         write(6,*)'weird ith',ith
         stop
      endif

      lenarg=lastch(argum,30)

      write(6,'(a)')argum(1:lenarg)

c ... file 2 stores the regular output table
      open(2,file=argum(1:lenarg),status='unknown')

c ... file 3 stores the rotational propagator in the data block format
      open(3,file=argum(1:lenarg)//'_rho',status='unknown')
c ... file 4 stores the rotational energy estimator in the data block format
      open(4,file=argum(1:lenarg)//'_eng',status='unknown')
c ... file 7 stores the heat capacity estimator in the data block format
      open(7,file=argum(1:lenarg)//'_esq',status='unknown')

      if(ith.eq.0) then
        write(2,'(''# T='',f10.5,'' NSLICE='',I5,'' KMOD='',I5)')
     +      temprt,nslice,kmod
        write(2,'(a)')'# the  phi  chi       rho            engrot'
      endif

      do iph=0,360
        do ich=0,360
          ind=ith*361*361+iph*361+ich
          write(2,'(3(I5),3(1x,E15.8))')ith,iph,ich,rhopro(ind),
     +                                erotpr(ind),erotsq(ind)
          write(3,'(E15.8)')rhopro(ind)
          write(4,'(E15.8)')erotpr(ind)
          write(7,'(E15.8)')erotsq(ind)
        enddo
      enddo

      close(2,status='keep')
      close(3,status='keep')
      close(4,status='keep')
      close(7,status='keep')

   20 continue


      end
c------------------------------------------------------------------
      integer function kdel(i,j)
      implicit double precision(a-h,o-z)

c ... shift-up for the situation of i=j=0
c ... because of the shifting, the delta function
c ... is ILL-DEFINED for the case of i=j=-1000
      ii=i+1000
      jj=j+1000
      kdel=((ii+jj)-iabs(ii-jj))/((ii+jj)+iabs(ii-jj))

      return
      end
c------------------------------------------------------------------

      double precision function cplus(j,k)
      implicit double precision(a-h,o-z)

      parameter(one=1.0d+00,zero=0.0d+00)

c ... in case k runs out of the range
      if(k.ge.j.or.k.lt.(-j)) then
        cplus=zero
        return
      endif

      dj=dfloat(j)
      dk=dfloat(k)

      cplus=sqrt(dj*(dj+one)-dk*(dk+one))

      return
      end
c------------------------------------------------------------------

      double precision function cminus(j,k)
      implicit double precision(a-h,o-z)

      parameter(one=1.0d+00,zero=0.0d+00)

c ... in case k runs out of the range
      if(k.le.(-j).or.k.gt.j)then
        cminus=zero
        return
      endif

      dj=dfloat(j)
      dk=dfloat(k)

      cminus=sqrt(dj*(dj+one)-dk*(dk-one))

      return
      end

c------------------------------------------------------------------
      double precision function rotmat(j,k,kp,A,B,C)
      implicit double precision(a-h,o-z)

      rotmat=0.d0
      if(iabs(k-kp).ne.0.and.(iabs(k-kp).ne.2)) then
        return
      elseif(k.eq.kp) then
        rotmat=0.5d0*(A+C)*dfloat(j*(j+1))+(B-0.5d0*(A+C))*dfloat(k*k)
        return
      elseif(k.eq.(kp+2))then
        rotmat=0.25d0*(A-C)*cplus(j,kp)*cplus(j,kp+1)
        return
      elseif(k.eq.(kp-2))then
        rotmat=0.25d0*(A-C)*cminus(j,kp)*cminus(j,kp-1)
        return
      else
        stop 'wrong with rotmat'
      endif


      return
      end
c------------------------------------------------------------------
      double precision function rotma2(j,k,kp,A,B,C)
      implicit double precision(a-h,o-z)

      rotma2=0.d0
      if(iabs(k-kp).ne.0.and.(iabs(k-kp).ne.2)) then
        return
      elseif(k.eq.kp) then
        rotma2=0.5d0*(B+C)*dfloat(j*(j+1)-k*k)+A*dfloat(k*k)
        return
      elseif(k.eq.(kp+2))then
        rotma2=0.25d0*(B-C)*cplus(j,kp)*cplus(j,kp+1)
        return
      elseif(k.eq.(kp-2))then
        rotma2=0.25d0*(B-C)*cminus(j,kp)*cminus(j,kp-1)
        return
      else
        stop 'wrong with rotma2'
      endif


      return
      end

c------------------------------------------------------------------
      SUBROUTINE TQL(MD,N,Z,D,E)
C       Z(-1) A  Z  =D                                    
C       A = Z
C       EIGENVALUE         D(I)
C       EIGENFUNCTION      Z(J,I),J=1,N
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION   D(MD),E(MD),Z(MD,MD)
      EPS=1D-12
      NITER=50
      CALL TRED2(MD,N,Z,D,E)
      DO 10 I=2,N
  10  E(I-1)=E(I)
      F=0.0D0
      B=0.0D0
      E(N)=0.0D0
      DO 20 L=1,N
      J=0
      H=EPS*(DABS(D(L))+DABS(E(L)))
      LP1=L+1
      IF (B-H) 30,40,40
  30  B=H
  40  DO 50 M=L,N
      IF (DABS(E(M))-B) 60,60,50
  50  CONTINUE
  60  IF (M-L) 70,80,70
  70  IF (J-NITER) 90,100,90
  90  J=J+1
      P=(D(LP1)-D(L))/(2*E(L))
      R=DSQRT(P*P+1)
      IF (P) 110,111,111
  110 H=D(L)-E(L)/(P-R)
      GOTO 130
  111 H=D(L)-E(L)/(P+R)
  130 DO 140 I=L,N
  140 D(I)=D(I)-H
      F=F+H
      P=D(M)
      C=1.0D0
      S=0.0D0
      MM1=M-1
      IF (MM1-L) 270,280,280
  280 DO 120 LMIP=L,MM1
      I=L+MM1-LMIP
      IP1=I+1
      G=C*E(I)
      H=C*P
      IF (DABS(P)-DABS(E(I))) 160,170,170
  170 C=E(I)/P
      R=DSQRT(C*C+1.0D0)
      E(IP1)=S*P*R
      S=C/R
      C=1.0D0/R
      GOTO 180
  160 C=P/E(I)
      R=DSQRT(C*C+1)
      E(IP1)=S*E(I)*R
      S=1/R
      C=C/R
  180 P=C*D(I)-S*G
      D(IP1)=H+S*(C*G+S*D(I))
      DO 190 K=1,N
      H=Z(K,IP1)
      Z(K,IP1)=S*Z(K,I)+C*H
  190 Z(K,I)=C*Z(K,I)-S*H
  120 CONTINUE
  270 E(L)=S*P
      D(L)=C*P
      IF (DABS(E(L))-B) 80,80,70
  80  D(L)=D(L)+F
  20  CONTINUE
      DO 112 I=1,N
      IP1=I+1
      K=I
      P=D(I)
      IF (N-I) 230,230,300
  300 DO 210 J=IP1,N
      IF (D(J)-P) 220,210,210
  220 K=J
      P=D(J)
  210 CONTINUE
  230 IF (K-I) 240,112,240
  240 D(K)=D(I)
      D(I)=P
      DO 260 J=1,N
      P=Z(J,I)
      Z(J,I)=Z(J,K)
  260 Z(J,K)=P
  112 CONTINUE
      RETURN
  100 STOP '  FAIL'
      END
c------------------------------------------------------------------
      SUBROUTINE TRED2(MD,N,Z,D,E)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  D(MD),E(MD),Z(MD,MD)
      BETA=1D-20
      DO 20 NMIP2=2,N
      I=N+2-NMIP2
      IM1=I-1
      IM2=I-2
      L=IM2
      F=Z(I,IM1)
      G=0.0D0
      IF (L) 30,30,40
  40  DO 50 K=1,L
  50  G=G+Z(I,K)*Z(I,K)
  30  H=G+F*F
      IF (G-BETA) 60,60,70
  60  E(I)=F
      H=0.0D0
      GOTO 180
  70  L=L+1
      IF (F) 80,90,90
  90  E(I)=-DSQRT(H)
      G=E(I)
      GOTO 100
  80  E(I)=DSQRT(H)
      G=E(I)
 100  H=H-F*G
      Z(I,IM1)=F-G
      F=0.0D0
      DO 110 J=1,L
      Z(J,I)=Z(I,J)/H
      G=0.0D0
      DO 201 K=1,J
  201 G=G+Z(J,K)*Z(I,K)
      JP1=J+1
      IF (JP1-L) 130,130,140
  130 DO 120 K=JP1,L
  120 G=G+Z(K,J)*Z(I,K)
  140 E(J)=G/H
      F=F+G*Z(J,I)
  110 CONTINUE
      HH=F/(H+H)
      DO 160    J=1,L
      F=Z(I,J)
      E(J)=E(J)-HH*F
      G=E(J)
      DO 170 K=1,J
  170 Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  160 CONTINUE
  180 D(I)=H
  20  CONTINUE
      D(1)=0.0D0
      E(1)=0.0D0
      DO 190 I=1,N
      L=I-1
      IF (D(I)) 202,210,202
  202 IF (L) 210,210,220
  220 DO 230 J=1,L
      G=0.0D0
      DO 240 K=1,L
  240 G=G+Z(I,K)*Z(K,J)
      DO 250 K=1,L
  250 Z(K,J)=Z(K,J)-G*Z(K,I)
  230 CONTINUE
  210 D(I)=Z(I,I)
      Z(I,I)=1.0D0
      IF (L) 260,260,270
  270 DO 280 J=1,L
      Z(I,J)=0.0D0
  280 Z(J,I)=0.0D0
  260 CONTINUE
  190 CONTINUE
      return
      END
c------------------------------------------------------------------
      double precision function wigd(j,m,k,theta,maxfac,fact)
      implicit real*16(a-h,o-z)
      real*16 fact(0:maxfac)
      double precision theta
      parameter(eps=1.q-16)
c ... this function calculates the wigner d-matrix element.
c ... It takes Eq. 3.57 of Zare, 1988.

      pre1=sqrt(fact(j+k))*sqrt(fact(j-k))*sqrt(fact(j+m))*
     +     sqrt(fact(j-m))

c ... judge the upper bound of nu, the summing index
      nulow=max(0,k-m)
      nuup=min(j+k,j-m)
      thehlf=0.5q+00*theta
      wigd_temp=0.0q+00
c     wigd=0.0d+00

c ... summation over nu
      do nu=nulow,nuup
        denorm=fact(j-m-nu)*fact(j+k-nu)*fact(nu+m-k)*fact(nu)
     +          *(-1)**(nu)
        pre2=pre1/denorm
        cosfac=cos(thehlf)
        sinfac=-sin(thehlf)
        cosfac=cosfac**(2*j+k-m-2*nu)
        sinfac=sinfac**(m-k+2*nu)
        wigd_temp=wigd_temp+pre2*cosfac*sinfac
      enddo

      wigd=dble(wigd_temp)

      return
      end
c-----------------------------------------------------------------
      subroutine calfac(fact,maxfac)
      implicit double precision(a-h,o-z)
      real*16 fact(0:maxfac)

      fact(0)=1.0q+00

      do i=1,maxfac
        fact(i)=fact(i-1)*dfloat(i)
      enddo
      return
      end
c     Subroutine bubble_sort
c       this routine sorts the given data
c
      subroutine bubble_sort(data,count)
c     
c     argument:  count is a positive integer
      integer count
c     argument:  data is an array of size count
      double precision data(count)
c
c     local variables:
      integer i
c       how many times we have passed through the array
      integer pass
c       flag variable: 1 if sorted; 0 if not sorted  
      integer sorted
c       temporary variable used for swapping       
      double precision temp

      pass = 1
 1    continue
      sorted = 1
      do 2 i = 1,count-pass
        if(data(i) .gt. data(i+1)) then
          temp = data(i)
          data(i) = data(i+1)
          data(i+1) = temp
          sorted = 0
        endif
 2    continue
      pass = pass +1
      if(sorted .eq. 0) goto 1
      return
      end
c---------------------------------------------------------------------
      integer function lastch(line,len)
      implicit double precision(a-h,o-z)
      character*(*) line

      do i=len,1,-1
         if(line(i:i).ne.' ') then
           lastch=i
           return
         endif
      enddo
      return
      end
