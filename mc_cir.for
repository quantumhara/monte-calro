c  *********************************************************************
c  for monte calro method
c  calculating the pi with fractuation: f*r^2=pi*r^2
c  counting up the points in the area of circle
c  using the only one seed : dependent running with previous seed
c  *********************************************************************
      implicit real*8 (a-h,o-z)
      real*8 ran0,am
c  dimension  **********************************************************
      dimension farea(100000),avefarea(100000),variance(100000),
     1     standev(100000),difarea(100000)
c  open unit  **********************************************************
      open(unit=10,file='10.dat')
      open(unit=100,file='dfor100.dat')
c  reading starting values  ********************************************
      rewind 10
      read(10,*)
      read(10,*)seed
      read(10,*)n,m,l
c  set starting values  ************************************************
      pi=3.14159265358979323d0
      eps=10.d0**(-8.d0)
      zero=0.d0
c  main  ***************************************************************
      idum=seed
      totm=m
      totl=l

      do i3=1,l,1
         do i1=1,m,1
            incir=zero
            do i2=1,n,1
               coords1=ran0(idum)
               coords2=ran0(idum)
               slength=sqrt(coords1*coords1+coords2*coords2)
               if (slength<=1) then
                  incir=incir+1
               end if
            end do
            farea(i1)=4.d0*real(incir)/real(n)
         end do

         avefarea(i3)=zero
         do i1=1,m,1
            avefarea(i3)=avefarea(i3)+farea(i1)/totm
         end do

         variance(i3)=zero
         do i1=1,m,1
            variance(i3)=variance(i3)
     1           +(avefarea(i3)-farea(i1))*(avefarea(i3)-farea(i1))
     1           /(totm-1.d0)
         end do

         standev(i3)=sqrt(variance(i3))
         difarea(i3)=abs(avefarea(i3)-pi)

c  write data  *********************************************************
      write(6,*) avefarea(i3),difarea(i3),standev(i3),variance(i3),i3
      write(100,*) avefarea(i3),difarea(i3),standev(i3),variance(i3),i3
c      pause'**************************'
      end do

      estave=zero
      do i3=1,l,1
         estave=estave+avefarea(i3)/totl
      end do

      estvariance=zero
      do i3=1,l,1
         estvariance=estvariance
     1        +(estave-avefarea(i3))*(estave-avefarea(i3))/(totl-1.d0)
      end do

      eststandev=sqrt(estvariance)
      estdif=abs(estave-pi)

c  write data  *********************************************************
      write(6,700) n,m,l,seed,estave,estdif,eststandev,estvariance
      write(100,700) n,m,l,seed,estave,estdif,eststandev,estvariance
c      pause'**************************'

601   format(1(1pe15.7))
602   format(2(1pe15.7))
604   format(4(1pe15.7))
700   format('#',1x,'monte calro calculation',/,
     1 '#',1x,'random points=',i10,3x,'ind trial=',i10,3x,
     1 'total traial=',i10,3x,'seed=',1pe15.7,/,
     1 '#',1x,'ave_area ',8x,'ab_df_area',11x,'stan_dev ',10x,
     1 'variance',/,
     1  1pe15.7,2x,1pe15.7,6x,1pe15.7,4x,1pe15.7)

      stop
      end



c  function ************************************************************
      function ran0(idum)
      integer idum,ia,im,iq,ir,mask
      real*8 ran0,am
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     *mask=123459876)
      integer k
      idum=ieor(idum,mask)
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      ran0=am*idum
      idum=ieor(idum,mask)
      return
      end
