      function ran3(idum)
      integer idum
      integer mbig,mseed,mz
c     real mbig,mseed,mz
      real ran3,fac
      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
c     parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig)
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
c     real mj,mk,ma(55)
      save iff,inext,inextp,ma
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=abs(mseed-abs(idum))
        mj=mod(mj,mbig)
        ma(55)=mj

        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz)mk=mk+mbig
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)

      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end
