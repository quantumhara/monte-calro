      function ran4(idum)
      integer idum
      real ran4
cu    uses psdes
      integer idums,irword,itemp,jflmsk,jflone,lword
      real ftemp
      equivalence (itemp,ftemp)
      save idums,jflone,jflmsk
      data idums /0/, jflone /16#3f800000/, jflmsk /16#007fffff/
      if(idum.lt.0)then
        idums=-idum
        idum=1
      endif
      irword=idum
      lword=idums
      call psdes(lword,irword)
      itemp=ior(jflone,iand(jflmsk,irword))
      ran4=ftemp-1.0
      idum=idum+1
      return
      end
