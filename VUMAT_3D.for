c-----------------------------------------------------------------
c TO USE THIS SUBROUTINE, USE SOLID ELEMENT C3D8R
c USE 24 USER DEFINE MATERIAL PROPERTIES 
c-----------------------------------------------------------------
      subroutine vumat (
C Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension coordMp(nblock,*), charLength(nblock), props(nprops),
     1     density(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr), 
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),
     1     fieldNew(nblock,nfieldv),
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
C     
      character*80 cmname

c
C        e11=stateNew(k,1)=strainNew(k,1)
C        e22=stateNew(k,2)=strainNew(k,2)
C        e33=stateNew(k,3)=strainNew(k,3)
C        e12=stateNew(k,4)=strainNew(k,4)
C        e23=stateNew(k,5)=strainNew(k,5)
C        e13=stateNew(k,6)=strainNew(k,6)
c
C        dn1=stateNew(k,7)=damageNew(k,1)
C        dn2=stateNew(k,8)=damageNew(k,2)
C        dn3=stateNew(k,9)=damageNew(k,3)
C        ds12=stateNew(k,10)=damageNew(k,4)
C        ds23=stateNew(k,11)=damageNew(k,5)
C        ds13=stateNew(k,12)=damageNew(k,6)
C        d=stateNew(k,13)=damageNew(k,7)
C        d=0 or 1
c
C        e11c=stateNew(k,14)=cri_strainNew(k,1)
C        e22c=stateNew(k,15)=cri_strainNew(k,2)
C        e33c=stateNew(k,16)=cri_strainNew(k,3)
C        e12c=stateNew(k,17)=cri_strainNew(k,4)
C        e23c=stateNew(k,18)=cri_strainNew(k,5)
C        e13c=stateNew(k,19)=cri_strainNew(k,6)
c
C        t1=stateNew(k,20)=cri_stressNew(k,1)
C        t2=stateNew(k,21)=cri_stressNew(k,2)
C        t3=stateNew(k,22)=cri_stressNew(k,3)
C        t12=stateNew(k,23)=cri_stressNew(k,4)
C        t23=stateNew(k,24)=cri_stressNew(k,5)
C        t13=stateNew(k,25)=cri_stressNew(k,6)
c
c-----------------------------------------------------------------
c Extract material properties from props = Engineering Constant
c-----------------------------------------------------------------
      parameter (zero=0.d0,one=1.d0,two=2.d0)
C Temporary arrays
      real stiffmat(9)
      common stiffmat
      real nu12,nu13,nu23,nu21,nu31,nu32,Lamina_thickness,density
      Lamina_thickness=props(1)
      Integral_points=props(2)
      E110=props(3)
      E220=props(4)
      E330=props(5)
      G120=props(6)
      G130=props(7)
      G230=props(8)
      nu12=props(9)
      nu13=props(10)
      nu23=props(11)
      nu21=nu12*(E220/E110)
      nu31=nu13*(E330/E110)
      nu32=nu23*(E330/E220)
      X_t=props(12)
      X_c=props(13)
      Y_t=props(14)
      y_c=props(15)
      S_l=props(16)
      S_t=props(17)
      G_ft=props(18)
      G_fc=props(19)
      G_mt=props(20)
      G_mc=props(21)
      G_II=props(22)
      alpha=props(23)
      beta=props(24)
      Sbeta=props(25)
c-----------------------------------------------------------------
c-----------------------------------------------------------------
        delta=one - (nu32*nu23) - (nu31*nu13) - 
     1              (nu12*nu21)-(two*nu13*nu21*nu32)
        stiffmat(1)=(E110*(one - nu32*nu23)) / delta
        stiffmat(2)=(E220*(one - nu13*nu31)) / delta
        stiffmat(3)=(E330*(one - nu12*nu21)) / delta
c      
        stiffmat(4)=(E110*((nu31*nu23)+nu21)) / delta
        stiffmat(5)=(E220*((nu12*nu31)+nu32)) / delta
        stiffmat(6)=(E110*((nu21*nu32)+nu31)) / delta
c      
        stiffmat(7)=G120
        stiffmat(8)=G230
        stiffmat(9)=G130
c-----------------------------------------------------------------
c-----------------------------------------------------------------
        do k=1,nblock
          if ( stepTime .eq. zero ) then    
            call updatestrain(nblock,strainInc,stateOld(:,1:6),
     1                    stateNew(:,1:6))
c-----------------------------------------------------------------
            call stress_update(stateNew(:,1:6),
     1         stressNew,stateNew(:,14:19),stateNew(:,20:25),
     2         stateOld(:,13))
          else ! step time .eq. zero
            stateNew(k,1)=stateold(k,1)
            stateNew(k,2)=stateold(k,2)
            stateNew(k,3)=stateold(k,3)
            stateNew(k,4)=stateold(k,4)
            stateNew(k,5)=stateold(k,5)
            stateNew(k,6)=stateold(k,6)
c
            stateNew(k,7)=stateold(k,7)
            stateNew(k,8)=stateold(k,8)
            stateNew(k,9)=stateold(k,9)
            stateNew(k,10)=stateold(k,10)
            stateNew(k,11)=stateold(k,11)
            stateNew(k,12)=stateold(k,12)
            stateNew(k,13)=stateold(k,13)
c
            stateNew(k,14)=stateold(k,14)
            stateNew(k,15)=stateold(k,15)
            stateNew(k,16)=stateold(k,16)
            stateNew(k,17)=stateold(k,17)
            stateNew(k,18)=stateold(k,18)
            stateNew(k,19)=stateold(k,19)
c
            stateNew(k,20)=stateold(k,20)
            stateNew(k,21)=stateold(k,21)
            stateNew(k,22)=stateold(k,22)
            stateNew(k,23)=stateold(k,23)
            stateNew(k,24)=stateold(k,24)
            stateNew(k,25)=stateold(k,25)
c
            stateNew(k,26)=stateold(k,26)! d1
            stateNew(k,27)=stateold(k,27)! d2
            stateNew(k,28)=stateold(k,28)! d3
c-----------------------------------------------------------------
            call updatestrain(nblock,strainInc,stateOld(:,1:6),
     1                      stateNew(:,1:6))
c-----------------------------------------------------------------
            if (stateOld(k,13).eq.zero) then ! if no damage
c
              print *, 'damage=0'
c-----------------------------------------------------------------
              call stress_update(stateNew(k,1:6),
     1             stressNew,stateNew(k,14:19),stateNew(k,20:25),
     2             stateNew(k,26:28))
c-----------------------------------------------------------------
              call damage_check(stateNew(k,13),stateNew(k,26:28),
     3                         stateNew(k,1:6), stressNew,
     4                         E110, E220, E330, G120, G230, G130,
     5                        X_c, X_t, Y_c, Y_t, S_l, S_t)
            end if
c-----------------------------------------------------------------
            if (stateNew(k,13).ne.zero) then
              print *,''
              print *,'***************damage occured************'
              print *,''
C-----------------------------------------------------------------
              print *,'comp dmg var before cri_stress',stateNew(k,20:25)
              call compute_damage_variable(stateNew(k,26:28),
     1             stateOld(k,7:13),stateNew(k,7:13),stateNew(k,14:19),
     2             stateNew(k,20:25),stateNew(k,1:6),E110,E220,E330,
     3             G120,G230,G130,alpha,beta,charLength,
     4             G_ft,G_fc,G_mt,G_mc,G_II)
              print *,'comp dmg var after damageNew',stateNew(k,7:13)
              print *,''
c-----------------------------------------------------------------
              print *,'comp Stif bfore damageNew',stateNew(k,7:13)
              print *,'comp Stif bfore nu',nu12,nu13,nu23,nu21,nu31,nu32
              print *,'comp Stif bfore stiffmat',stiffmat
              call compute_stiffmat(stateNew(k,7:13),E110,E220,E330,
     1             G120,G130,G230,nu12,nu13,nu23,nu21,nu31,nu32)
              print *,'comp Stif after stiffmat',stiffmat
c-----------------------------------------------------------------
              call stress_update(stateNew(k,1:6),
     1             stressNew,stateNew(k,14:19),stateNew(k,20:25),
     2             stateOld(k,26:28))
c-----------------------------------------------------------------
            end if
          end if 
        end do
c-----------------------------------------------------------------
        call energy_update(nblock,stressOld,stressNew,
     1       strainInc, enerInternNew, enerInternOld,density)
      return
      end 
c
c
      subroutine updatestrain(nblock,strainInc,
     1 strainOld,strainNew)
      include 'vaba_param.inc'
      integer nblock
      dimension strainOld(nblock,6),strainNew(nblock,6),
     *          strainInc(nblock,6)
      do k=1,nblock 
        strainNew(k,1)=strainInc(k,1)+strainOld(k,1)
        strainNew(k,2)=strainInc(k,2)+strainOld(k,2)
        strainNew(k,3)=strainInc(k,3)+strainOld(k,3)
        strainNew(k,4)=strainInc(k,4)+strainOld(k,4)
        strainNew(k,5)=strainInc(k,5)+strainOld(k,5)
        strainNew(k,6)=strainInc(k,6)+strainOld(k,6)
      end do
      return
      end
c
c
      subroutine damage_check(d,damage,
     1 strainNew,stressNew,
     2 E110,E220,E330,G120,G230,G130,
     3 X_c,X_t,Y_c,Y_t,S_l,S_t)
      include 'vaba_param.inc'
c
      dimension strainNew(6),stressNew(6),
     2  cri_strainNew(6),
     3  cri_stressNew(6), damage(3)
c
      parameter (zero=0.d0, one=1.d0,two=2.d0)
c
      real e11,e22,e33,e12,e23,e13,s11,s22,s33,s12,s23,s13,d,
     2 E110,E220,E330,G120,G230,G130,X_c,X_t,Y_c,Y_t,S_l,S_t
c
        e11=strainNew(1)
        e22=strainNew(2)
        e33=strainNew(3)
        e12=strainNew(4)
        e23=strainNew(5)
        e13=strainNew(6)
c
        s11=stressNew(1)
        s22=stressNew(2)
        s33=stressNew(3)
        s12=stressNew(4)
        s23=stressNew(5)
        s13=stressNew(6)
c
        if (e11 .gt. zero) then
          if (((e11/(X_t/E110))**two) .gt. one) then
            d=one
            damage(1)=one
          end if
        else
          if (((e11/(X_c/E110))**two) .gt. one) then
            d=one
            damage(1)=one
          end if
        end if
c
        if (e22 .gt. zero) then
          if ((((e22/(y_t/E220))**two) + ((two*e12/(S_l/G120))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            d=one
            damage(2)=one
          end if 
        else
          if ((((e22/(y_c/E220))**two) + ((two*e12/(S_l/G120))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            d=one
            damage(2)=one
          end if 
        end if
c
        if (e33 .gt. zero) then
          if ((((e33/(y_t/E330))**two) + ((two*e13/(S_t/G130))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            d=one
            damage(3)=one
          end if
        else
          if ((((e33/(y_c/E330))**two) + ((two*e13/(S_t/G130))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            d=one
            damage(3)=one
          end if
        end if
c
      return
      end
c
c
      subroutine compute_damage_variable(damage,damageOld,
     1    damageNew,cri_strainNew,cri_stressNew,strainNew,
     2    E110,E220,E330,G120,G230,G130,alpha,beta,lc,
     3    G_ft,G_fc,G_mt,G_mc,G_II)
c
      include 'vaba_param.inc'
      dimension damage(3),damageOld(7),damageNew(7),
     1    cri_strainNew(6),cri_stressNew(6),strainNew(6)
c
      real dn1,dn2,dn3,ds12,ds13,ds23_22,ds23_33,ds23,temp,E110,E220,
     2    E330,G120,G230,G130,alpha,beta,lc,G_ft,G_fc,G_mt,G_mc,G_II,
     3    de11,de22,de33,de12,de23,de13,t11,t22,t33,t12,t23,t13
c
      parameter (one=1.d0,zero=0.d0,half=0.5d0,small=1.d-8, two=2.d0)

        dn1=zero
        dn2=zero
        dn3=zero
        ds12=zero
        ds23_22=zero
        ds23_33=zero
        ds23=zero
        ds13=zero
        print *, 'all inside:'
        print *, 'damage:',damage
        print *, 'cri_strainNew:',cri_strainNew
        print *, 'cri_stressNew:',cri_stressNew
C         print *, 'prop:',E110,E220,E330,G120,G230,G130,alpha,
C      1       beta,lc,G_ft,G_fc,G_mt,G_mc,G_II
        if (damage(1) .eq. one )then ! first failure criteria
          if (strainNew(1).gt.zero) then
            G=G_ft
          else
            G=G_fc
          end if
          t11=cri_stressNew(1)
          de11=abs(strainNew(1) - cri_strainNew(1))
          print *, 'de11:',de11
          if (de11 .gt. small) then
            temp=one+(t11/(de11*E110))-((lc*t11*t11*half)/(G*E110))
            dn1=one/temp
            print *, 'temp,dn1:',temp,dn1
          end if
        end if
c
        if (damage(2) .eq. one )then
          if (strainNew(2).gt.zero) then
            G=G_mt
          else
            G=G_mc
          end if
          t22=cri_stressNew(2)
          t12=cri_stressNew(4)
          t23=cri_stressNew(5)
          de22=abs(strainNew(2) - cri_strainNew(2))
          de12=abs(strainNew(4) - cri_strainNew(4))
          de23=abs(strainNew(5) - cri_strainNew(5))
          print *, 'de22,de12,de23:',de22,de12,de23
          if (de22 .gt. small) then
            temp=one+(t22/(de22*E220))-((lc*t22*t22*half)/(G*E220))
            dn2=one/temp
            print *, 'temp,dn2:',temp,dn2
          end if
          if (de12 .gt. small) then
            temp=one+(t12/(de12*G120))-((lc*t12*t12*half)/(G_II*G120))
            ds12=one/temp
            print *, 'temp,ds12:',temp,ds12
          end if
          if (de23 .gt. small) then
            temp=one+(t23/(de23*G230))-((lc*t23*t23*half)/(G_II*G230))
            ds23_22=one/temp
            print *, 'temp,ds23_22:',temp,ds23_22
          end if
        end if
c
        if (damage(3) .eq. one )then
          if (strainNew(3).gt.zero) then
            G=G_mt
          else
            G=G_mc
          end if
          t33=cri_stressNew(3)
          t23=cri_stressNew(5)
          t13=cri_stressNew(6)
          de33=strainNew(3) - cri_strainNew(3)
          de23=strainNew(5) - cri_strainNew(5)
          de13=strainNew(6) - cri_strainNew(6)
          print *, 'de33,de23,de13:',de33,de23,de13
          if (abs(de33) .gt. small) then
            temp=one+(t33/(de33*E330))-((lc*t33*t33*half)/(G*E330))
            dn3=one/temp
            print *, 'temp,dn3:',temp,dn3
          end if
          if ((abs(de23) .gt. small) .and. (de23 .lt. (two*G/lc*t23)))then
            temp=one+(t23/(de23*G230))-((lc*t23*t23*half)/(G_II*G230))
            ds23_33=one/temp
            print *, 'temp,ds23_33:',temp,ds23_33
          end if
          if ((abs(de13) .gt. small) .and. (de13 .lt. (two*G/lc*t13)))then
            temp=one+(t13/(de13*G130))-((lc*t13*t13*half)/(G_II*G130))
            ds13=one/temp
            print *, 'temp,ds13:',temp,ds13
          end if
        end if
        ds23=max(ds23_22,ds23_33)
        print *, 'ds23:',ds23
c
        damageNew(1)=max(dn1,damageOld(1))
        damageNew(2)=max(dn2,damageOld(2))
        damageNew(3)=max(dn3,damageOld(3))
        temp=one-damageNew(1)
c
        damageNew(4)=max(damageOld(4),(one-(temp*
     1    (one-alpha*damageNew(2))*(one-beta*ds12))))
c
        damageNew(5)=max(damageOld(5),
     1                 (one-(temp*(one-alpha*damageNew(2))*
     2                 (one-alpha*damageNew(3))*(one-beta*ds23))))
c     
        damageNew(6)=max(damageOld(6),
     1                 (one-(temp*
     2                 (one-alpha*damageNew(3))*
     3                 (one-beta*ds13))))
        damageNew(7)=one
      return
      end
c
c
      subroutine compute_stiffmat(damageNew,
     1 E110,E220,E330,G120,G130,G230,nu12,nu13,nu23,nu21,nu31,nu32)
c      
      include 'vaba_param.inc'
      dimension damageNew(7)
      real stiff(9)
      real delta,dd1,dd2,dd3,dd4,dd5,dd6,E110,E220,E330,G120,G130,G230,
     1 nu12,nu13,nu23,nu21,nu31,nu32
      parameter (one=1.d0,two=2.d0)
      common stiff
c      
        dd1=one-damageNew(1)
        dd2=one-damageNew(2)
        dd3=one-damageNew(3)
        dd4=one-damageNew(4)
        dd5=one-damageNew(5)
        dd6=one-damageNew(6)
c
        print *, 'dd',dd1,dd2,dd3,dd4,dd5,dd6
        delta=one - (dd2*dd3*nu32*nu23) - (dd1*dd3*nu31*nu13) 
     2        - (dd1*dd2*nu12*nu21) - (two*dd1*dd2*dd3*nu13*nu21*nu32) 
c
        print *, 'inside delta=',delta
        print *, 'inside stiff=',stiff
        stiff(1)=(dd1*E110*(one-(dd2*dd3*nu32*nu23)))/delta
        stiff(2)=(dd2*E220*(one-(dd1*dd3*nu13*nu31)))/delta
        stiff(3)=(dd3*E330*(one-(dd1*dd2*nu12*nu21)))/delta
c      
        stiff(4)=(dd1*dd2*E110*((dd3*nu31*nu23)+nu21))/delta
        stiff(5)=(dd2*dd3*E220*((dd1*nu12*nu31)+nu32))/delta
        stiff(6)=(dd1*dd3*E110*((dd2*nu21*nu32)+nu31))/delta
c      
        stiff(7)=dd4*G120
        stiff(8)=dd5*G230
        stiff(9)=dd6*G130
        print *,'after stiff=',stiff
c
      return 
      end
c
c
      subroutine stress_update(strainNew,stressNew,
     1  cri_strainNew,cri_stressNew,damage)
c
      include 'vaba_param.inc'
      parameter(zero=0.d0)
      dimension damage(3),
     1  strainNew(6),
     2  stressNew(6),
     3  cri_strainNew(6),
     4  cri_stressNew(6)
c
      real stiff(9)
      real d1,d2,d3,c11,c22,c33,c12,c23,c13,c44,c55,c66
      real e11,e22,e33,e12,e23,e13
      common stiff
c
        d1=damage(1)
        d2=damage(2)
        d3=damage(3)
c
        c11=stiff(1)
        c22=stiff(2)
        c33=stiff(3)
        c12=stiff(4)
        c23=stiff(5)
        c13=stiff(6)
        c44=stiff(7)
        c55=stiff(8)
        c66=stiff(9)
c
        e11=strainNew(1)
        e22=strainNew(2)
        e33=strainNew(3)
        e12=strainNew(4)
        e23=strainNew(5)
        e13=strainNew(6)
c
        stressNew(1)= (c11*e11) + (c12*e22) +(c13*e33) 
        stressNew(2)= (c12*e11) + (c22*e22) +(c23*e33) 
        stressNew(3)= (c13*e11) + (c23*e22) +(c33*e33)
c      
        stressNew(4)= c44*e12
        stressNew(5)= c55*e23 
        stressNew(6)= c66*e13
c-----------------------------------------------------------------
        ! update critical strain and stress here
c-----------------------------------------------------------------
        if (d1 .eq. zero) then ! if first criteria fails
          cri_strainNew(1)=e11
          cri_stressNew(1)=stressNew(1)
        end if
c
        if (d2 .eq. zero) then ! if second criteria fails
          cri_strainNew(2)=e22
          cri_strainNew(4)=e12
          cri_strainNew(5)=e23
c
          cri_stressNew(2)=stressNew(2)
          cri_stressNew(4)=stressNew(4)
          cri_stressNew(5)=stressNew(5)
        end if
c
        if (d3 .eq. zero) then ! if thired criteria fails
          cri_strainNew(3)=e33
          cri_strainNew(6)=e13
          cri_strainNew(5)=e23
c
          cri_stressNew(3)=stressNew(3)
          cri_stressNew(6)=stressNew(6)
          cri_stressNew(5)=stressNew(5)
        end if
      return
      end
c
c
      subroutine energy_update(nblock,stressOld,stressNew,
     1 strainInc, enerInternNew, enerInternOld,density)
c      
      include 'vaba_param.inc'
      real new_energy,density
      parameter (half=0.5d0,two=2.d0)
      dimension stressOld(nblock,6),
     2  stressNew(nblock,6),
     3  strainInc(nblock,6),
     4  enerInternNew(nblock),enerInternOld(nblock)
c
      do k=1,nblock
       new_energy=half*((stressOld(k,1)+stressNew(k,1))*strainInc(k,1)+
     3   (stressOld(k,2) + stressNew(k,2))*strainInc(k,2) +
     4   (stressOld(k,3) + stressNew(k,3))*strainInc(k,3) +
     5   (stressOld(k,4) + stressNew(k,4))*strainInc(k,4)*two +
     6   (stressOld(k,5) + stressNew(k,5))*strainInc(k,5)*two +
     7   (stressOld(k,6) + stressNew(k,6))*strainInc(k,6)*two )
c
      enerInternNew(k)=enerInternOld(k) + (new_energy/density)
      end do
c
      return
      end
