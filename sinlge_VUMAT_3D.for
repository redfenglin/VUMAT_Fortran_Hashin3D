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
c-----------------------------------------------------------------
c Extract material properties from props = Engineering Constant
c-----------------------------------------------------------------
      parameter (zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,small=1.d-12)
C Temporary arrays
      real stiff(9)
      real E110,E220,E330,G120,G130,G230,X_t,X_c,Y_t,Y_
      real nu12,nu13,nu23,nu21,nu31,nu32,Lamina_thickness,density
      real e11,e22,e33,e12,e23,e13
      real delta,dd1,dd2,dd3,dd4,dd5,dd6, new_energy
      real dn1,dn2,dn3,ds12,ds13,ds23_22,ds23_33,ds23,temp,
     1    de11,de22,de33,de12,de23,de13,t11,t22,t33,t12,t23,t13,G
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
      Y_c=props(15)
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
        delta=one - (nu32*nu23) - (nu31*nu13) - 
     1              (nu12*nu21)-(two*nu13*nu21*nu32)
        stiff(1)=(E110*(one - nu32*nu23)) / delta
        stiff(2)=(E220*(one - nu13*nu31)) / delta
        stiff(3)=(E330*(one - nu12*nu21)) / delta
c      
        stiff(4)=(E110*((nu31*nu23)+nu21)) / delta
        stiff(5)=(E220*((nu12*nu31)+nu32)) / delta
        stiff(6)=(E110*((nu21*nu32)+nu31)) / delta
c      
        stiff(7)=G120
        stiff(8)=G230
        stiff(9)=G130
c-----------------------------------------------------------------
c-----------------------------------------------------------------
        do k=1,nblock
          if ( stepTime .eq. zero ) then  ! step time =0
c-----------------------------------------------------------------
c t=0, update strain 
c
           stateNew(k,1)=strainInc(k,1)+stateOld(k,1)
           stateNew(k,2)=strainInc(k,2)+stateOld(k,2)
           stateNew(k,3)=strainInc(k,3)+stateOld(k,3)
           stateNew(k,4)=strainInc(k,4)+stateOld(k,4)
           stateNew(k,5)=strainInc(k,5)+stateOld(k,5)
           stateNew(k,6)=strainInc(k,6)+stateOld(k,6)
c
c-----------------------------------------------------------------
c t=0  stress update
c
            e11=stateNew(k,1)
            e22=stateNew(k,2)
            e33=stateNew(k,3)
            e12=stateNew(k,4)
            e23=stateNew(k,5)
            e13=stateNew(k,6)
c
         stressNew(k,1)= (stiff(1)*e11) + (stiff(4)*e22) +(stiff(6)*e33) 
         stressNew(k,2)= (stiff(4)*e11) + (stiff(2)*e22) +(stiff(5)*e33) 
         stressNew(k,3)= (stiff(6)*e11) + (stiff(5)*e22) +(stiff(3)*e33)
c
           stressNew(k,4)= stiff(7)*e12
           stressNew(k,5)= stiff(8)*e23 
           stressNew(k,6)= stiff(9)*e13
c t=0 end of if
c-----------------------------------------------------------------
          else ! step time .ne. zero
c
c  t != 0 all SDV's update
            stateNew(k,1)=stateold(k,1)! strain
            stateNew(k,2)=stateold(k,2)
            stateNew(k,3)=stateold(k,3)
            stateNew(k,4)=stateold(k,4)
            stateNew(k,5)=stateold(k,5)
            stateNew(k,6)=stateold(k,6)
c
            stateNew(k,7)=stateold(k,7)!dn1 damage variables
            stateNew(k,8)=stateold(k,8)!dn2
            stateNew(k,9)=stateold(k,9)!dn3
            stateNew(k,10)=stateold(k,10)!ds12
            stateNew(k,11)=stateold(k,11)!ds22
            stateNew(k,12)=stateold(k,12)!ds13
            stateNew(k,13)=stateold(k,13)! d <- overall damage
c
            stateNew(k,14)=stateold(k,14)! cri_strain
            stateNew(k,15)=stateold(k,15)
            stateNew(k,16)=stateold(k,16)
            stateNew(k,17)=stateold(k,17)
            stateNew(k,18)=stateold(k,18)
            stateNew(k,19)=stateold(k,19)
c
            stateNew(k,20)=stateold(k,20)! cri_stress
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
c  t!= 0   update strain 
c
           stateNew(k,1)=strainInc(k,1)+stateOld(k,1)
           stateNew(k,2)=strainInc(k,2)+stateOld(k,2)
           stateNew(k,3)=strainInc(k,3)+stateOld(k,3)
           stateNew(k,4)=strainInc(k,4)+stateOld(k,4)
           stateNew(k,5)=strainInc(k,5)+stateOld(k,5)
           stateNew(k,6)=strainInc(k,6)+stateOld(k,6)
c
c-----------------------------------------------------------------
           if (stateOld(k,13).eq.zero) then ! if no damage
c
              print *, 'damage=0'
c-----------------------------------------------------------------
c  t!=0, d=0   stress update
c
            e11=stateNew(k,1)
            e22=stateNew(k,2)
            e33=stateNew(k,3)
            e12=stateNew(k,4)
            e23=stateNew(k,5)
            e13=stateNew(k,6)
c
         stressNew(k,1)= (stiff(1)*e11) + (stiff(4)*e22) +(stiff(6)*e33) 
         stressNew(k,2)= (stiff(4)*e11) + (stiff(2)*e22) +(stiff(5)*e33) 
         stressNew(k,3)= (stiff(6)*e11) + (stiff(5)*e22) +(stiff(3)*e33)
c
           stressNew(k,4)= stiff(7)*e12
           stressNew(k,5)= stiff(8)*e23 
           stressNew(k,6)= stiff(9)*e13
c
c-----------------------------------------------------------------
c t != 0 d=0    damage_check
       if (e11 .gt. zero) then
          if (((e11/(X_t/E110))**two) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,26)=one
          end if
        else
          if (((e11/(X_c/E110))**two) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,26)=one
          end if
        end if
c
        if (e22 .gt. zero) then
          if ((((e22/(Y_t/E220))**two) + ((two*e12/(S_l/G120))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,27)=one
          end if 
        else
          if ((((e22/(Y_c/E220))**two) + ((two*e12/(S_l/G120))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,27)=one
          end if 
        end if
c
        if (e33 .gt. zero) then
          if ((((e33/(Y_t/E330))**two) + ((two*e13/(S_t/G130))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,28)=one
          end if
        else
          if ((((e33/(Y_c/E330))**two) + ((two*e13/(S_t/G130))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,28)=one
          end if
        end if
c
c-----------------------------------------------------------------
c t!=0 d=0    update cri_strain and cri_stress
        if (stateNew(k,26) .eq. zero) then ! if first criteria fails
          stateNew(k,14)=e11
          stateNew(k,20)=stressNew(k,1)
        end if
c
        if (stateNew(k,27) .eq. zero) then ! if second criteria fails
          stateNew(k,15)=e22
          stateNew(k,17)=e12
          stateNew(k,18)=e23
c
          stateNew(k,21)=stressNew(k,2)
          stateNew(k,23)=stressNew(k,4)
          stateNew(k,24)=stressNew(k,5)
        end if
c
        if (stateNew(k,28) .eq. zero) then ! if thired criteria fails
          stateNew(k,16)=e33
          stateNew(k,19)=e13
          stateNew(k,18)=e23
c
          stateNew(k,22)=stressNew(k,3)
          stateNew(k,25)=stressNew(k,6)
          stateNew(k,24)=stressNew(k,5)
        end if
      end if ! till here damage .eq. zero
c t!=0 and (d=0 end of if loop) 
c-----------------------------------------------------------------
        if (stateNew(k,13).ne.zero) then ! damage occures
          print *,''
          print *,'***************damage occured************'
c d!=0  damage variables computation 
c      real dn1,dn2,dn3,ds12,ds13,ds23_22,ds23_33,ds23,temp,
c     3    de11,de22,de33,de12,de23,de13,t11,t22,t33,t12,t23,t13,G
        dn1=zero
        dn2=zero
        dn3=zero
        ds12=zero
        ds23_22=zero
        ds23_33=zero
        ds23=zero
        ds13=zero
        if (stateNew(k,26) .eq. one )then ! first failure criteria
          if (stateNew(k,1).gt.zero) then
            G=G_ft
          else
            G=G_fc
          end if
          t11=stateNew(k,20)
          de11=stateNew(k,1) - stateNew(k,14)
          if (de11 .gt. small) then
            if (abs(de11) .lt. (two*G/(charLength(k)*abs(t11)))) then
            temp=(t11/(de11*E110))
     1        -((charLength(k)*t11*t11*half)/(G*E110))
            dn1=one/(one+temp)
            else
            dn1=one
            end if
          end if
        end if
c
        if (stateNew(k,27) .eq. one )then 
          if (stateNew(k,2).gt.zero) then
            G=G_mt
          else
            G=G_mc
          end if
          t22=stateNew(k,21)
          t12=stateNew(k,23)
          t23=stateNew(k,24)
          de22=stateNew(k,2) - stateNew(k,15)
          de12=stateNew(k,4) - stateNew(k,17)
          de23=stateNew(k,5) - stateNew(k,18)
          if (de22 .gt. small) then
            if (abs(de22) .lt. (two*G/(charLength(k)*abs(t22)))) then
            temp=(t22/(de22*E220))
     1       -((charLength(k)*t22*t22*half)/(G*E220))
            dn2=one/(one+temp)
            else
            dn2=one
            end if
          end if
          if (de12 .gt. small) then
           if (abs(de12) .lt. (two*G_II/(charLength(k)*abs(t12)))) then
             temp=(t12/(de12*G120))
     1       -((charLength(k)*t12*t12*half)/(G_II*G120))
            ds12=one/(one+temp)
            else
            ds12=one
            end if
          end if
          if (de23 .gt. small) then
           if (abs(de23) .lt. (two*G_II/(charLength(k)*abs(t23)))) then
            temp=(t23/(de23*G230))
     1      - ((charLength(k)*t23*t23*half)/(G_II*G230))
            ds23_22=one/(one+temp)
            else
            ds23_22=one
            end if
          end if
        end if
c
        if (stateNew(k,28) .eq. one )then
          if (stateNew(k,3).gt.zero) then
            G=G_mt
          else
            G=G_mc
          end if
          t33=stateNew(k,22)
          t23=stateNew(k,24)
          t13=stateNew(k,25)
          de33=stateNew(k,3) - stateNew(k,16)
          de23=stateNew(k,5) - stateNew(k,18)
          de13=stateNew(k,6) - stateNew(k,19)
          if (abs(de33) .gt. small) then 
            if (abs(de33) .lt. (two*G/(charLength(k)*abs(t33)))) then
            temp=(t33/(de33*E330))
     1      -((charLength(k)*t33*t33*half)/(G*E330))
            dn3=one/(one+temp)
            else
            dn3=one
            end if
          end if
          if (abs(de23).gt.small) then
           if (abs(de23) .lt. (two*G_II/(charLength(k)*abs(t23)))) then
            temp=(t23/(de23*G230))
     1       -((charLength(k)*t23*t23*half)/(G_II*G230))
            ds23_33=one/(one+temp)
           else
            ds23_33=one
           end if
          end if
          if ((abs(de13) .gt. small)) then
           if (abs(de13) .lt. (two*G_II/(charLength(k)*abs(t13)))) then
            temp=(t13/(de13*G130))
     1        -((charLength(k)*t13*t13*half)/(G_II*G130))
            ds13=one/(one+temp)
           else
            ds13=one
           end if
          end if
        end if
        ds23=max(ds23_22,ds23_33)
c
        stateNew(k,7)=max(dn1,stateOld(k,7))
        stateNew(k,8)=max(dn2,stateOld(k,8))
        stateNew(k,9)=max(dn3,stateOld(k,9))
        temp=one-stateNew(k,7)
c
        stateNew(k,10)=max(stateOld(k,10),(one-(temp*
     1    (one-alpha*stateNew(k,8))*(one - beta*ds12))))
c
        stateNew(k,11)=max(stateOld(k,11),
     1                 (one-(temp*(one-alpha*stateNew(k,8))*
     2                 (one-alpha*stateNew(k,9))*(one-beta*ds23))))
c
        stateNew(k,12)=max(stateOld(k,12),
     1                 (one-(temp*
     2                 (one-alpha*stateNew(k,9))*
     3                 (one-beta*ds13))))
        stateNew(k,13)=one
c-----------------------------------------------------------------
c d=0    update stiffness matrix
c      real delta,dd1,dd2,dd3,dd4,dd5,dd6
c      
        dd1=one-stateNew(k,7)
        dd2=one-stateNew(k,8)
        dd3=one-stateNew(k,9)
        dd4=one-stateNew(k,10)
        dd5=one-stateNew(k,11)
        dd6=one-stateNew(k,12)
c
        delta=one - (dd2*dd3*nu32*nu23) - (dd1*dd3*nu31*nu13) 
     2        - (dd1*dd2*nu12*nu21) - (two*dd1*dd2*dd3*nu13*nu21*nu32) 
c
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
        print *, 'stiff',stiff
        print *, 'damage', stateNew(k,7:12)
c
c-----------------------------------------------------------------
c    stress update
            e11=stateNew(k,1)
            e22=stateNew(k,2)
            e33=stateNew(k,3)
            e12=stateNew(k,4)
            e23=stateNew(k,5)
            e13=stateNew(k,6)
c
           stressNew(k,1)= (stiff(1)*e11) + (stiff(4)*e22) +(stiff(6)*e33) 
           stressNew(k,2)= (stiff(4)*e11) + (stiff(2)*e22) +(stiff(5)*e33) 
           stressNew(k,3)= (stiff(6)*e11) + (stiff(5)*e22) +(stiff(3)*e33)
c
           stressNew(k,4)= stiff(7)*e12
           stressNew(k,5)= stiff(8)*e23 
           stressNew(k,6)= stiff(9)*e13
c
c-----------------------------------------------------------------
c again check the damage criterial failure 
c     damage_check  
        if (stateNew(k,26) .eq. zero) then ! if first criteria fails
            if (e11 .gt. zero) then
               if (((e11/(X_t/E110))**two) .gt. one) then
                 stateNew(k,13)=one
                 stateNew(k,26)=one
                end if
            else
               if (((e11/(X_c/E110))**two) .gt. one) then
                stateNew(k,13)=one
                stateNew(k,26)=one
               end if
            end if
        end if
c
      if (stateNew(k,27) .eq. zero) then ! if first criteria fails
        if (e22 .gt. zero) then
          if ((((e22/(Y_t/E220))**two) + ((two*e12/(S_l/G120))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,27)=one
          end if 
        else
          if ((((e22/(Y_c/E220))**two) + ((two*e12/(S_l/G120))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,27)=one
          end if 
        end if
      end if
c
      if (stateNew(k,28) .eq. zero) then ! if first criteria fails
        if (e33 .gt. zero) then
          if ((((e33/(Y_t/E330))**two) + ((two*e13/(S_t/G130))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,28)=one
          end if
        else
          if ((((e33/(Y_c/E330))**two) + ((two*e13/(S_t/G130))**two) 
     1       + ((two*e23/(S_t/G230))**two)) .gt. one) then
            stateNew(k,13)=one
            stateNew(k,28)=one
          end if
        end if
      end if
c-----------------------------------------------------------------
c    update cri_strain and cri_stress
        if (stateNew(k,26) .eq. zero) then ! if first criteria fails
          stateNew(k,14)=e11
          stateNew(k,20)=stressNew(k,1)
        end if
c
        if (stateNew(k,27) .eq. zero) then ! if second criteria fails
          stateNew(k,15)=e22
          stateNew(k,17)=e12
          stateNew(k,18)=e23
c
          stateNew(k,21)=stressNew(k,2)
          stateNew(k,23)=stressNew(k,4)
          stateNew(k,24)=stressNew(k,5)
        end if
c
        if (stateNew(k,28) .eq. zero) then ! if thired criteria fails
          stateNew(k,16)=e33
          stateNew(k,19)=e13
          stateNew(k,18)=e23
c
          stateNew(k,22)=stressNew(k,3)
          stateNew(k,25)=stressNew(k,6)
          stateNew(k,24)=stressNew(k,5)
        end if

      end if ! damage occured
c
c    Eneergy update
c         real new_energy
       new_energy=half*((stressOld(k,1)+stressNew(k,1))*strainInc(k,1)+
     3   (stressOld(k,2) + stressNew(k,2))*strainInc(k,2) +
     4   (stressOld(k,3) + stressNew(k,3))*strainInc(k,3) +
     5   (stressOld(k,4) + stressNew(k,4))*strainInc(k,4)*two +
     6   (stressOld(k,5) + stressNew(k,5))*strainInc(k,5)*two +
     7   (stressOld(k,6) + stressNew(k,6))*strainInc(k,6)*two )
c
      enerInternNew(k)=enerInternOld(k) + (new_energy/density(k))
c
c-----------------------------------------------------------------
          end if ! time .ne. zero
        end do ! end do loop
      return
      end 