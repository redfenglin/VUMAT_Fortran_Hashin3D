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
c
      parameter (zero=0.d0,one=1.d0,two=2.d0)
C Temporary arrays
      real stiffmat(9),EnggConstant(6)
      real nu12,nu13,nu23,nu21,nu31,nu32,Lamina_thickness,density
      Lamina_thickness=props(1) ! lamina thickness
      Integral_points=props(2)! either 3,5 or 7
      ele_thickness=Lamina_thickness/Integral_points
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
      G_II=props(22)! mode 2 energy release rate
    
      alpha=props(23)! for model given in NASA paper = 0.5
      beta=props(24)! for model given in NASA paper = 0.5
      Sbeta=props(25)! defined by us for matrix compression 

c
        delta=one - (nu32*nu23) - (nu31*nu13) - 
     1               (nu12*nu21)-(two*nu13*nu21*nu32)
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
        print * ,'t=0, Stiffmat',stiffmat
        print * ,'t=0, all properties', E110,E220,E330,G120,
     1   G130,G230, nu12, nu21,nu13,nu31,nu23,nu32, G_ft, G_fc, G_mt,
     2   G_mc,G_II, alpha, beta, Sbeta
c
c-----------------------------------------------------------------
        do k=1,nblock
          if ( stepTime .eq. zero ) then    
C             print *,'t=0,strain update,before strainOld', stateOld(:,1:6)
C             print *,'t=0,strain update,before strainNew', stateNew(:,1:6)
            call updatestrain(nblock,strainInc,stateOld(:,1:6),
     1                    stateNew(:,1:6))
C             print *,'t=0,strain update,after strainOld', stateOld(:,1:6)
C             print *,'t=0,strain update,after strainNew', stateNew(:,1:6)
c
c-----------------------------------------------------------------
c
C             print *,'t=0,stress update before Stiffmat',stiffmat
C             print *,'t=0,stress update before stressOld', stressOld
C             print *,'t=0,stress update before damage(0/1)=',stateOld(:,13)
C             print *,'t=0,stress update before strain ',stateNew(:,1:6)
C             print *,'t=0,stress update before stressOld', stressOld
C             print *,'t=0,stress update before cri strian ',stateNew(:,14:19)
C             print *,'t=0,stress update before cri stress ',stateNew(:,20:25)
            call stress_update(stiffmat,stateNew(:,1:6),
     1         stressNew,stateNew(:,14:19),stateNew(:,20:25),
     2         stateOld(:,13))
C             print *,'t=0,stress update after Stiffmat',stiffmat
C             print *,'t=0,stress update after stressOld', stressOld
C             print *,'t=0,stress update after damage(0/1)=',stateOld(:,13)
C             print *,'t=0,stress update after strain ',stateNew(:,1:6)
C             print *,'t=0,stress update after stressNew', stressNew
C             print *,'t=0,stress update after cri strian ',stateNew(:,14:19)
C             print *,'t=0,stress update after cri stress ',stateNew(:,20:25)
C             print *, 't=0 complete'
C         print *, 'updated stress: ', stressNew
          else
c
c-----------------------------------------------------------------
c

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
c
c-----------------------------------------------------------------
C             print *,'strain update,before strainOld', stateOld(:,1:6)
C             print *,'strain update,before strainNew', stateNew(:,1:6)
            call updatestrain(nblock,strainInc,stateOld(:,1:6),
     1                      stateNew(:,1:6))
C             print *,'strain update,after strainOld', stateOld(:,1:6)
C             print *,'strain update,after strainNew', stateNew(:,1:6)
c      ! since strainOld is not given, we store strains in SDV
c
            if (stateOld(k,13).eq.zero) then
c
              print *, 'damage=0'
              delta=one - (nu32*nu23) - (nu31*nu13) - 
     1             (nu12*nu21) - 2*nu13*nu21*nu32
              stiffmat(1)=E110*(one - nu32*nu23)/delta
              stiffmat(2)=E220*(one - nu13*nu31)/delta
              stiffmat(3)=E330*(one - nu12*nu21)/delta
c      
              stiffmat(4)=(E110*((nu31*nu23)+nu21)) / delta
              stiffmat(5)=(E220*((nu12*nu31)+nu32)) / delta
              stiffmat(6)=(E110*((nu21*nu32)+nu31)) / delta
c      
              stiffmat(7)=G120
              stiffmat(8)=G230
              stiffmat(9)=G130
c-----------------------------------------------------------------
              print *,'stress update before Stiffmat',stiffmat
              print *,'stress update before stressOld', stressOld
              print *,'stress update before damage(0/1)=',stateOld(k,13)
              print *,'stress update before strain ',stateNew(k,1:6)
              print *,'stress update before stressOld', stressOld
              print *,'stress update before cristrian',stateNew(k,14:19)
              print *,'stress update before cristress',stateNew(k,20:25)
              call stress_update(stiffmat,stateNew(k,1:6),
     1             stressNew,stateNew(k,14:19),stateNew(k,20:25),
     2             stateNew(k,26:28))
              print *,'stress update after Stiffmat',stiffmat
              print *,'stress update after stressOld', stressOld
              print *,'stress update after damage(0/1)=',stateOld(k,13)
              print *,'stress update after strain ',stateNew(k,1:6)
              print *,'stress update after stressNew', stressNew
              print *,'stress update after cri strian',stateNew(k,14:19)
              print *,'stress update after cri stress',stateNew(k,20:25)
            ! stress again if damage happened
c
c-----------------------------------------------------------------
              print *,'damage check before d=',stateNew(k,13)
              print *,'damage check before d1,d2,d3', stateNew(k,26:28)
              print *,'damage check before stress=',stressNew
              print *,'damage check before strainNew',stateNew(k,1:6)
              print *,'damage check before prop', E110,E220,E330,G120, 
     5               G230, G130,X_c, X_t, Y_c, Y_t, S_l, S_t
              call damage_check(stateNew(k,13),stateNew(k,26:28),
     3                         stateNew(k,1:6), stressNew,
     4                         E110, E220, E330, G120, G230, G130,
     5                        X_c, X_t, Y_c, Y_t, S_l, S_t)
              print *,'damage check after damage',stateNew(k,13)
              print *,'damage check after d1,d2,d3', stateNew(k,26:28)
              print *,'damage check after stress=',stressNew
              print *,'damage check after strain ',stateNew(k,1:6)
              print *,'damage check after prop', E110,E220,E330,G120, 
     5               G230, G130,X_c, X_t, Y_c, Y_t, S_l, S_t

            end if
c-----------------------------------------------------------------
            if (stateNew(k,13).ne.zero) then
              print *,'damage occured'
c-----------------------------------------------------------------
c
              print *,'comp engg const before strainNew',stateNew(k,1:6)
              print *,'comp engg const before d1-d3',stateNew(k,26:28)
              print *,'comp engg const before charLength',charLength(k)
              print *,'comp engg const before el_thick',ele_thickness
              print *,'comp engg const before c_strn',stateNew(k,14:19)
              print *,'comp engg const before c_strss',stateNew(k,20:25)
              print *,'comp engg const before prop',E110,E220,E330,G120,
     5               G230, G130,G_ft,G_fc,G_mt,G_mc,G_II
              print *,'comp engg const before EnggConstant',EnggConstant
              call compute_Engg_constant(stateNew(k,1:6),
     1             stateNew(k,26:28),charLength(k),ele_thickness,
     2             stateNew(k,14:19),stateNew(k,20:25),
     2             E110,E220,E330,G120,G230,G130,
     3             G_ft,G_fc,G_mt,G_mc,G_II,
     4             EnggConstant)
              print *,'comp engg const after strainNew',stateNew(k,1:6)
              print *,'comp engg const after d1-d3',stateNew(k,26:28)
              print *,'comp engg const after charLength',charLength(k)
              print *,'comp engg const after el_thick',ele_thickness
              print *,'comp engg const after c_strn',stateNew(k,14:19)
              print *,'comp engg const after c_strss',stateNew(k,20:25)
              print *,'comp engg const after prop',E110,E220,E330,G120,
     5               G230, G130,G_ft,G_fc,G_mt,G_mc,G_II
              print *,'comp engg const after EnggConstant',EnggConstant
c
c-----------------------------------------------------------------
c
              print *,'comp dmg var before damageOld',stateOld(k,7:13)
              print *,'comp dmg var before damageNew',stateNew(k,7:13)
              print *,'comp dmg var before prop',E110,E220,E330,G120,
     5               G230, G130,alpha,beta
              print *,'comp dmg var before EnggConstant',EnggConstant
              call compute_damage_variable(stateOld(k,7:13),
     1             stateNew(k,7:13),EnggConstant,E110,E220,E330
     2             G120,G130,G230,alpha,beta)
              print *,'comp dmg var after damageOld',stateOld(k,7:13)
              print *,'comp dmg var after damageNew',stateNew(k,7:13)
              print *,'comp dmg var after prop',E110,E220,E330,G120,
     5               G230, G130,alpha,beta
              print *,'comp dmg var after EnggConstant',EnggConstant
c
c-----------------------------------------------------------------
c
              print *,'comp Stif bfore damageNew',stateNew(k,7:13)
              print *,'comp Stif bfore EnggConstant',EnggConstant
              print *,'comp Stif bfore nu',nu12,nu13,nu23,nu21,nu31,nu32
              print *,'comp Stif bfore stiffmat',stiffmat
              call compute_stiffmat(stateNew(k,7:13),EnggConstant,
     1             nu12,nu13,nu23,nu21,nu31,nu32,stiffmat)
              print *,'comp Stif after damageNew',stateNew(k,7:13)
              print *,'comp Stif after EnggConstant',EnggConstant
              print *,'comp Stif after nu',nu12,nu13,nu23,nu21,nu31,nu32
              print *,'comp Stif after stiffmat',stiffmat
c
c-----------------------------------------------------------------
c
              print *,'stress update before Stiffmat',stiffmat
              print *,'stress update before stressOld', stressOld
              print *,'stress update before damage(0/1)=',stateOld(k,13)
              print *,'stress update before strain ',stateNew(k,1:6)
              print *,'stress update before stressOld', stressOld
              print *,'stress update before cristrian',stateNew(k,14:19)
              print *,'stress update before cristress',stateNew(k,20:25)
              call stress_update(stiffmat,stateNew(k,1:6),
     1             stressNew,stateNew(k,14:19),stateNew(k,20:25),
     2             stateOld(k,13))
              print *,'stress update after Stiffmat',stiffmat
              print *,'stress update after stressOld', stressOld
              print *,'stress update after damage(0/1)=',stateOld(k,13)
              print *,'stress update after strain ',stateNew(k,1:6)
              print *,'stress update after stressOld', stressOld
              print *,'stress update after cri strian',stateNew(k,14:19)
              print *,'stress update after cri stress',stateNew(k,20:25)
c
            end if
          end if 
        end do
c          
        call energy_update(nblock,stressOld,stressNew,
     1       strainInc, enerInternNew, enerInternOld,density)
      return
      end 
c
C--------------------------------------------------------
C update strain values, stored in SDV.
C--------------------------------------------------------
      subroutine updatestrain(nblock,strainInc,
     1 strainOld,strainNew)
      include 'vaba_param.inc'
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
C--------------------------------------------------------
C Damage initiation check and critical stress and strain
C--------------------------------------------------------
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
C--------------------------------------------------------
C compute Engg elastic constants 
C--------------------------------------------------------
      subroutine compute_Engg_constant(strainNew,damage
     1 L_c,ele_thickness,
     2 cri_strainNew,cri_stressNew,
     2 E110,E220,E330,G120,G230,G130,
     3 G_ft,G_fc,G_mt,G_mc,G_II,
     4 EnggConstant)
c
      include 'vaba_param.inc'
      real temp,G,EnggConstant(6)
      dimension stiffmat(9),
     2  strainNew(6),
     3  cri_strainNew(6),
     4  cri_stressNew(6),
     5  damage(3)
      parameter (zero=0.d0, half=0.5d0)
c
      e11=strainNew(1)
      e22=strainNew(2)
      e33=strainNew(3)
      e12=strainNew(4)
      e23=strainNew(5)
      e13=strainNew(6)
c
      e11c=cri_strainNew(1)
      e22c=cri_strainNew(2)
      e33c=cri_strainNew(3)
      e12c=cri_strainNew(4)
      e23c=cri_strainNew(5)
      e13c=cri_strainNew(6)
c
      t1=cri_stressNew(1)
      t2=cri_stressNew(2)
      t3=cri_stressNew(3)
      t12=cri_stressNew(4)
      t23=cri_stressNew(5)
      t13=cri_stressNew(6)
c
      if (damage(1) .ne. zero) then
        if (e11.gt.zero) then
          G=G_ft
        else
          G=G_fc
        end if
        de11=abs(e11) - e11c
        temp=(L_c*t1*t1*half)/G
        print *,'for d1, c_strss,G,E110,temp,de11',t1,G,E110,temp,de11
        EnggConstant(1)=(E110*(t1-(temp*de11)))/
     *    (t1-((temp-E110)*de11))
      else 
        EnggConstant(1)=E110
      end if
c
      if (damage(2) .ne. zero) then
        if (e22.gt.zero) then
          G=G_mt
        else
          G=G_mc
        end if
        de22=abs(e22) - e22c
        temp=(L_c*t2*t2*half)/G
        print *,'for d2, c_strss,G,E220,temp,de22',t2,G,E220,temp,de22
        EnggConstant(2)=(E220*(t2-(temp*de22)))/
     *    (t2-((temp-E220)*de22))
c
        de12=abs(e12) - e12c
        temp=(L_c*t12*t12*half)/G_II
        EnggConstant(4)=(G120*(t12-(temp*de12)))/
     *    (t12-((temp-G120)*de12))
c      
        de23=abs(e23) - e23c
        temp=(L_c*t23*t23*half)/G_II
        EnggConstant(5)=(G230*(t23-(temp*de23)))/
     *  (t23-((temp-G230)*de23))
      else
        EnggConstant(2)=E220
        EnggConstant(4)=G120
        EnggConstant(5)=G230
      end if 
c
      print *, 'd3',damage(3)
      if (damage(3) .ne. zero) then
        if (e33.gt.zero) then
          G=G_mt
        else
          G=G_mc
        end if
        de33=abs(e33) - e33c
        temp=(L_c*t3*t3*half)/G
        print *,'for d2, c_strss,G,E220,temp,de22',t3,G,E330,temp,de33
        EnggConstant(3)=(E330*(t3-(temp*de33)))/
     *    (t3-((temp-E330)*de33))
c     
        de23=abs(e23) - e23c
        temp=(L_c*t23*t23*half)/G_II
        EnggConstant(5)=(G230*(t23-(temp*de23)))/
     *    (t23-((temp-G230)*de23))
c
        de13=abs(e13) - e13c
        temp=(L_c*t13*t13*half)/G_II
        EnggConstant(6)=(G130*(t13-(temp*de13)))/
     *    (t13-((temp-G130)*de13))
c
      else
        print *,'no d3'
        EnggConstant(3)=E330
        EnggConstant(5)=G230
        EnggConstant(6)=G130
      end if 
      return
      end
c
C--------------------------------------------------------
C compute damage variable 
C--------------------------------------------------------
      subroutine compute_damage_variable(damageOld,
     1    damageNew,EnggConstant,E110,E220,E330,
     2    G120,G130,G230,alpha,beta)
c      
      include 'vaba_param.inc'
      real temp
      dimension EnggConstant(6),
     2  damageOld(7),damageNew(7)

      parameter (one=1.d0)
c
        E11=EnggConstant(1)
        E22=EnggConstant(2)
        E33=EnggConstant(3)
        G12=EnggConstant(4)
        G23=EnggConstant(5)
        G13=EnggConstant(6)
c
        dn1=one-(E11/E110)
        dn2=one-(E22/E220)
        dn3=one-(E33/E330) 
        ds12=one-(G12/G120)
        ds13=one-(G13/G130)
        ds23=one-(G23/G230)
c
        damageNew(1)=max(dn1,damageOld(1))
        damageNew(2)=max(dn2,damageOld(2))
        damageNew(3)=max(dn3,damageOld(3))
        temp=one-damageNew(1)
c
        damageNew(4)=max(damageOld(4),(one-temp*
     1    (one-alpha*damageNew(2))*(one-beta*ds12)))
c      
        damageNew(5)=max(damageOld(5),
     1                 one-temp*(one-alpha*damageNew(2))*
     2                 (one-alpha*damageNew(3))*(one-beta*ds23))
c     
        damageNew(6)=max(damageOld(6),
     1                 one-temp*
     2                 (one-alpha*damageNew(3))*
     3                 (one-beta*ds13))
        damageNew(7)=one
      return
      end
c      
C--------------------------------------------------------
C Compute Stiffness matrix
C--------------------------------------------------------
      subroutine compute_stiffmat(damageNew,
     1 EnggConstant,nu12,nu13,nu23,nu21,nu31,nu32,
     2 stiffmat)
c      
      include 'vaba_param.inc'
      dimension stiffmat(9),
     2  damageNew(7),
     5  EnggConstant(6)
        E11=EnggConstant(1)
        E22=EnggConstant(2)
        E33=EnggConstant(3)
        G12=EnggConstant(4)
        G23=EnggConstant(5)
        G13=EnggConstant(6)
c      
        dd1=one-damageNew(1)
        dd2=one-damageNew(2)
        dd3=one-damageNew(3)
        dd4=one-damageNew(4)
        dd5=one-damageNew(5)
        dd6=one-damageNew(6)
c
        delta=one - (dd2*dd3*nu32*nu23) - 
     1               (dd1*dd3*nu31*nu13) - 
     2               (dd1*dd2*nu12*nu21) -
     3               (two*dd1*dd2*dd3*nu13*nu21*nu32) 
c
        stiffmat(1)=(dd1*E11*(one-(dd2*dd3*nu32*nu23)))/delta
        stiffmat(2)=(dd2*E22*(one-(dd1*dd3*nu13*nu31)))/delta
        stiffmat(3)=(dd3*E33*(one-(dd1*dd2*nu12*nu21)))/delta
c      
        stiffmat(4)=(dd1*dd2*E11*((dd3*nu31*nu23)+nu21))/delta
        stiffmat(5)=(dd2*dd3*E22*((dd1*nu12*nu31)+nu32))/delta
        stiffmat(6)=(dd1*dd3*E11*((dd2*nu21*nu32)+nu31))/delta
c      
        stiffmat(7)=dd4*G12
        stiffmat(8)=dd5*G23
        stiffmat(9)=dd6*G13
c
      return 
      end
c      
C----------------------------------------------------------
C update stresses
C----------------------------------------------------------
      subroutine stress_update(stiffmat,strainNew,stressNew,
     1  cri_strainNew,cri_stressNew,damage)
c
      include 'vaba_param.inc'
      parameter(zero=0.d0)
      dimension stiffmat(9),damage(3),
     1  strainNew(6),
     2  stressNew(6),
     3  cri_strainNew(6),
     4  cri_stressNew(6)
c
        d1=damage(1)
        d2=damage(2)
        d3=damage(3)
c
        c11=stiffmat(1)
        c22=stiffmat(2)
        c33=stiffmat(3)
        c12=stiffmat(4)
        c23=stiffmat(5)
        c13=stiffmat(6)
        c44=stiffmat(7)
        c55=stiffmat(8)
        c66=stiffmat(9)
c
        e11=strainNew(1)
        e22=strainNew(2)
        e33=strainNew(3)
        e12=strainNew(4)
        e23=strainNew(5)
        e13=strainNew(6)

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
        cri_strainNew(1)=abs(e11)
c
        cri_stressNew(1)=abs(stressNew(1))
      end if

      if (d2 .eq. zero) then ! if second criteria fails
        cri_strainNew(2)=abs(e22)
        cri_strainNew(4)=abs(e12)
        cri_strainNew(5)=abs(e23)
c
        cri_stressNew(2)=abs(stressNew(2))
        cri_stressNew(4)=abs(stressNew(4))
        cri_stressNew(5)=abs(stressNew(5))
      end if

      if (d3 .eq. zero) then ! if thired criteria fails
        cri_strainNew(3)=abs(e33)
        cri_strainNew(6)=abs(e13)
        cri_strainNew(5)=abs(e23)
c
        cri_stressNew(3)=abs(stressNew(3))
        cri_stressNew(6)=abs(stressNew(6))
        cri_stressNew(5)=abs(stressNew(5))
      end if
      return
      end
c
C----------------------------------------------------------
C update internal energy and dissipation energy
C----------------------------------------------------------
      subroutine energy_update(nblock,stressOld,stressNew,
     1 strainInc, enerInternNew, enerInternOld,density)
c      
      include 'vaba_param.inc'
      parameter (half=0.5d0,two=2.d0)
      dimension stressOld(nblock,6),
     2  stressNew(nblock,6),
     3  strainInc(nblock,6),
     4  enerInternNew(nblock),enerInternOld(nblock)
      real new_energy
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
