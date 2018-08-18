!    Copyright (C) 2018 Quentin Vanhaelen

subroutine Gillespie_LDM()
use  parameters
use  variables
implicit none
real    ::  r1,r2, a_t, F, lambda1, lambda2,x, tau
real    ::  check_t
integer ::   j, n, k,i,l
integer ::  ongoing, index,x_min
 mol_added(1:1) = .TRUE. 
ongoing  = 0
lambda1 = 0.0
lambda2 = 0.0
i = 0
do n = 1, ntotal
if  (t >= 600) then
   Y(20) =60
       if (mol_added(1)) then
      a_nu(36) =  c_nu(36)*Y(1)*Y(20)
          do j=36,M_react_channel 
           p_sum(j) = a_nu(j)+p_sum(j-1)
          end do
          mol_added(1) = .FALSE.
        end if
end if
 CALL RANDOM_NUMBER(r1)
 CALL RANDOM_NUMBER(r2)
if (ongoing .eq. 0) then
tau = (1./p_sum(M_react_channel))*log(1./r1)
  t = t+  tau
else 
  lambda1 = t
  lambda2 =  react_t_delay(1) 
 a_t = p_sum(M_react_channel)*(lambda2-lambda1)
  F = 1 - EXP(-a_t)
 do while( F .lt. r1) 
     mu = mu_delayed(1)
    if ( consuming(mu))  then
       do l = 1, stochio_number_product(mu)
        number_species_product(l,mu)%species =   number_species_product(l,mu)%species+   update_product(l,mu) 
       end do
    else
      do l = 1, stochio_number(mu)
        number_species(l,mu)%species =   number_species(l,mu)%species+   update(l,mu) 
      end do
   end if
  mu_yes = vector_mu_tilde(mu)
if (mu_yes .le.3) then
if (mu_yes.eq.1) then
!Equation number=6
              a_nu(5)=c_nu(5)*Y(15)
              a_nu(6)=c_nu(6)*Y(15)
              a_nu(4)=c_nu(4)*Y(16)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
             x_min=4
else if (mu_yes.eq.2) then
!Equation number=7
              a_nu(7)=c_nu(7)*Y(14)
              a_nu(5)=c_nu(5)*Y(15)
              a_nu(6)=c_nu(6)*Y(15)
             x_min=5
else
!Equation number=9
              a_nu(9)=c_nu(9)*Y(12)
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
             x_min=9
end if
else
if (mu_yes.eq.4) then
!Equation number=10
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
              a_nu(7)=c_nu(7)*Y(14)
             x_min=7
else if (mu_yes.eq.5) then
!Equation number=16
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(9)=c_nu(9)*Y(12)
             x_min=9
else
end if
end if











































































do j=x_min,M_react_channel 
      p_sum(j) = a_nu(j)+p_sum(j-1)
 end do
 lambda1 = react_t_delay(1)
  if (ongoing.gt.1) then 
   do j = 2,ongoing
    react_t_delay(j-1)  =  react_t_delay(j) 
    mu_delayed(j-1)  =  mu_delayed(j) 
   end do
   react_t_delay(ongoing)  =  10000.0 
   mu_delayed(ongoing)  =  0 
  end if
  if (ongoing.eq.1) then
    react_t_delay(1)  =  10000.0
    mu_delayed(1)  =  0 
  end if
   ongoing  = ongoing -1
  if (ongoing == 0) exit
  if (ongoing > 0) then
    lambda2  =react_t_delay(1) 
    a_t = a_t + p_sum(M_react_channel)*(lambda2-lambda1)
     F = 1 - EXP(-a_t)
   end if
end do  
if (ongoing.gt.0) t = lambda1+( - log(1-r1)-a_t+p_sum(M_react_channel)*(lambda2-lambda1))/p_sum(M_react_channel)
if (ongoing.eq.0) t = lambda1+( - log(1-r1)-a_t)/p_sum(M_react_channel)
end if
  sum = r2*p_sum(M_react_channel)
x= p_sum(24)
   if ( sum .lt. x) then
   index = 1
     x =  tree_structure(index)%psum 
   else
    index = 2
     x =  tree_structure(index)%psum 
  end if
  do while (  tree_structure(index)%childgauche /= 0)
   if ( sum .lt. x) then
        index = tree_structure(index)%childgauche
        x = tree_structure(index)%psum
   else
        index = tree_structure(index)%childdroit
        x = tree_structure(index)%psum
   end if
  end do
 if  ( sum .lt. x) then
   mu = tree_structure(index)%label 
 else
   mu = tree_structure(index)%parents 
 end if 
  if (d_mu(mu).eq.0.0) then
  do l = 1, stochio_number(mu)
   number_species(l,mu)%species =   number_species(l,mu)%species+   update(l,mu) 
 end do
  mu_no = vector_mu_tilde(mu)
if (mu_no .le.24) then
if (mu_no .le.12) then
if (mu_no .le.6) then
if (mu_no .le.3) then
if (mu_no.eq.1) then
!Equation number=1
              a_nu(30)=c_nu(30)*Y(1)
              a_nu(36)=c_nu(36)*Y(1)*Y(20)
             x_min=30
else if (mu_no.eq.2) then
!Equation number=2
              a_nu(27)=c_nu(27)*Y(3)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
             x_min=27
else
!Equation number=3
              a_nu(3)=c_nu(3)*Y(19)
             x_min=3
end if
else
if (mu_no.eq.4) then
!Equation number=4
              a_nu(4)=c_nu(4)*Y(16)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
             x_min=4
else if (mu_no.eq.5) then
!Equation number=5
              a_nu(5)=c_nu(5)*Y(15)
              a_nu(6)=c_nu(6)*Y(15)
             x_min=5
else
!Equation number=8
              a_nu(8)=c_nu(8)*Y(13)
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
             x_min=8
end if
end if
else
if (mu_no .le.9) then
if (mu_no.eq.7) then
!Equation number=11
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
              a_nu(14)=c_nu(14)*Y(9)
              a_nu(32)=c_nu(32)*Y(9)*Y(17)
             x_min=10
else if (mu_no.eq.8) then
!Equation number=12
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
              a_nu(9)=c_nu(9)*Y(12)
             x_min=9
else
!Equation number=13
              a_nu(13)=c_nu(13)*Y(10)
              a_nu(14)=c_nu(14)*Y(9)
              a_nu(32)=c_nu(32)*Y(9)*Y(17)
             x_min=13
end if
else
if (mu_no.eq.10) then
!Equation number=14
              a_nu(14)=c_nu(14)*Y(9)
              a_nu(32)=c_nu(32)*Y(9)*Y(17)
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
             x_min=10
else if (mu_no.eq.11) then
!Equation number=15
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(14)=c_nu(14)*Y(9)
              a_nu(32)=c_nu(32)*Y(9)*Y(17)
             x_min=14
else
!Equation number=17
              a_nu(17)=c_nu(17)*Y(7)
              a_nu(18)=c_nu(18)*Y(7)
              a_nu(19)=c_nu(19)*Y(7)
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
              a_nu(4)=c_nu(4)*Y(16)
             x_min=4
end if
end if
end if
else
if (mu_no .le.18) then
if (mu_no .le.15) then
if (mu_no.eq.13) then
!Equation number=18
              a_nu(17)=c_nu(17)*Y(7)
              a_nu(18)=c_nu(18)*Y(7)
              a_nu(19)=c_nu(19)*Y(7)
             x_min=17
else if (mu_no.eq.14) then
!Equation number=19
              a_nu(17)=c_nu(17)*Y(7)
              a_nu(18)=c_nu(18)*Y(7)
              a_nu(19)=c_nu(19)*Y(7)
              a_nu(3)=c_nu(3)*Y(19)
             x_min=3
else
!Equation number=20
              a_nu(20)=c_nu(20)*Y(6)
              a_nu(21)=c_nu(21)*Y(6)
              a_nu(22)=c_nu(22)*Y(6)
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
              a_nu(14)=c_nu(14)*Y(9)
              a_nu(32)=c_nu(32)*Y(9)*Y(17)
             x_min=14
end if
else
if (mu_no.eq.16) then
!Equation number=21
              a_nu(20)=c_nu(20)*Y(6)
              a_nu(21)=c_nu(21)*Y(6)
              a_nu(22)=c_nu(22)*Y(6)
             x_min=20
else if (mu_no.eq.17) then
!Equation number=22
              a_nu(20)=c_nu(20)*Y(6)
              a_nu(21)=c_nu(21)*Y(6)
              a_nu(22)=c_nu(22)*Y(6)
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
             x_min=15
else
!Equation number=23
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
             x_min=23
end if
end if
else
if (mu_no .le.21) then
if (mu_no.eq.19) then
!Equation number=24
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
              a_nu(25)=c_nu(25)*Y(4)
              a_nu(26)=c_nu(26)*Y(4)
             x_min=23
else if (mu_no.eq.20) then
!Equation number=25
              a_nu(25)=c_nu(25)*Y(4)
              a_nu(26)=c_nu(26)*Y(4)
              a_nu(27)=c_nu(27)*Y(3)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
              a_nu(28)=c_nu(28)*Y(2)
              a_nu(29)=c_nu(29)*Y(2)
             x_min=25
else
!Equation number=26
              a_nu(25)=c_nu(25)*Y(4)
              a_nu(26)=c_nu(26)*Y(4)
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
             x_min=23
end if
else
if (mu_no.eq.22) then
!Equation number=27
              a_nu(27)=c_nu(27)*Y(3)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
             x_min=27
else if (mu_no.eq.23) then
!Equation number=28
              a_nu(28)=c_nu(28)*Y(2)
              a_nu(29)=c_nu(29)*Y(2)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
              a_nu(30)=c_nu(30)*Y(1)
              a_nu(36)=c_nu(36)*Y(1)*Y(20)
             x_min=28
else
!Equation number=29
              a_nu(28)=c_nu(28)*Y(2)
              a_nu(29)=c_nu(29)*Y(2)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
             x_min=28
end if
end if
end if
end if
else
if (mu_no .le.31) then
if (mu_no .le.30) then
if (mu_no .le.27) then
if (mu_no.eq.25) then
!Equation number=30
              a_nu(30)=c_nu(30)*Y(1)
              a_nu(36)=c_nu(36)*Y(1)*Y(20)
             x_min=30
else if (mu_no.eq.26) then
!Equation number=31
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
              a_nu(8)=c_nu(8)*Y(13)
             x_min=8
else
!Equation number=32
              a_nu(14)=c_nu(14)*Y(9)
              a_nu(32)=c_nu(32)*Y(9)*Y(17)
              a_nu(13)=c_nu(13)*Y(10)
             x_min=13
end if
else
if (mu_no.eq.28) then
!Equation number=33
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
              a_nu(20)=c_nu(20)*Y(6)
              a_nu(21)=c_nu(21)*Y(6)
              a_nu(22)=c_nu(22)*Y(6)
             x_min=15
else if (mu_no.eq.29) then
!Equation number=34
              a_nu(23)=c_nu(23)*Y(5)
              a_nu(24)=c_nu(24)*Y(5)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
              a_nu(4)=c_nu(4)*Y(16)
              a_nu(17)=c_nu(17)*Y(7)
              a_nu(18)=c_nu(18)*Y(7)
              a_nu(19)=c_nu(19)*Y(7)
             x_min=4
else
!Equation number=35
              a_nu(28)=c_nu(28)*Y(2)
              a_nu(29)=c_nu(29)*Y(2)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
              a_nu(27)=c_nu(27)*Y(3)
              a_nu(25)=c_nu(25)*Y(4)
              a_nu(26)=c_nu(26)*Y(4)
             x_min=25
end if
end if
else
if (mu_no .le.31) then
if (mu_no.eq.31) then
!Equation number=36
              a_nu(30)=c_nu(30)*Y(1)
              a_nu(36)=c_nu(36)*Y(1)*Y(20)
              a_nu(28)=c_nu(28)*Y(2)
              a_nu(29)=c_nu(29)*Y(2)
              a_nu(35)=c_nu(35)*Y(2)*Y(3)
             x_min=28
else if (mu_no.eq.32) then
else
end if
else
if (mu_no.eq.34) then
else if (mu_no.eq.35) then
else
end if
end if
end if
else
if (mu_no .le.31) then
if (mu_no .le.31) then
if (mu_no.eq.37) then
else if (mu_no.eq.38) then
else
end if
else
if (mu_no.eq.40) then
else if (mu_no.eq.41) then
else
end if
end if
else
if (mu_no .le.31) then
if (mu_no.eq.43) then
else if (mu_no.eq.44) then
else
end if
else
if (mu_no.eq.46) then
else if (mu_no.eq.47) then
else
end if
end if
end if
end if
end if






























































































































































































do j=x_min,M_react_channel 
      p_sum(j) = a_nu(j)+p_sum(j-1)
 end do
 else
  !ongoing = ongoing+1
   if (consuming(mu)) then
     do l = 1, stochio_number_tilde(mu)
       number_species_tilde(l,mu)%species =   number_species_tilde(l,mu)%species+   update_tilde(l,mu) 
     end do
     mu_yes = vector_mu_tilde(mu)
if (mu_yes .le.3) then
if (mu_yes.eq.1) then
!Equation number=6
              a_nu(5)=c_nu(5)*Y(15)
              a_nu(6)=c_nu(6)*Y(15)
              a_nu(4)=c_nu(4)*Y(16)
              a_nu(34)=c_nu(34)*Y(5)*Y(16)
             x_min=4
else if (mu_yes.eq.2) then
!Equation number=7
              a_nu(7)=c_nu(7)*Y(14)
              a_nu(5)=c_nu(5)*Y(15)
              a_nu(6)=c_nu(6)*Y(15)
             x_min=5
else
!Equation number=9
              a_nu(9)=c_nu(9)*Y(12)
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
             x_min=9
end if
else
if (mu_yes.eq.4) then
!Equation number=10
              a_nu(10)=c_nu(10)*Y(11)
              a_nu(11)=c_nu(11)*Y(11)
              a_nu(12)=c_nu(12)*Y(11)
              a_nu(31)=c_nu(31)*Y(11)*Y(18)
              a_nu(7)=c_nu(7)*Y(14)
             x_min=7
else if (mu_yes.eq.5) then
!Equation number=16
              a_nu(15)=c_nu(15)*Y(8)
              a_nu(16)=c_nu(16)*Y(8)
              a_nu(33)=c_nu(33)*Y(5)*Y(8)
              a_nu(9)=c_nu(9)*Y(12)
             x_min=9
else
end if
end if











































































   do j=x_min,M_react_channel 
      p_sum(j) = a_nu(j)+p_sum(j-1)
   end do
 end if
 ongoing = ongoing+1
  if (ongoing.gt.1) then
  k = 1
 check_t = react_t_delay(k)
   do while (t+ d_mu(mu) >= check_t ) 
      k = k+1
 check_t = react_t_delay(k)
   end do  
     react_t_delay(k+1:ongoing) =  react_t_delay(k:ongoing-1)
     mu_delayed(k+1:ongoing)   = mu_delayed(k:ongoing-1)
    react_t_delay(k) = t+ d_mu(mu)
    mu_delayed(k) = mu
  end if 
  if (ongoing.eq.1) then
  react_t_delay(1) = t+ d_mu(mu)
    mu_delayed(1) = mu
  end if  
  end if  
  if (n >=10.0) then
   matrix_nbre_transition(etat_initial,mu) = matrix_nbre_transition(etat_initial,mu)+1
   matrix_A_transition(etat_initial,mu) = matrix_A_transition(etat_initial,mu)+a_nu(mu)
   matrix_Adt_transition(etat_initial,mu) = matrix_Adt_transition(etat_initial,mu)+a_nu(mu)*tau
end if 
  etat_initial = mu
  if (mod(n,saving_frequency).eq.0) then
   i = i+1
   matrix_result(i,1:N_species) =  matrix_result(i,1:N_species) +Y(1:N_species)
    matrix_result(i,N_species+1) = matrix_result(i,N_species+1) + t
 end if
end do
end subroutine Gillespie_LDM








































































































































































































































































































































































































































































































































