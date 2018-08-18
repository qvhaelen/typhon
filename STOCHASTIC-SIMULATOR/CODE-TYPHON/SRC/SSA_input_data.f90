!    Copyright (C) 2018 Quentin Vanhaelen

subroutine  initial_condition()
use variables
use  parameters
implicit none
integer :: j,i
!  initial  amount for each chemical  species
!species name: LIFR
Y(1)=int(scaling*Y_concentration(1))
!species name: LIFR-LIF
Y(2)=int(scaling*Y_concentration(2))
!species name: GP130
Y(3)=int(scaling*Y_concentration(3))
!species name: LIFR-LIF-GP130-JAK
Y(4)=int(scaling*Y_concentration(4))
!species name: LIFR-LIF-GP130-JAK-PP
Y(5)=int(scaling*Y_concentration(5))
!species name: LIFR-LIF-GP130-JAK-PP-2STAT3
Y(6)=int(scaling*Y_concentration(6))
!species name: LIFR-LIF-GP130-JAK-PP-SOCS3
Y(7)=int(scaling*Y_concentration(7))
!species name: 2STAT3 
Y(8)=int(scaling*Y_concentration(8))
!species name: p2STAT3
Y(9)=int(scaling*Y_concentration(9))
!species name: p2STAT3-PIAS3 
Y(10)=int(scaling*Y_concentration(10))
!species name: p2STAT3 
Y(11)=int(scaling*Y_concentration(11))
!species name: 2STAT3
Y(12)=int(scaling*Y_concentration(12))
!species name: p2STAT3-PIAS3 
Y(13)=int(scaling*Y_concentration(13))
!species name: SOCS3-mRNA 
Y(14)=int(scaling*Y_concentration(14))
!species name: SOCS3-mRNA 
Y(15)=int(scaling*Y_concentration(15))
!species name: SOCS3 
Y(16)=int(scaling*Y_concentration(16))
!species name: PIAS3 
Y(17)=int(scaling*Y_concentration(17))
!species name: PIAS3 
Y(18)=int(scaling*Y_concentration(18))
!species name: LIFR-LIF-GP130-JAK-SOCS3
Y(19)=int(scaling*Y_concentration(19))
!species name: LIF INDUCER
Y(20)=int(scaling*Y_concentration(20))

h_nu( 1 )=1.0
h_nu( 2 )=1.0
h_nu(3)=Y(19)
h_nu(4)=Y(16)
h_nu(5)=Y(15)
h_nu(6)=Y(15)
h_nu(7)=Y(14)
h_nu(8)=Y(13)
h_nu(9)=Y(12)
h_nu(10)=Y(11)
h_nu(11)=Y(11)
h_nu(12)=Y(11)
h_nu(13)=Y(10)
h_nu(14)=Y(9)
h_nu(15)=Y(8)
h_nu(16)=Y(8)
h_nu(17)=Y(7)
h_nu(18)=Y(7)
h_nu(19)=Y(7)
h_nu(20)=Y(6)
h_nu(21)=Y(6)
h_nu(22)=Y(6)
h_nu(23)=Y(5)
h_nu(24)=Y(5)
h_nu(25)=Y(4)
h_nu(26)=Y(4)
h_nu(27)=Y(3)
h_nu(28)=Y(2)
h_nu(29)=Y(2)
h_nu(30)=Y(1)
h_nu(31)=Y(11)*Y(18)
h_nu(32)=Y(9)*Y(17)
h_nu(33)=Y(5)*Y(8)
h_nu(34)=Y(5)*Y(16)
h_nu(35)=Y(2)*Y(3)
h_nu(36)=Y(1)*Y(20)
p_sum(1) =  h_nu(1)*c_nu(1)
a_nu(1) =  p_sum(1)
 do j=2,M_react_channel
      a_nu(j) = h_nu(j)*c_nu(j)
     p_sum(j) = a_nu(j)+p_sum(j-1)
 end do
end subroutine  initial_condition
subroutine  stochiometric_initialization()
use variables
use parameters
implicit none
integer :: j,i,k
kinetic_rate(1)=60
kinetic_rate(2)=50
kinetic_rate(3)=0.0099
kinetic_rate(4)=30
kinetic_rate(5)=30
kinetic_rate(6)=30
kinetic_rate(7)=30
kinetic_rate(8)=30
kinetic_rate(9)=6/scaling
kinetic_rate(10)=0.0099
kinetic_rate(11)=6
kinetic_rate(12)=0.0099
kinetic_rate(13)=5
kinetic_rate(14)=3/scaling
kinetic_rate(15)=0.6/scaling
kinetic_rate(16)=0.3/scaling
kinetic_rate(17)=0.3/scaling
kinetic_rate(18)=0.2/scaling
kinetic_rate(19)=0.1
kinetic_rate(20)=0.04
kinetic_rate(21)=0.0263
kinetic_rate(22)=0.0099
kinetic_rate(23)=0.0099
kinetic_rate(24)=0.0099
kinetic_rate(25)=9.3*scaling
kinetic_rate(26)=10.23*scaling
kinetic_rate(27)=0.007
kinetic_rate(28)=0.0099
kinetic_rate(29)=0.0011
kinetic_rate(30)=0.63
kinetic_rate(31)=0.01
kinetic_rate(32)=2
kinetic_rate(33)=2
kinetic_rate(34)=2
kinetic_rate(35)=0.1
kinetic_rate(36)=0.28
update_matrix(1,1)=1
update_matrix(3,2)=1
update_matrix(19,3)=-1
update_matrix(16,4)=-1
update_matrix(15,5)=-1
update_matrix(15,6)=-1
update_matrix(15,6)=1
update_matrix(16,6)=1
update_matrix(14,7)=-1
update_matrix(15,7)=1
update_matrix(13,8)=-1
update_matrix(11,8)=1
update_matrix(18,8)=1
update_matrix(12,9)=-1
update_matrix(8,9)=1
update_matrix(11,10)=-1
update_matrix(11,10)=1
update_matrix(14,10)=1
update_matrix(11,11)=-1
update_matrix(9,11)=1
update_matrix(11,12)=-1
update_matrix(12,12)=1
update_matrix(10,13)=-1
update_matrix(9,13)=1
update_matrix(17,13)=1
update_matrix(9,14)=-1
update_matrix(11,14)=1
update_matrix(8,15)=-1
update_matrix(9,15)=1
update_matrix(8,16)=-1
update_matrix(12,16)=1
update_matrix(7,17)=-1
update_matrix(5,17)=1
update_matrix(16,17)=1
update_matrix(7,18)=-1
update_matrix(7,19)=-1
update_matrix(19,19)=1
update_matrix(6,20)=-1
update_matrix(5,20)=1
update_matrix(9,20)=1
update_matrix(6,21)=-1
update_matrix(6,22)=-1
update_matrix(5,22)=1
update_matrix(8,22)=1
update_matrix(5,23)=-1
update_matrix(5,24)=-1
update_matrix(4,24)=1
update_matrix(4,25)=-1
update_matrix(3,25)=1
update_matrix(2,25)=1
update_matrix(4,26)=-1
update_matrix(5,26)=1
update_matrix(3,27)=-1
update_matrix(2,28)=-1
update_matrix(1,28)=1
update_matrix(20,28)=1
update_matrix(2,29)=-1
update_matrix(1,30)=-1
update_matrix(11,31)=-1
update_matrix(18,31)=-1
update_matrix(13,31)=1
update_matrix(9,32)=-1
update_matrix(17,32)=-1
update_matrix(10,32)=1
update_matrix(5,33)=-1
update_matrix(8,33)=-1
update_matrix(6,33)=1
update_matrix(5,34)=-1
update_matrix(16,34)=-1
update_matrix(7,34)=1
update_matrix(2,35)=-1
update_matrix(3,35)=-1
update_matrix(4,35)=1
update_matrix(1,36)=-1
update_matrix(20,36)=-1
update_matrix(2,36)=1
update_matrix(15,6)=0
update_matrix(11,10)=0






























do i=1, M_react_channel
k = 0
  do j=1,N_species
   if (update_matrix(j,i) /= 0) then
    k = k+1
    update(k,i) = update_matrix(j,i)
     number_species(k,i)%species => Y(j)
    end if
  end do
 stochio_number(i) = k
  end do
do i=1,31
  leftchild(i) = 2*i+1
  rightchild(i) = 2*i+2
  tree_structure(i)%childgauche =  leftchild(i) 
 tree_structure(i)%childdroit =  rightchild(i) 
end do
do i=31,62
  leftchild(i) = 0
  rightchild(i) = 0
  tree_structure(i)%childgauche =  leftchild(i) 
 tree_structure(i)%childdroit =  rightchild(i)
end do

!level 1 of the binary tree
          binary_tree(1)=12
          binary_tree(2)=36
!level 2 of the binary tree
          binary_tree(3)=6
          binary_tree(4)=18
          binary_tree(5)=30
          binary_tree(6)=36
!level 3 of the binary tree
          binary_tree(7)=3
          binary_tree(8)=9
          binary_tree(9)=15
          binary_tree(10)=21
          binary_tree(11)=27
          binary_tree(12)=33
          binary_tree(13)=36
          binary_tree(14)=36
!level 4 of the binary tree
          binary_tree(15)=2
          binary_tree(16)=5
          binary_tree(17)=8
          binary_tree(18)=11
          binary_tree(19)=14
          binary_tree(20)=17
          binary_tree(21)=20
          binary_tree(22)=23
          binary_tree(23)=26
          binary_tree(24)=29
          binary_tree(25)=32
          binary_tree(26)=35
          binary_tree(27)=36
          binary_tree(28)=36
          binary_tree(29)=36
          binary_tree(30)=36
!level 5 of the binary tree
          binary_tree(31)=1
          binary_tree(32)=3
          binary_tree(33)=4
          binary_tree(34)=6
          binary_tree(35)=7
          binary_tree(36)=9
          binary_tree(37)=10
          binary_tree(38)=12
          binary_tree(39)=13
          binary_tree(40)=15
          binary_tree(41)=16
          binary_tree(42)=18
          binary_tree(43)=19
          binary_tree(44)=21
          binary_tree(45)=22
          binary_tree(46)=24
          binary_tree(47)=25
          binary_tree(48)=27
          binary_tree(49)=28
          binary_tree(50)=30
          binary_tree(51)=31
          binary_tree(52)=33
          binary_tree(53)=34
          binary_tree(54)=36
          binary_tree(55)=36
          binary_tree(56)=36
          binary_tree(57)=36
          binary_tree(58)=36
          binary_tree(59)=36
          binary_tree(60)=36
          binary_tree(61)=36
          binary_tree(62)=36
do i=1,62
  tree_structure(i)%parents =binary_tree(floor((i-1)/2.0)) 
   tree_structure(i)%label =binary_tree(i)
   tree_structure(i)%psum => p_sum( binary_tree(i) )
end do
  do i = 1, M_delayed
   react_t_delay(i) = 10000.0
    mu_delayed(i) = 0
   end do
   do i=1,N_species
   do  j=1, M_react_channel
    if (update_matrix(i,j).le. 0) then
            update_matrix_tilde(i,j) =  update_matrix(i,j)
             update_matrix_product(i,j) =  0
      end if  
    if (update_matrix(i,j).gt. 0) then 
             update_matrix_tilde(i,j) =  0
             update_matrix_product(i,j) =  update_matrix(i,j)
     end if 
   end do
   end do
do i=1, M_react_channel
k = 0
 do j=1,N_species
    if (update_matrix_tilde(j,i) /= 0) then
      k = k+1
     update_tilde(k,i) = update_matrix_tilde(j,i)
     number_species_tilde(k,i)%species => Y(j)
    end if
  end do
 stochio_number_tilde(i) = k
end do
do i=1, M_react_channel
k = 0
  do j=1,N_species
   if (update_matrix_product(j,i) /= 0) then
     k = k+1
      update_product(k,i) = update_matrix_product(j,i)
    number_species_product(k,i)%species => Y(j)
   end if
  end do
  stochio_number_product(i) = k
end do
d_mu(1)= 0*1.0
d_mu(2)= 0*1.0
d_mu(3)= 0*1.0
d_mu(4)= 0*1.0
d_mu(5)= 0*1.0
d_mu(6)= 100*1.0 !d1-d4 = 0.00001! 100*1.0 
d_mu(7)= 1.9215*1.0!d1-d4 = 0.00001!0.00001! 1.9215*1.0
d_mu(8)= 0*1.0
d_mu(9)= 1.4*1.0!d1-d4 = 0.00001!0.00001!1.4*1.0
d_mu(10)=  1.5875*1.0  !d4 = 1.5875*1.0!d3 = 0.1 !d2 = 0.001 !d1 = 0.00001!transcription 1.5875*1.0
d_mu(11)= 0*1.0
d_mu(12)= 0*1.0
d_mu(13)= 0*1.0
d_mu(14)= 0*1.0
d_mu(15)= 0*1.0
d_mu(16)= 1.2*1.0!d1-d4 = 0.00001! 0.00001!1.2*1.0
d_mu(17)= 0*1.0
d_mu(18)= 0*1.0
d_mu(19)= 0*1.0
d_mu(20)= 0*1.0
d_mu(21)= 0*1.0
d_mu(22)= 0*1.0
d_mu(23)= 0*1.0
d_mu(24)= 0*1.0
d_mu(25)= 0*1.0
d_mu(26)= 0*1.0
d_mu(27)= 0*1.0
d_mu(28)= 0*1.0
d_mu(29)= 0*1.0
d_mu(30)= 0*1.0
d_mu(31)= 0*1.0
d_mu(32)= 0*1.0
d_mu(33)= 0*1.0
d_mu(34)= 0*1.0
d_mu(35)= 0*1.0
d_mu(36)= 0*1.0
consuming(1)= .TRUE.
consuming(2)= .TRUE.
consuming(3)= .TRUE.
consuming(4)= .TRUE.
consuming(5)= .TRUE.
consuming(6)= .FALSE.
consuming(7)= .TRUE.
consuming(8)= .TRUE.
consuming(9)= .TRUE.
consuming(10)= .FALSE.
consuming(11)= .TRUE.
consuming(12)= .TRUE.
consuming(13)= .TRUE.
consuming(14)= .TRUE.
consuming(15)= .TRUE.
consuming(16)= .TRUE.
consuming(17)= .TRUE.
consuming(18)= .TRUE.
consuming(19)= .TRUE.
consuming(20)= .TRUE.
consuming(21)= .TRUE.
consuming(22)= .TRUE.
consuming(23)= .TRUE.
consuming(24)= .TRUE.
consuming(25)= .TRUE.
consuming(26)= .TRUE.
consuming(27)= .TRUE.
consuming(28)= .TRUE.
consuming(29)= .TRUE.
consuming(30)= .TRUE.
consuming(31)= .TRUE.
consuming(32)= .TRUE.
consuming(33)= .TRUE.
consuming(34)= .TRUE.
consuming(35)= .TRUE.
consuming(36)= .TRUE.
c_nu(1)= kinetic_rate(25)
c_nu(2)= kinetic_rate(26)
c_nu(3)= kinetic_rate(24)
c_nu(4)= kinetic_rate(27)
c_nu(5)= kinetic_rate(19)
c_nu(6)= kinetic_rate(31)
c_nu(7)= kinetic_rate(36)
c_nu(8)= kinetic_rate(8)
c_nu(9)= kinetic_rate(33)
c_nu(10)= kinetic_rate(30)
c_nu(11)= kinetic_rate(35)
c_nu(12)= kinetic_rate(20)
c_nu(13)= kinetic_rate(7)
c_nu(14)= kinetic_rate(34)
c_nu(15)= kinetic_rate(21)
c_nu(16)= kinetic_rate(32)
c_nu(17)= kinetic_rate(6)
c_nu(18)= kinetic_rate(10)
c_nu(19)= kinetic_rate(11)
c_nu(20)= kinetic_rate(5)
c_nu(21)= kinetic_rate(3)
c_nu(22)= kinetic_rate(4)
c_nu(23)= kinetic_rate(12)
c_nu(24)= kinetic_rate(13)
c_nu(25)= kinetic_rate(1)
c_nu(26)= kinetic_rate(2)
c_nu(27)= kinetic_rate(23)
c_nu(28)= kinetic_rate(29)
c_nu(29)= kinetic_rate(28)
c_nu(30)= kinetic_rate(22)
c_nu(31)= kinetic_rate(17)
c_nu(32)= kinetic_rate(16)
c_nu(33)= kinetic_rate(14)
c_nu(34)= kinetic_rate(15)
c_nu(35)= kinetic_rate(9)
c_nu(36)= kinetic_rate(18)
vector_mu_tilde(1)=1
vector_mu_tilde(2)=2
vector_mu_tilde(3)=3
vector_mu_tilde(4)=4
vector_mu_tilde(5)=5
vector_mu_tilde(6)=1
vector_mu_tilde(7)=2
vector_mu_tilde(8)=6
vector_mu_tilde(9)=3
vector_mu_tilde(10)=4
vector_mu_tilde(11)=7
vector_mu_tilde(12)=8
vector_mu_tilde(13)=9
vector_mu_tilde(14)=10
vector_mu_tilde(15)=11
vector_mu_tilde(16)=5
vector_mu_tilde(17)=12
vector_mu_tilde(18)=13
vector_mu_tilde(19)=14
vector_mu_tilde(20)=15
vector_mu_tilde(21)=16
vector_mu_tilde(22)=17
vector_mu_tilde(23)=18
vector_mu_tilde(24)=19
vector_mu_tilde(25)=20
vector_mu_tilde(26)=21
vector_mu_tilde(27)=22
vector_mu_tilde(28)=23
vector_mu_tilde(29)=24
vector_mu_tilde(30)=25
vector_mu_tilde(31)=26
vector_mu_tilde(32)=27
vector_mu_tilde(33)=28
vector_mu_tilde(34)=29
vector_mu_tilde(35)=30
vector_mu_tilde(36)=31
end subroutine  stochiometric_initialization























