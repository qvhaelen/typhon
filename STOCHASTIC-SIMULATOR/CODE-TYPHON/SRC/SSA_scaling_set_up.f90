!    Copyright (C) 2018 Quentin Vanhaelen

subroutine scaling_set_up()
use variables
use  parameters
implicit none
integer :: j,i,index_smallest
real :: intermediaire

!species name: LIFR
Y_concentration(1)=1000
!species name: LIFR-LIF
Y_concentration(2)=0
!species name: GP130
Y_concentration(3)=1000
!species name: LIFR-LIF-GP130-JAK
Y_concentration(4)=0
!species name: LIFR-LIF-GP130-JAK-PP
Y_concentration(5)=0
!species name: LIFR-LIF-GP130-JAK-PP-2STAT3
Y_concentration(6)=0
!species name: LIFR-LIF-GP130-JAK-PP-SOCS3
Y_concentration(7)=0
!species name: 2STAT3 
Y_concentration(8)=1000
!species name: p2STAT3
Y_concentration(9)=0
!species name: p2STAT3-PIAS3 
Y_concentration(10)=0
!species name: p2STAT3 
Y_concentration(11)=0
!species name: 2STAT3
Y_concentration(12)=0
!species name: p2STAT3-PIAS3 
Y_concentration(13)=0
!species name: SOCS3-mRNA 
Y_concentration(14)=0
!species name: SOCS3-mRNA 
Y_concentration(15)=0
!species name: SOCS3 
Y_concentration(16)=0
!species name: PIAS3 
Y_concentration(17)=0
!species name: PIAS3 
Y_concentration(18)=0
!species name: LIFR-LIF-GP130-JAK-SOCS3
Y_concentration(19)=0
!species name: LIF INDUCER
Y_concentration(20)=0
intermediaire = 100.0
do i=1,N_species
  if ((Y_concentration(i) <= intermediaire ) .AND.  (Y_concentration(i)> 0.0) ) then
     intermediaire = Y_concentration(i) 
     index_smallest = i
   end if
end do
scaling = 0.0
do while ( floor(intermediaire) < 100 )
      scaling = scaling+10.0
      intermediaire = scaling*Y_concentration( index_smallest)
end do
if (scaling == 0.0) then
  scaling = 1.0
end if
end subroutine scaling_set_up







