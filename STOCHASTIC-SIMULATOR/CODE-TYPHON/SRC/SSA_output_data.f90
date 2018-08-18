!    Copyright (C) 2018 Quentin Vanhaelen

subroutine output_data()
use parameters
use variables
implicit none

character(len=70)  :: filename, filename2
character(len=4)  :: rnum
integer :: i ,j
write (rnum,'(i4.4)') nrun 

 filename = '../CODE-TYPHON-2/RESULTS-SSA/'//rnum//'-TEMPORAL-DYNAMICS'
 filename2 = '../CODE-TYPHON-2/RESULTS-SSA/'//rnum//'-TEMPORAL-DYNAMICS-LOW-RESOLUTION'
 open(1,file=filename, position='append') 
 open(2,file=filename2, position='append')  
do  i =1,number_events_saved
 write(1,1006) (  matrix_result_global(i,j), j=1,21)
  if (mod(i,10).eq.0) then
     write(2,1006) (  matrix_result_global(i,j),j=1,21)
  end if
 end do
close(1)
close(2)
 filename = '../CODE-TYPHON-2/RESULTS-SSA/'//rnum//'-M_nbre_transition'
 open(2,file=filename, position='append') 
do  i =1,M_react_channel
 write(2,1007) (  matrix_nbre_transition_global(i,j), j=1,36)
 end do
close(2)
 filename = '../CODE-TYPHON-2/RESULTS-SSA/'//rnum//'-M_A_transition'
 open(2,file=filename, position='append') 
do  i =1,M_react_channel
 write(2,1007) (  matrix_A_transition_global(i,j), j=1,36)
 end do
close(2)

 filename = '../CODE-TYPHON-2/RESULTS-SSA/'//rnum//'-M_Adt_transition'
 open(2,file=filename, position='append') 
do  i =1,M_react_channel
 write(2,1007) (  matrix_Adt_transition_global(i,j), j=1,36)
 end do
close(2)
1006  format(21(2x, e14.8))
1007  format(36(2x, e14.8))
end subroutine output_data






