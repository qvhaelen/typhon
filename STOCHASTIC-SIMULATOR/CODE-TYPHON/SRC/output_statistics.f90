!    Copyright (C) 2018 Quentin Vanhaelen

subroutine output_statistics()


use  statistics
use  parameters
implicit none
character(len=70)  :: filename
character(8) :: date
character(10) :: time
character(5) ::  zone
integer, dimension(8) ::  value
CALL date_and_time(date,time,zone,value)
filename = '../CODE-TYPHON-2/RESULTS-SSA/SSA_statistics'  
  open(1, file = filename, position = 'append')
   write(1,*) '***************  Date of the simulation completion*********************************'  
   write(1,*) '  Year = ', value(1)  
   write(1,*) '  Month = ', value(2)  
   write(1,*) '  Day of the month = ', value(3)  
   write(1,*) '  Time difference with UTC in minutes = ', value(4)  
   write(1,*) '  Hour of the day = ', value(5)  
   write(1,*) '  The minutes of the hour = ', value(6)  
   write(1,*) '  The seconds of the minute = ', value(7) 
   write(1,*) '***************Statistics for SSA simulation **************************************'  
   write(1,*) 'Value of the scaling factor used = ', scaling  
   write(1,*) '  A number of',nruntotal,' run(s) of',ntotal,'reactions has been  done.'  
   write(1,*) '  Average physical duration of one run is t = ',int(SSA_physical_time_av), 'minutes'  
   write(1,*) '  Total  number of fired reactions is ',nruntotal* ntotal  
   write(1,*) '  The Total Elapsed CPU time for the complete set of runs is ', int(SSA_duration), 'seconds'  
   write(1,*) '  Approximately',SSA_duration/60.,'minutes'  
   write(1,*)'   Average speed of the SSA is ',  floor(nruntotal* ntotal/SSA_duration) , ' reactions/s'  
   write(1,*)'************************************************************************************'  
 close(1) 

end subroutine output_statistics


































































