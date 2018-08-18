!    Copyright (C) 2018 Quentin Vanhaelen


subroutine initialisation_setup()
use parameters

implicit none


 ntotal = 3000000   ! Total  number of chemical  reactions to  be fired for each run
 nruntotal = 200 !  Number of run

 ! Saving frequecncy for the output files. freq_run_saving must be a dvisor of nruntotal 
 freq_run_saving = 50 
 !  saving_frequency is computed automatically in a way to  the system will save exactly 50 000 time points of the simulation in full resolution and 5000 otherwise
 saving_frequency = floor(ntotal/(number_events_saved*1.0))
 
 


end subroutine initialisation_setup 