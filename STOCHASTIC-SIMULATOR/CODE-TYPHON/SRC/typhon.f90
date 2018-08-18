!    Copyright (C) 2018 Quentin Vanhaelen

program  typhon


    use  variables
    use parameters
	use statistics
    
    
    implicit none
    real     :: start_time, stop_time
    integer ::  begin_step, end_step
    integer ::  i,j,ms
 
   
      write(*,*)'           ***********************************************************'
      write(*,*)'           ***********************************************************'
      write(*,*)'           *****             WELCOME IN TYPHON CODE           ********'
      write(*,*)'           *****  A PROGRAM TO PERFORM STOCHASTIC SIMULATIONS ********'
      write(*,*)'           ***** OF BIOLOGICAL PATHWAY WITH DELAYED REACTIONS ********'
      write(*,*)'           ***********************************************************'
      write(*,*)'           ***********************************************************'
      write(*,*)''
      write(*,*)''
	  CALL initialisation_setup()

	 

  last_node_tree = 2**(deep_tree+1)-2


  allocate (c_nu(M_react_channel), h_nu(M_react_channel), a_nu(M_react_channel),  p_sum(M_react_channel))
  allocate (Y(N_species ) )
  allocate (react_t_delay(M_delayed))
  allocate ( mu_delayed(M_delayed))
  allocate ( d_mu(M_react_channel))
  allocate (consuming(M_react_channel))
  allocate (update_matrix(N_species,  M_react_channel), update_matrix_tilde(N_species,  M_react_channel), update_matrix_product(N_species,  M_react_channel))
  allocate( matrix_result(number_events_saved, N_species+1),matrix_result_global(number_events_saved,N_species+1))
  allocate( binary_tree(last_node_tree))
  allocate(leftchild(last_node_tree),  rightchild(last_node_tree))
  allocate(  tree_structure(last_node_tree))
  allocate( number_species(3, M_react_channel),number_species_tilde(3, M_react_channel),number_species_product(3, M_react_channel))
  allocate( update(3,M_react_channel), update_tilde(3,M_react_channel), update_product(3,M_react_channel))
  allocate(stochio_number(M_react_channel),stochio_number_product(M_react_channel),stochio_number_tilde(M_react_channel))
  allocate(kinetic_rate(M_react_channel),vector_mu_tilde(M_react_channel))
  allocate(matrix_nbre_transition(M_react_channel,M_react_channel))
  allocate(matrix_A_transition(M_react_channel,M_react_channel))
  allocate(matrix_Adt_transition(M_react_channel,M_react_channel))
  allocate(matrix_nbre_transition_global(M_react_channel,M_react_channel))
  allocate(matrix_A_transition_global(M_react_channel,M_react_channel))
  allocate(matrix_Adt_transition_global(M_react_channel,M_react_channel))
  allocate(Y_concentration(N_species))

  Y_concentration(:) = 0.0
  CALL scaling_set_up()
  CALL stochiometric_initialization()
	
  matrix_result(:,:) = 0.0 
  matrix_nbre_transition(:,:) = 0
  matrix_A_transition(:,:) = 0.0
  matrix_Adt_transition(:,:) = 0.0
  
 
   call CPU_TIME(start_time)
    CALL init_random_seed() 
   do  nrun=1,  nruntotal
  !    CALL init_random_seed() 
	  CALL initial_condition()
     !  set the initial  time to zero
     t = 0.0
     !  call  the main routine to  compute the stochastic evolution of the system
     CALL Gillespie_LDM()
   
   if (mod(nrun,10).eq.0) then 
	  write(*,*)''
      write(*,*)'*****************************************************************************************************************************************'
      write(*,*)'  The run number', nrun,'   is finished. A number of',nrun,' runs of',ntotal,'reactions has been  done.'
      write(*,*)'                       Average physical duration of the simulation is t=',t, 'minutes'
      write(*,*)'*****************************************************************************************************************************************'
   end if
 
      
    if (mod(nrun,freq_run_saving).eq.0 .OR. nrun.eq. 1) then	
         matrix_result_global(:,:) = 0.0
		 matrix_result_global(:,:)= matrix_result(:,:)/real(nrun)
		 matrix_result_global(:,1:N_species)=matrix_result_global(:,1:N_species)/scaling
		 matrix_nbre_transition_global(:,:) = 0.0
	     matrix_A_transition_global(:,:) = 0.0
	     matrix_Adt_transition_global(:,:) = 0.0
		 matrix_nbre_transition_global(:,:) = matrix_nbre_transition(:,:)/real(nrun)
	     matrix_A_transition_global(:,:) = matrix_A_transition(:,:)/real(nrun)
	     matrix_Adt_transition_global(:,:) =  matrix_Adt_transition(:,:)/real(nrun)
		 do i=1,M_react_channel
		   do j=1,M_react_channel
		     if (  matrix_nbre_transition_global(i,j) > 0.0) then
			     matrix_A_transition_global(i,j) = matrix_A_transition_global(i,j)/matrix_nbre_transition_global(i,j) 
				 matrix_Adt_transition_global(i,j) =matrix_Adt_transition_global(i,j)/matrix_nbre_transition_global(i,j)
				 matrix_nbre_transition_global(i,j) =  matrix_nbre_transition_global(i,j)/real(ntotal)
			 end if
		   end do
		 end do
		 call  output_data()
 	end if
 
   
     
   end do ! fin du  do  sur les nrun
  
   call CPU_TIME(stop_time)
      SSA_duration = stop_time-start_time
	  SSA_physical_time_av =  matrix_result_global(number_events_saved,N_species+1)
     



  
 CALL output_statistics() 
  
 deallocate (c_nu, h_nu, a_nu,  p_sum)
 deallocate (Y)
 deallocate (update_matrix, update_matrix_tilde,update_matrix_product)
 deallocate (react_t_delay)
 deallocate ( d_mu)
 deallocate ( mu_delayed)
 deallocate(consuming)
 deallocate(matrix_result,matrix_result_global)
 deallocate( binary_tree)
 deallocate(leftchild,  rightchild)
 deallocate(  tree_structure)
 deallocate( number_species,number_species_product,number_species_tilde)
 deallocate(update,update_tilde,update_product)
 deallocate(stochio_number,stochio_number_tilde,stochio_number_product)

 deallocate(kinetic_rate,vector_mu_tilde)
 deallocate(matrix_nbre_transition)
 deallocate(matrix_A_transition)
 deallocate(matrix_Adt_transition)
 deallocate(matrix_nbre_transition_global)
 deallocate(matrix_A_transition_global)
 deallocate(matrix_Adt_transition_global)

 deallocate(Y_concentration)
 
end ! fin du programme typhon

!===================================================================================================

subroutine init_random_seed()
  
  
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
          
   CALL RANDOM_SEED(size = n)
   ALLOCATE(seed(n))
          
   CALL SYSTEM_CLOCK(COUNT=clock)
          
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)
          
   DEALLOCATE(seed)
end subroutine init_random_seed

!===================================================================================================



