subroutine ImEx2fp

  use arrays
  use global_numbers
  use lu  

  implicit none

  integer k,j, i,l,m, p
  real(kind=8) dt_temp, igamma, dum
  integer d, code, salto
  real(kind=8) NormaerrorSR, NormaerrorER, Norma_beforeSR, Norma_beforeER
  real(kind=8) NormaerrorE, NormaerrorM, Norma_beforeM, Norma_beforeE, Norma_beforePg
  real (kind=8) tol_down, NormaerrorPg
  igamma = 1.0d0 - 1.0d0/sqrt(2.0d0)

  tol_down = 1.d-15

  u_pp(:,:) = 0.0d0
  u_p = 0.0d0
  u_p(:,:) = u(:,:)  


  do m=1,3


     
     call rhs

     if(m.eq.1) then

        do p=1,1000
           

           ErrorER(:) = 0.0d0
           ErrorE(:) = 0.0d0
           ErrorM(:) = 0.0d0
           ErrorSR(:) = 0.0d0
           ErrorPg(:) = 0.0d0 
           !**********************************************************
           !**********************************************************
           M_0(:,:) = 0.0d0
           u_0(:,:) = 0.0d0
           M_f(:,:) = 0.0d0
           b(:) = 0.0d0
           chi_t(:) = opacity_a * w(1,:) !pii * cte_rad * Temp**4 *opacity_a   
           chi_s(:) = 0.0d0  !opacity_s * w(1,:)
           Temp(:) = w(3,:)/w(1,:)
           !**********************************************************

           do i=0,Nr


              M_f(1,1) = (chi_t(i) + chi_s(i))/LF(i) + LF(i) * w(2,i)**2 &
                   * (chi_t(i) * (2.0d0/(1.0d0 + 2.0d0 * LF(i)**2) -1.0d0) + chi_s(i) &
                   * (2.0d0 - 1.0d0/(1.0d0 + 2.0d0 * LF(i)**2)))
              M_f(1,2) =  LF(i) * w(2,i) * (chi_t(i) *(1.0d0 - (4.0d0/(1.0d0 + 2.0d0 * LF(i)**2))) + 2.0d0 * chi_s(i) &
                   * (1.0d0/(1.0d0 + 2.0d0 * LF(i)**2) - 1.0d0))
              M_f(2,1) =  -LF(i)* w(2,i) * (chi_t(i) + chi_s(i) &
                   * (3.0d0/(1.0d0 + 2.0d0 * LF(i)**2) - 2.0d0))
              M_f(2,2) =  -LF(i)*(2.0d0 * chi_s(i) * (1.0d0 - 3.0d0/(1.0d0 + 2.0d0 * LF(i)**2)) - chi_t(i))

              !***********************

              b(1) =  - chi_t(i) * cte_rad * Temp(i)**4 * LF(i) * w(2,i)
              b(2) =  - LF(i)*chi_t(i) * cte_rad * Temp(i)**4

              !***********************
              I_r(1,1) = 1.0d0
              I_r(1,2) = 0.0d0
              I_r(2,1) = 0.0d0
              I_r(2,2) = 1.0d0
              u_Ir(:,1) = 0.0d0  

              !***************************+ solves the system u^(i) = u^(') + beta * (M_tilde*u^(i)+b) *********************
              !*********** M_tilde comes from the linealization of M^(i) given by M^(i) = M_tilde*u^(i) +b *****************
              !******** and we obtain the solution for the radiation energy and radiation momentum for the first stage******
              !*************************************************************************************************************
              do j=1,2

                 u_Ir(j,1) = u_p(j+3,i)  + dt*igamma*(-b(j))

                 do k=1,2

                    I_r(j,k) = I_r(j,k) - dt*igamma*(-M_f(j,k))

                 end do

              end do

!!$           call Gaussjordan(I_r,2,u_Ir,1)           
              call LUDCMP(I_r,2,INDX,D,CODE)
              call LUBKSB(I_r,2,INDX,u_Ir)

              do l=1,2 ! se construye U0
                 u_0(l+3,i) = u_Ir(l,1)             
              end do

              !*************************************************************************************************************
              !*************************************************************************************************************


              !******************** We obtain the implicit terms for the 1st stage M^(i) for the radiatio ******************
              !*************************************************************************************************************              
              do j=1,2
                 do k=1,2
                    M_0(j,i) = M_0(j,i) - M_f(j,k)*u_0(k+3,i)
                 end do
                 M_0(j,i) = M_0(j,i)-b(j)
              end do
           end do

           !********** Hidrodynamic solution for the 1st supposing M_rad = -M_fluid **********************************
           !**********************************************************************************************************           

           u_0(1,:)=u_p(1,:)

           do j=1,2

              u_0(j+1,:) = u_p(j+1,:) - dt*igamma*(M_0(j,:))

           end do
           !************************************************************************************************************           


           !**************** Error norms of Total energy, momentum, and radiation energy and radiation momentum ********
           !****************                             for iterations                           **********************


           if (p.gt.1) then

              ErrorE(:) = (u_0(3,:)+u_0(5,:))-(u(3,:)+u(5,:))
              ErrorM(:) = (u_0(2,:)+u_0(4,:))-(u(2,:)+u(4,:))
              ErrorER(:) = u_0(5,:)-u(5,:)
              ErrorSR(:) = u_0(4,:)-u(4,:)
              ErrorPg(:) = w(3,:) - pg_before(:)

              NormaerrorE = 0.0d0
              NormaerrorM = 0.0d0
              NormaerrorER = 0.0d0
              NormaerrorSR = 0.0d0
              NormaerrorPg = 0.0d0

              do i=0,Nr
                 NormaerrorER = NormaerrorER +abs(ErrorER(i))                
                 NormaerrorSR = NormaerrorSR +abs(ErrorSR(i))
                 NormaerrorPg = NormaerrorPg +abs(ErrorPg(i))                  
              end do

              do i=0,Nr-1
                 NormaerrorE = NormaerrorE + 0.5d0*(abs(ErrorE(i+1))+abs(ErrorE(i)))*dr                
                 NormaerrorM = NormaerrorM + 0.5d0*(abs(ErrorM(i+1))+abs(ErrorM(i)))*dr           
              end do
!!$              print *, '    Norma ER ', NormaerrorER, '    Norma SR ', NormaerrorSR            
!!$              print *, '    Norma E ', NormaerrorE, '    Norma M ', NormaerrorM
!!$              print *, '    Norma Pg ', NormaerrorPg              
!!$              print *, 'convergenciaER', NormaerrorER - Norma_beforeER
!!$              print *, 'convergenciaSR', NormaerrorSR - Norma_beforeSR
!!$              print *, 'convergenciaPg', NormaerrorPg - Norma_beforePg                
!!$              print *, '**********************************************************'

              if(abs(NormaerrorER - Norma_beforeER).lt.tol_down.or.abs(NormaerrorSR - Norma_beforeSR).lt.tol_down.or. &
                   abs(NormaerrorPg - Norma_beforePg).lt.tol_down ) then

!!$                 print*, "converge al orden", abs(NormaerrorER - Norma_beforeER), p
!!$                 print*, "converge al orden", abs(NormaerrorSR - Norma_beforeSR)
!!$                 print*, "converge al orden", abs(NormaerrorPg - Norma_beforePg)
!                 print*, 'converge1', p
                 
                 u(:,:) = 0.0d0

                 u(:,:) = u_0(:,:)
                 
                 go to 1
              end if
           end if            

           u(:,:) = 0.0d0

           u(:,:) = u_0(:,:)

           Norma_beforeER = NormaerrorER
           Norma_beforeSR = NormaerrorSR
           Norma_beforePg = NormaerrorPg           
           pg_before(:) = w(3,:)

           call BC
           call recv_primitives

        end do

1       salto=0.0d0

        !******************************************************
        !******************************************************
     else if (m.eq.2) then

        u_pp(:,:) = (3.0d0*igamma - 1.0d0)/igamma *u_p(:,:) &
             + (1.0d0 - 2.0d0*igamma)/igamma * u_0(:,:) + dt * rhs_u(:,:)

        do p=1,1000

           

           ErrorER(:) = 0.0d0
           ErrorE(:) = 0.0d0
           ErrorM(:) = 0.0d0
           ErrorSR(:) = 0.0d0
           ErrorPg(:) = 0.0d0           

           M_1(:,:) = 0.0d0
           u_1(:,:) = 0.0d0 
           M_f(:,:) = 0.0d0
           b(:) = 0.0d0
           chi_t(:) = opacity_a * w(1,:) !pii * cte_rad * Temp**4 *opacity_a   
           chi_s(:) = 0.0d0  !opacity_s !* w(1,:)
           Temp(:) = w(3,:)/w(1,:)            




           do i=0,Nr             


              M_f(1,1) = (chi_t(i) + chi_s(i))/LF(i) + LF(i) * w(2,i)**2 &
                   * (chi_t(i) * (2.0d0/(1.0d0 + 2.0d0 * LF(i)**2) -1.0d0) + chi_s(i) &
                   * (2.0d0 - 1.0d0/(1.0d0 + 2.0d0 * LF(i)**2)))
              M_f(1,2) =  LF(i) * w(2,i) * (chi_t(i) *(1.0d0 - (4.0d0/(1.0d0 + 2.0d0 * LF(i)**2))) + 2.0d0 * chi_s(i) &
                   * (1.0d0/(1.0d0 + 2.0d0 * LF(i)**2) - 1.0d0))
              M_f(2,1) =  -LF(i)* w(2,i) * (chi_t(i) + chi_s(i) &
                   * (3.0d0/(1.0d0 + 2.0d0 * LF(i)**2) - 2.0d0))
              M_f(2,2) =  -LF(i)*(2.0d0 * chi_s(i) * (1.0d0 - 3.0d0/(1.0d0 + 2.0d0 * LF(i)**2)) - chi_t(i))

              !***********************

              b(1) =  - chi_t(i) * cte_rad * Temp(i)**4 * LF(i) * w(2,i)
              b(2) =  - LF(i)*chi_t(i) * cte_rad * Temp(i)**4

              !***********************
              I_r(1,1) = 1.0d0
              I_r(1,2) = 0.0d0
              I_r(2,1) = 0.0d0
              I_r(2,2) = 1.0d0

              u_Ir(:,1) = 0.0d0           


              do j=1,2

                 u_Ir(j,1) = u_pp(j+3,i)  + dt*igamma*(-b(j))
                 do k=1,2
                    I_r(j,k) = I_r(j,k) - dt*igamma*(-M_f(j,k))
                 end do
              end do

              !           call Gaussjordan(I_r,2,u_Ir,1)
              call LUDCMP(I_r,2,INDX,D,CODE)
              call LUBKSB(I_r,2,INDX,u_Ir)           

              do l=1,2 ! se construye U1

                 u_1(l+3,i) = u_Ir(l,1)              

              end do



              do j=1,2
                 do k=1,2

                    M_1(j,i) = M_1(j,i) - M_f(j,k)*u_1(k+3,i)

                 end do

                 M_1(j,i) = M_1(j,i)-b(j)

              end do


           end do


           u_1(1,:) = u_pp(1,:)

           do l=1,2

              u_1(l+1,:) = u_pp(l+1,:) - dt*igamma*(M_1(l,:))

           end do

           if (p.gt.1) then

              ErrorE(:) = (u_1(3,:)+u_0(5,:))-(u(3,:)+u(5,:))
              ErrorM(:) = (u_1(2,:)+u_0(4,:))-(u(2,:)+u(4,:))
              ErrorER(:) = u_1(5,:)-u(5,:)
              ErrorSR(:) = u_1(4,:)-u(4,:)
              ErrorPg(:) = w(3,:) - pg_before(:)

              NormaerrorE = 0.0d0
              NormaerrorM = 0.0d0
              NormaerrorER = 0.0d0
              NormaerrorSR = 0.0d0
              NormaerrorPg = 0.0d0

              do i=0,Nr
                 NormaerrorER = NormaerrorER +abs(ErrorER(i))                
                 NormaerrorSR = NormaerrorSR +abs(ErrorSR(i))
                 NormaerrorPg = NormaerrorPg +abs(ErrorPg(i))                  
              end do

              do i=0,Nr-1
                 NormaerrorE = NormaerrorE + 0.5d0*(abs(ErrorE(i+1))+abs(ErrorE(i)))*dr                
                 NormaerrorM = NormaerrorM + 0.5d0*(abs(ErrorM(i+1))+abs(ErrorM(i)))*dr           
              end do
!!$                 print *, '    Norma ER ', NormaerrorER, '    Norma SR ', NormaerrorSR            
!!$                 print *, '    Norma E ', NormaerrorE, '    Norma M ', NormaerrorM
!!$                 print *, '    Norma Pg ', NormaerrorPg              
!!$                 print *, 'convergenciaER', NormaerrorER - Norma_beforeER
!!$                 print *, 'convergenciaSR', NormaerrorSR - Norma_beforeSR
!!$                 print *, 'convergenciaPg', NormaerrorPg - Norma_beforePg              
!!$                 print *, '**********************************************************'

              if(abs(NormaerrorER - Norma_beforeER).lt.tol_down.or.abs(NormaerrorSR - Norma_beforeSR).lt.tol_down.or. &
                   abs(NormaerrorPg - Norma_beforePg).lt.tol_down ) then
!!$                 print*, "converge al orden", abs(NormaerrorER - Norma_beforeER), p
!!$                 print*, "converge al orden", abs(NormaerrorSR - Norma_beforeSR)
!!$                 print*, "converge al orden", abs(NormaerrorPg - Norma_beforePg)
!                 print*, 'converge2', p                 

                 u(:,:) = 0.0d0
                 u(:,:) = u_1(:,:)

                 go to 2
              end if

           end if

           u(:,:) = 0.0d0
           u(:,:) = u_1(:,:)

           Norma_beforeER = NormaerrorER
           Norma_beforeSR = NormaerrorSR
           Norma_beforePg = NormaerrorPg           
           pg_before(:) = w(3,:)
           
           call BC
           call recv_primitives
         
        end do
        
2       salto=0.0d0


     else if (m.eq.3) then


        u(:,:) = 0.0d0

        do j=1,nvars

           if (j.eq.1) then

              u(j,:) = 1.0d0/2.0d0 * (u_p(j,:)+ u_1(j,:) + dt * rhs_u(j,:))

           else if (j.le.3) then

              u(j,:) = 1.0d0/2.0d0 * (u_p(j,:)+ u_1(j,:) + dt * rhs_u(j,:)) &
                   + dt * (igamma*(-M_0(j-1,:)) + (1.0d0-igamma)/2.0d0 *(- M_1(j-1,:)))
           else
              u(j,:) = 1.0d0/2.0d0 * (u_p(j,:)+ u_1(j,:) + dt * rhs_u(j,:)) &
                   + dt * (igamma*(M_0(j-3,:)) + (1.0d0-igamma)/2.0d0 *(M_1(j-3,:)))
           end if

        end do

     end if


     call BC
     call recv_primitives

  end do
  









end subroutine ImEx2fp
