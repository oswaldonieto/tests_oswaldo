

subroutine recv_primitives

  use arrays
  use global_numbers
  use lu

  implicit none

  integer i,l
  integer d, code  

  real(kind=8) fun, funp
  real(kind=8) tol 

  SS = sqrt(u(2,:)**2)

  do i = 0,Nr

     P_G(i) = 0.5d0 * (piL + piR)
     P_S(i) = P_G(i)

     do l = 1,10000000

        P_S_p(i) = P_S(i)

        fun = - P_S_p(i) + (gamma - 1.0d0)*u(3,i) &
             - (SS(i)**2*(gamma - 1.0d0))/(P_S_p(i) + u(3,i)) &
             - u(1,i)*(gamma - 1.0d0)*(sqrt(1.0d0 - SS(i)**2/(( &
             + P_S_p(i) + u(3,i))**2)))

        funp = - 1.0d0 + (SS(i)**2*(gamma - 1.0d0))/((P_S_p(i) + u(3,i))**2) &
             - (u(1,i)*SS(i)**2*(gamma - 1.0d0))/((P_S_p(i) &
             + u(3,i))**3*sqrt(1.0d0 - SS(i)**2/((P_S_p(i) + u(3,i))**2)))

!!$        fun = - P_S_p(i) + (gamma - 1.0d0)*u(3,i) &
!!$             - (SS(i)**2*(gamma - 1.0d0))/(P_S_p(i) + u(3,i) + u(1,i)) &
!!$             - u(1,i)*(gamma - 1.0d0)*(-1.0d0 + sqrt(1.0d0 - SS(i)**2/(( u(3,i)&
!!$             + P_S_p(i) + u(1,i))**2)))
!!$
!!$        funp = - 1.0d0 + (SS(i)**2*(gamma - 1.0d0))/((u(3,i) + P_S_p(i) + u(1,i))**2) &
!!$             - (u(1,i)*SS(i)**2*(gamma - 1.0d0))/((u(3,i) + P_S_p(i) &
!!$             + u(1,i))**3*sqrt(1.0d0 - SS(i)**2/((u(3,i) + P_S_p(i) + u(1,i))**2)))        
        

        P_S(i) = P_S_p(i) - fun/funp

        tol = 2.0d0*abs(P_S(i) - P_S_p(i))/(P_S(i) + P_S_p(i))

        if (tol.lt.1.e-6) then
           w(3,i) = P_S(i)
           exit
        else
        end if
        if (l.gt.1000000) then
           print *, 'The Newton_Raphson routine did not converge'
           print *, 'l=',l
           print *, 'i=',i
           print *, 'r=',r(i)
           stop
        end if

     end do
  end do

  w(3,:) = max(w(3,:),press_atm)
    
  w(2,:) = u(2,:)/(w(3,:) + u(3,:))
 ! w(2,:) = u(2,:)/(u(3,:) + w(3,:) + u(1,:))  
  
  LF  = 1.0d0/sqrt(1.0d0 - w(2,:)**2)
  
  w(1,:) = max(u(1,:)/LF,rho_atm) 

  CS = sqrt(gamma*w(3,:)/(w(1,:) + w(3,:)*gamma/(gamma - 1.0d0)))

!!$  w(4,:) = 3.0d0*( 2.0d0*u(4,:)*w(2,:) - u(5,:)*(w(2,:)**2 + 1.0d0) )
!!$  w(4,:) = w(4,:)/( 4.0d0*w(2,:)**2*LF**2 - 4.0d0*LF**2 + w(2,:)**2 + 1.0d0 )
!!$ 
!!$  w(5,:) = (1.0d0 - 4.0d0*LF**2)*u(4,:) + 4.0d0*LF**2*w(2,:)*u(5,:)
!!$  w(5,:) = w(5,:)/( LF*( 4.0d0*w(2,:)**2*LF**2 - 4.0d0*LF**2 + w(2,:)**2 + 1.0d0 ) )


!!$  w(4,:) = 2.0d0*u(4,:)*w(2,:) + u(5,:)*(1.0d0/LF(:)**2 -2.0d0)
!!$  w(4,:) = w(4,:) * ( -3.0d0*LF(:)**2/(1+2*LF(:)**2))
!!$
!!$  w(5,:) = u(4,:)/LF(:) - 4.0d0/3.0d0 * (w(4,:)*LF(:)*w(2,:))
!!$  w(5,:) = w(5,:) - (LF(:)/(1.0d0+2.0d0*LF(:)**2))*(-4.0d0*u(5,:)*(LF(:)**2-1)+ &
!!$       (4.0d0*LF(:)**2-1.0d0)*u(4,:)*w(2,:))

  w(4,:) = 0.0d0
  w(5,:) = 0.0d0
  w_rad(:,:)= 0.0d0
  u_rad(:,:)= 0.0d0  

  do i = 0, Nr
     
     w_rad(1,1) = (4.0d0/3.0d0)*LF(i)**2*w(2,i)
     w_rad(1,2) = LF(i) * (w(2,i)+1.0d0)
     w_rad(2,1) = (4.0d0*LF(i)**2-1.0d0)/3.0d0
     w_rad(2,2) = 2.0d0*LF(i)*w(2,i)

     u_rad(1,1) = u(4,i)
     u_rad(2,1) = u(5,i)

     call LUDCMP(w_rad,2,INDX,D,CODE)
     call LUBKSB(w_rad,2,INDX,u_rad)

     w(4,i) = u_rad(1,1)
     w(5,i) = u_rad(2,1)

  end do
  
         
end subroutine recv_primitives
