
subroutine BC

  use arrays
  use global_numbers

  implicit none

  integer i

!  if (IData.eq.'Rad_Shock_Tube') then

     ! Outflow boundary conditions

     u(:,0) = u(:,1)

     u(:,Nr) = u(:,Nr-1)

!!$  else if (IData.eq.'Rad_Jet') then
!!$
!!$     do i =0,Nr
!!$
!!$        if (r(i).lt.r0) then
!!$
!!$           u(1,i) = u_p(1,i)
!!$           u(2,i) = u_p(2,i)
!!$           u(3,i) = u_p(3,i)
!!$           u(4,i) = u_p(4,i)
!!$           u(5,i) = u_p(5,i)
!!$
!!$        end if
!!$
!!$     end do
!!$
!!$     u(1,Nr) = u(1,Nr-1)
!!$     u(2,Nr) = u(2,Nr-1)
!!$     u(3,Nr) = u(3,Nr-1)
!!$     u(4,Nr) = u(4,Nr-1)
!!$     u(5,Nr) = u(5,Nr-1)
!!$
!!$  end if
  


end subroutine BC
