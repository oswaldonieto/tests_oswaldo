
subroutine flux_hlle

  use arrays
  use global_numbers

  implicit none

  integer i
  
  real(kind=8) l_min,l_max
  real(kind=8) vL,vR
  real(kind=8) LFL,LFR
  real(kind=8) CSL,CSR
  real(kind=8) Cr
 
  real(kind=8), dimension (1:nvars) :: FrL, FrR
  real(kind=8), dimension (1:nvars) :: uL, uR
  real(kind=8), dimension (1:nvars) :: wL, wR
  real(kind=8), dimension (1:nvars) :: lL, lR

  Cr = sqrt(1.0d0/3.0d0)

  do i=0,Nr-1

     if (Limiter.eq.'mm' .or. Limiter.eq.'mc') then

        call lpm(i,Nr,dr,w(1,i-1),w(1,i),w(1,i+1),w(1,i+2),wL(1),wR(1))
        call lpm(i,Nr,dr,w(2,i-1),w(2,i),w(2,i+1),w(2,i+2),wL(2),wR(2))
        call lpm(i,Nr,dr,w(3,i-1),w(3,i),w(3,i+1),w(3,i+2),wL(3),wR(3))
        call lpm(i,Nr,dr,w(4,i-1),w(4,i),w(4,i+1),w(4,i+2),wL(4),wR(4))
        call lpm(i,Nr,dr,w(5,i-1),w(5,i),w(5,i+1),w(5,i+2),wL(5),wR(5))
        call lpm(i,Nr,dr,CS(i-1),CS(i),CS(i+1),CS(i+2),CSL,CSR)

     end if

     ! Lorentz Factor a izquierda y a derecha

     vL = sqrt(wL(2)**2)
     vR = sqrt(wR(2)**2)

     LFL = 1.0d0/sqrt(1.0d0 - vL**2)
     LFR = 1.0d0/sqrt(1.0d0 - vR**2)

!!$     CSL = sqrt(gamma*wL(3)/(wL(1) + wL(3)*gamma/(gamma - 1.0d0)))
!!$     CSR = sqrt(gamma*wR(3)/(wR(1) + wR(3)*gamma/(gamma - 1.0d0)))     

     ! Variables conservativas a izquierda y a derecha

     uL(1) = wL(1)*LFL
     uR(1) = wR(1)*LFR

     uL(2) = (wL(1) + gamma*wL(3)/(gamma - 1.0d0))*LFL**2*wL(2)
     uR(2) = (wR(1) + gamma*wR(3)/(gamma - 1.0d0))*LFR**2*wR(2)

     uL(3) = (wL(1) + gamma*wL(3)/(gamma - 1.0d0))*LFL**2 - wL(3) 
     uR(3) = (wR(1) + gamma*wR(3)/(gamma - 1.0d0))*LFR**2 - wR(3)

!!$     uL(3) = (wL(1) + gamma * wL(3)/(gamma - 1.0d0)) * LFL**2 - wL(3) - wL(1)*LFL
!!$     uR(3) = (wR(1) + gamma * wR(3)/(gamma - 1.0d0)) * LFR**2 - wR(3) - wR(1)*LFR     

     uL(4) = 4.0d0*wL(4)*LFL**2*wL(2)/3.0d0 + LFL*wL(5)*(wL(2)**2 + 1.0d0)
     uR(4) = 4.0d0*wR(4)*LFR**2*wR(2)/3.0d0 + LFR*wR(5)*(wR(2)**2 + 1.0d0)     

     uL(5) = 4.0d0*wL(4)*LFL**2/3.0d0 + 2.0d0*LFL*wL(5)*wL(2) - wL(4)/3.0d0
     uR(5) = 4.0d0*wR(4)*LFR**2/3.0d0 + 2.0d0*LFR*wR(5)*wR(2) - wR(4)/3.0d0
          
     ! ***************************************************************************
     ! ***************************************************************************

     ! HLL STAFF

     ! Eigenvalues along x-direction

     lL(1) = wL(2)
     lR(1) = wR(2)

     lL(2) =  (1.0d0/(1.0d0 - vL**2*CSL**2))*(wL(2)*(1.0d0 - CSL**2) &
          + sqrt(CSL**2*(1.0d0 - vL**2)*((1.0d0 - vL**2*CSL**2) &
          - wL(2)**2*(1.0d0 - CSL**2))))

     lR(2) = (1.0d0/(1.0d0 - vR**2*CSR**2))*(wR(2)*(1.0d0 - CSR**2) &
          + sqrt(CSR**2*(1.0d0 - vR**2)*((1.0d0 - vR**2*CSR**2) &
          - wR(2)**2*(1.0d0 - CSR**2))))

     lL(3) = (1.0d0/(1.0d0 - vL**2*CSL**2))*(wL(2)*(1.0d0 - CSL**2) &
          - sqrt(CSL**2*(1.0d0 - vL**2)*((1.0d0 - vL**2*CSL**2) &
          - wL(2)**2*(1.0d0 - CSL**2))))

     lR(3) = (1.0d0/(1.0d0 - vR**2*CSR**2))*(wR(2)*(1.0d0 - CSR**2) &
          - sqrt(CSR**2*(1.0d0 - vR**2)*((1.0d0 - vR**2*CSR**2) &
          - wR(2)**2*(1.0d0 - CSR**2))))


     lL(4) = Cr
     lR(4) = Cr

     lL(5) = -Cr
     lR(5) = -Cr
   
     l_max = max(0.0d0,lL(1),lR(1),lL(2),lR(2),lL(3),lR(3),lL(4),lR(4),lL(5),lR(5))
     l_min = min(0.0d0,lL(1),lR(1),lL(2),lR(2),lL(3),lR(3),lL(4),lR(4),lL(5),lR(5))

     ! ::::::::::::::::::::::::::::::::

     FrL(1) = uL(1)*wL(2)
     FrR(1) = uR(1)*wR(2)

     FrL(2) = uL(2)*wL(2) + wL(3)
     FrR(2) = uR(2)*wR(2) + wR(3)

     FrL(3) = uL(2) 
     FrR(3) = uR(2)

!!$     FrL(3) = uL(2) - uL(1)*wL(2)
!!$     FrR(3) = uR(2) - uR(1)*wR(2)     

     FrL(4) = 2.0d0*LFL*wL(2)*wL(5) + wL(4)/3.0d0 + 4.0d0*wL(4)*LFL**2*wL(2)**2/3.0d0 
     FrR(4) = 2.0d0*LFR*wR(2)*wR(5) + wR(4)/3.0d0 + 4.0d0*wR(4)*LFR**2*wR(2)**2/3.0d0

     FrL(5) = uL(4)
     FrR(5) = uR(4)

     ! Fluxes at the each intercell
     ! ============================


     Fr(1,i) = (l_max*FrL(1) - l_min*FrR(1) + l_max*l_min*(uR(1)-uL(1)))/(l_max-l_min)
     Fr(2,i) = (l_max*FrL(2) - l_min*FrR(2) + l_max*l_min*(uR(2)-uL(2)))/(l_max-l_min)
     Fr(3,i) = (l_max*FrL(3) - l_min*FrR(3) + l_max*l_min*(uR(3)-uL(3)))/(l_max-l_min)
     Fr(4,i) = (l_max*FrL(4) - l_min*FrR(4) + l_max*l_min*(uR(4)-uL(4)))/(l_max-l_min)     
     Fr(5,i) = (l_max*FrL(5) - l_min*FrR(5) + l_max*l_min*(uR(5)-uL(5)))/(l_max-l_min)
     
  end do


end subroutine flux_hlle
