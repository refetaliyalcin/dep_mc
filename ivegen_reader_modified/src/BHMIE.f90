!********************************************************************
!
! Bohren Huffman mie coefficients
!
!Bohren, Craig F., and Donald R. Huffman. Absorption and scattering of light by small particles
!********************************************************************
SUBROUTINE BHMIE(NSTOP, X,REFREL,a_n,b_n)
      
    !Arguments:
    real(dp), intent(in) :: X
    COMPLEX(dp), intent(in) :: REFREL 
    ! Local variables:
    INTEGER :: N,NSTOP,NMX,NN
    real(dp) :: CHI,CHI0,CHI1,DX,EN,FN,P,PII,PSI,PSI0,PSI1,YMOD
     

    COMPLEX(dp) :: AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
    COMPLEX(dp) :: D(NSTOP+1)
    COMPLEX(dp), intent(out) :: a_n(NSTOP), b_n(NSTOP)


    !C*** Obtain pi:
    PII=4.*ATAN(1.D0)
    DX=X
    DREFRL=REFREL
    Y=X*DREFRL
    YMOD=ABS(Y)
    !
    !*** Series expansion terminated after NSTOP terms
    !    Logarithmic derivatives calculated from NMX on down
         
    NMX=NSTOP
       
    !*** Logarithmic derivative D(J) calculated by downward recurrence
    !    beginning with initial value (0.,0.) at J=NMX
    !
    D(NMX)=(0.0d0,0.0d0)
    NN=NMX-1
    DO  N=1,NN
       EN=NMX-N+1
       D(NMX-N)=(EN/Y)-(1./(D(NMX-N+1)+EN/Y))
    end do
    !
    !*** Riccati-Bessel functions with real argument X
    !    calculated by upward recurrence
    !
    PSI0=COS(DX)
    PSI1=SIN(DX)
    CHI0=-SIN(DX)
    CHI1=COS(DX)
    XI1=DCMPLX(PSI1,-CHI1)

    P=-1.
    DO N=1,NSTOP
       EN=N
       FN=(2.0d0*EN+1.0d0)/(EN*(EN+1.0d0))
       
       PSI=(2.0d0*EN-1.0d0)*PSI1/DX-PSI0
       CHI=(2.0d0*EN-1.0d0)*CHI1/DX-CHI0
       XI=DCMPLX(PSI,-CHI)
       
       IF(N.GT.1)THEN
          AN1=AN
          BN1=BN
          a_n(N-1) = -an1
          b_n(N-1) = -bn1
          
       ENDIF
       
       AN=(D(N)/DREFRL+EN/DX)*PSI-PSI1
       AN=AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
       BN=(DREFRL*D(N)+EN/DX)*PSI-PSI1
       BN=BN/((DREFRL*D(N)+EN/DX)*XI-XI1)
       
       PSI0=PSI1
       PSI1=PSI
       CHI0=CHI1
       CHI1=CHI
       XI1=DCMPLX(PSI1,-CHI1)
       
       

    END DO
end subroutine bhmie

