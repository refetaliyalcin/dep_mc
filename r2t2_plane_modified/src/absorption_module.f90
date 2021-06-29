#include "macros.inc"

module absorption_module
!!Module which handles absorption

    use constants
    use typedefinitions
    use error_handler
    implicit none




contains

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine init_absorption_module(dB)
    !!Default initialization routine for module
    type(dataBlock), intent(in) :: dB
end subroutine

!///////////////////////////////////////////////////////////////////////////////////////////
subroutine absorb(I,Aabs,varpi)
    !!Take absorption from the Stoke's vector I 
    real(kind=rk), intent(inout) :: I(4)        !!Stoke's vector
    real(kind=rk), intent(inout) :: Aabs        !!Cumulative absorbed energy
    real(kind=rk), intent(in) :: varpi          !!The albedo of the single scatterer

    if(debugging) then
        ASSERTC(varpi,>,0)
        ASSERTC(I(1),>,0)
        ASSERTC(Aabs,>=,0.0_rk)
        ASSERTC(varpi,<=,1.0_rk)
    endif

    Aabs=Aabs+(1.0_rk-varpi)*I(1)
    I(:)=varpi*I(:)

    if(debugging) then
        ASSERTC(I(1),>,0)
        ASSERTC(Aabs,>=,0.0_rk)
    endif
end subroutine

end module
