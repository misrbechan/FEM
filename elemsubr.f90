! File: elemsubr.f90
      subroutine elemsubr ( npelm, x, y, nunk_pel, elem_mat,  &
                            elem_vec, elem_mass, prevsolution, itype )
!
!                       INPUT / OUTPUT PARAMETERS
!
      implicit none
      integer, intent(in) :: npelm, nunk_pel, itype
      double precision, intent(in) :: x(1:npelm), y(1:npelm),  &
                                      prevsolution(1:nunk_pel)
      double precision, intent(out) :: elem_mat(1:nunk_pel,1:nunk_pel),  &
                                       elem_vec(1:nunk_pel),  &
                                       elem_mass(1:nunk_pel)

!
!!
!
!
!for example:
integer :: i, j
double precision :: beta(1:3), gamma(1:3), delta, h
real(8), parameter :: a = 0.01
real(8) :: k


delta = (x(2)-x(1))*(y(3)-y(1))-(y(2)-y(1))*(x(3)-x(1))
h = sqrt ( (x(2)-x(1))**2 + (y(2)-y(1))**2 )
beta(1) = (y(2)-y(3))/delta
beta(2) = (y(3)-y(1))/delta
beta(3) = (y(1)-y(2))/delta
gamma(1) = (x(3)-x(2))/delta
gamma(2) = (x(1)-x(3))/delta
gamma(3) = (x(2)-x(1))/delta

if ( itype==1 ) then
! statements to fill the arrays elem_mat and elem_vec
if (((x(1) < 0.2) .and. (y(1) < 0.7) .and. (y(1) > 0.3)) .or. &
((x(2) < 0.2) .and. (y(2) < 0.7) .and. (y(2) > 0.3)) .or. &
((x(3) < 0.2) .and. (y(3) < 0.7) .and. (y(3) > 0.3))) then

k = 1d0 

else

k = 0.1d0

endif
do j = 1, nunk_pel
do i = 1, nunk_pel
elem_mat(i,j) = 0.5d0 * k * abs(delta) * ( beta(i)*beta(j) + gamma(i)*gamma(j) )
end do
elem_vec(j) = 0d0
end do

else if ( itype==2 ) then
!the same type of statements for itype = 2, etcetera
do j = 1, nunk_pel
do i = 1, nunk_pel
if ( i==j ) then
elem_mat(i,j) = h*a/3d0
else
elem_mat(i,j) = h*a/6d0
end if
end do
elem_vec(j) = 0d0
end do

end if

      end subroutine elemsubr
