module vtk


  use constantes, only: PR
  use mesh


  implicit none

  private

  public :: ecrit_vtk_AMR, ecrit_vtk_2d

contains



  subroutine ecrit_vtk_2d(ch, q, Lx, Ly)

    real(PR), intent(in) :: Lx, Ly

    real(PR), dimension(:,:,:), intent(in) :: q

    character(len=*), intent(in) :: ch

    character(len=:), allocatable :: str

    integer :: unite, i_f, j_f, i, j

    i_f = size(q, 1)
    j_f = size(q, 2)

    str = trim(adjustl(ch))

    open (newunit=unite, file=str, form='formatted', access='sequential', position='rewind', action ='write')

    write (unite,'(A)') '# vtk DataFile Version 3.0'
    write (unite,'(A)') '# ' // str
    write (unite,'(A)') 'ASCII'
    write (unite,'(A)') 'DATASET STRUCTURED_POINTS'

    write (unite,'(A,2(1x,i5))') 'DIMENSIONS', i_f, j_f

    write (unite,'(A,2(1x,g16.8))') 'ORIGIN', 0.5 * Lx / i_f, 0.5 * Ly / j_f

    write (unite,'(A,2(1x,g16.8))') 'SPACING', Lx / i_f, Ly / j_f

    write (unite,'(A,1x,i8)') 'POINT_DATA', size(q)

    write (unite,'(A)') 'SCALARS ' // str // ' float'
    write (unite,'(A)') 'LOOKUP_TABLE default'

    do j = 1, j_f
      do i = 1, i_f
        write (unite,'(g16.8)') q(i,j,1)
      end do
    end do

    close (unite)

  end subroutine ecrit_vtk_2d

  subroutine ecrit_vtk_AMR(ch, q)

    real(PR), dimension(:), intent(in) :: q

    character(len=*), intent(in) :: ch

    character(len=:), allocatable :: str

    integer :: unite, n_cell, n_nodes, k, n_val

    character(len=80) :: lab

    character(len=2) :: num

    n_cell = size(q, 1)

    n_nodes = size(CoordNoeuds, 1)

    n_val = n_cell

    do k = 1, n_cell

      n_val = n_val + MaillageComplet(k)%nb_noeuds

    end do

    str = trim(adjustl(ch))

    open (newunit=unite, file=str, form='formatted', access='sequential', position='rewind', action ='write')

    write (unite,'(A)') '# vtk DataFile Version 3.0'
    write (unite,'(A)') '# ' // str
    write (unite,'(A)') 'ASCII'
    write (unite,'(A)') 'DATASET UNSTRUCTURED_GRID'

    write (unite,'(A,1x,i8,1x,A)') 'POINTS', n_nodes, 'float'

    do k = 1, n_nodes

      write (unite,'(g16.8,2(1x,g16.8))') CoordNoeuds(k,1), CoordNoeuds(k,2), 0._PR

    end do

    write (unite,'(A,2(1x,i8))') 'CELLS', n_cell, n_val

    do k = 1, n_cell

      write(num,'(i2.2)') MaillageComplet(k)%nb_noeuds

      lab = '(i1,' // num // '(1x,i8))'

      write (unite,lab) MaillageComplet(k)%nb_noeuds, MaillageComplet(k)%noeuds(1:MaillageComplet(k)%nb_noeuds)

    end do

    write (unite,'(A,1x,i8)') 'CELL_TYPES', n_cell

    do k = 1, n_cell

      write (unite,'(i1)') 9

    end do

    write (unite,'(A,1x,i8)') 'CELL_DATA', n_cell
    write (unite,'(A,1x,i1)') 'SCALARS ' // str // ' float', 1
    write (unite,'(A)') 'LOOKUP_TABLE default'

    do k = 1, n_cell
       write (unite,'(g16.8)') q(k)
    end do

    close (unite)

  end subroutine ecrit_vtk_AMR

end module vtk
