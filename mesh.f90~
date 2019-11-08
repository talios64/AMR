module mesh

  use constantes
  use p4est_binding
  use iso_fortran_env

  implicit none

  private


  type, public ::  face_cellule

    integer(kind=C_INT) :: n_face !numero de la face de la cellule a laquelle se trouve la liste voisins

    integer(kind=C_INT) :: nb_voisins !nombre de voisins de la cellule a cette face

    integer(kind=C_INT), allocatable, dimension(:) :: voisin !liste des voisins a cette face

  end type face_cellule


  type, public :: voisins

    integer(kind=C_INT) :: nb_faces !nombre de faces d'une cellule (2 en 1D, 4 en 2D et 6 en 3D)

    real(kind=PR) :: hx, hy, hz !pas d'espace dans chaque direction

    real(kind=PR) :: x, y, z !coordonnees du centre de la cellule

    type(face_cellule), allocatable, dimension(:) :: face !liste des voisins par face

    integer(kind=C_INT) :: nb_noeuds !nombre de noeuds constituant la cellule

    integer(kind=C_INT), allocatable, dimension(:) :: noeuds !liste des numeros des noeuds de la cellule

  end type voisins

  type(voisins), dimension(:), allocatable, save, public, protected :: MaillageComplet

  real(kind=C_DOUBLE), dimension(:,:), allocatable, save, public, protected :: CoordNoeuds

  public :: generate_mesh


contains

  subroutine generate_nodes_2d(p4est,noeuds,nb_noeuds,Lx,Ly,Lz)

    type(C_PTR), intent(in)               :: p4est, noeuds
    integer(kind=C_INT), intent(in)       :: nb_noeuds
    real(kind=PR), intent(in)             :: Lx,Ly,Lz

    integer :: i
    real(kind=C_DOUBLE) :: x,y

    allocate(CoordNoeuds(nb_noeuds,2))

    do i = 0, nb_noeuds-1

      call p4_get_node(p4est,noeuds,i,y,x)

      CoordNoeuds(i+1,1) = x * Lx
      CoordNoeuds(i+1,2) = y * Ly

    end do

  end subroutine generate_nodes_2d



  subroutine generate_mesh(maillage_en_espace,nb_elts,nb_MF,nb_noeuds,Lx,Ly,Lz)

    integer(kind=C_INT), intent(inout)    :: nb_elts, nb_MF, nb_noeuds
    integer(kind=C_INT), intent(in)       :: maillage_en_espace
    real(kind=PR), intent(in)             :: Lx,Ly,Lz

    type(C_PTR)              :: p4est
    type(C_PTR)              :: ghost
    type(C_PTR)              :: mesh
    type(C_PTR)              :: noeuds
    type(C_PTR)              :: connectivity
    integer(kind=C_INT)      :: recursive

    integer(kind=C_INT)       :: num_quad, nb_voisins, num_test, nb_noeuds_cel
    real(kind=C_DOUBLE)       :: x, y, z, area, length, h
    integer                   :: i, j, k, l
    integer(kind=C_INT), dimension(6) :: numero_faces
    integer(kind=C_INT), dimension(0:7) :: noeuds_cel
    integer(kind=C_INT), dimension(:,:), allocatable :: voisins_par_face, voisins_par_face_temp

    if (maillage_en_espace==2) then

      print*, '-----------Initialisation p4est------------'
      call p4_build_p4est(connectivity,p4est)
      print*, 'Fait.'

      print*, '-----------Raffinement-----------'
      recursive = 1
      call p4_raffiner(p4est, recursive, nb_elts, ghost, mesh, noeuds, nb_noeuds)
      print*, 'Fait.'
      print*, 'Nombre de cellules du domaine : ', nb_elts
      print*, 'Nombre de noeuds du domaine : ', nb_noeuds

      print*, '-----------Parcours du maillage (browse)-----------'
      call p4_browse_mesh(p4est, ghost, mesh, nb_MF)
      print*, 'Fait.'

      print*, '-----------Enregistrement du maillage-----------'
      print*, 'numero fichier '
      read(*,*) num_test
      call p4_save_vtk(num_test, p4est, connectivity)
      print*, 'Fait.'

    else

      ! print*, '-----------Initialisation p4est------------'
      ! call p8_build_p4est(connectivity,p4est)
      ! print*, 'Fait.'
      !
      ! print*, '-----------Raffinement-----------'
      ! recursive = 1
      ! call p8_raffiner(p4est, recursive, nb_elts, noeuds, nb_noeuds)
      ! print*, 'Fait.'
      ! print*, 'Nombre de cellules du domaine : ', nb_elts
      !
      ! print*, '-----------Parcours du maillage (browse)-----------'
      ! call p8_browse_mesh(p4est, nb_MF, ghost, mesh)
      ! print*, 'Fait.'
      !
      ! print*, '-----------Enregistrement du maillage-----------'
      ! print*, 'numero fichier '
      ! read(*,*) num_test
      ! call p8_save_vtk(num_test, p4est, connectivity)
      ! print*, 'Fait.'

    end if

    print*, '-----------Transfert du maillage (browse)-----------'
    allocate(MaillageComplet(nb_elts))
    allocate(voisins_par_face(0:2**(maillage_en_espace-1)-1,0:maillage_en_espace*2-1)) !lignes et colonnes echangees à cause du C
    allocate(voisins_par_face_temp(0:2**(maillage_en_espace-1)-1,0:maillage_en_espace*2-1))

    numero_faces(1)=1  !face +x
    numero_faces(2)=-1 !face -x
    numero_faces(3)=2  !face +y
    numero_faces(4)=-2 !face -y
    numero_faces(5)=3  !face +z
    numero_faces(6)=-3 !face -z

    l = 1

    do i = 1, nb_elts

      MaillageComplet(i)%nb_faces = maillage_en_espace*2

      allocate(MaillageComplet(i)%face(MaillageComplet(i)%nb_faces))

      num_quad = i-1

      if (maillage_en_espace==2) then

      call p4_iterate_mesh2(p4est,ghost,mesh,noeuds,num_quad,maillage_en_espace,area,length,y,x,z, &
                                    & nb_voisins,voisins_par_face_temp,h,nb_noeuds_cel,noeuds_cel)

      else

        ! call p8_iterate_mesh2(p4est,ghost,mesh,noeuds,num_quad,maillage_en_espace,area,length,y,x,z, &
        !                               & nb_voisins,voisins_par_face_temp,h)

      end if

      voisins_par_face(:,0) = voisins_par_face_temp(:,2)
      voisins_par_face(:,1) = voisins_par_face_temp(:,3)
      voisins_par_face(:,2) = voisins_par_face_temp(:,0)
      voisins_par_face(:,3) = voisins_par_face_temp(:,1)
      if (maillage_en_espace==3) then
        voisins_par_face(:,4) = voisins_par_face_temp(:,4)
        voisins_par_face(:,5) = voisins_par_face_temp(:,5)
      end if

      do j = 1, MaillageComplet(i)%nb_faces

        MaillageComplet(i)%face(j)%n_face = numero_faces(j)

        MaillageComplet(i)%face(j)%nb_voisins = 0

        do k = 1, 2**(maillage_en_espace-1)

          if (voisins_par_face(k-1,j-1)/=-1) then

            MaillageComplet(i)%face(j)%nb_voisins = MaillageComplet(i)%face(j)%nb_voisins + 1

          end if

        end do

        allocate(MaillageComplet(i)%face(j)%voisin(MaillageComplet(i)%face(j)%nb_voisins))

        do k = 1, MaillageComplet(i)%face(j)%nb_voisins

          if (voisins_par_face(k-1,j-1)/=num_quad) then

            MaillageComplet(i)%face(j)%voisin(k) = voisins_par_face(k-1,j-1) + 1

          else

            MaillageComplet(i)%face(j)%voisin(k) = nb_elts + l

            l = l + 1

          end if

          if (l > nb_MF+1) then

            print*, 'erreur mailles fictives, l = ', l

          end if

        end do

      end do

      MaillageComplet(i)%x = x * Lx

      MaillageComplet(i)%hx = h * Lx

      MaillageComplet(i)%y = y * Ly

      MaillageComplet(i)%hy = h * Ly

      if (maillage_en_espace==3) then

        MaillageComplet(i)%z = z * Lz

        MaillageComplet(i)%hz = h * Lz

      else

        MaillageComplet(i)%z = Lz / 2.0

        MaillageComplet(i)%hz = Lz

      end if

      if (maillage_en_espace==2) then

        MaillageComplet(i)%nb_noeuds = nb_noeuds_cel

        allocate(MaillageComplet(i)%noeuds(nb_noeuds_cel))

        do j = 1, nb_noeuds_cel

          MaillageComplet(i)%noeuds(j) = noeuds_cel(j-1) - 1

        end do

      end if

      ! print*, 'Info maille ', i, ' : '
      ! print*, ' coord centre = ', MaillageComplet(i)%x, MaillageComplet(i)%y, MaillageComplet(i)%z
      ! ! print*, ' aire = ', MaillageComplet(i)%area
      ! print*, ' dx = ', MaillageComplet(i)%hx
      ! print*, ' liste voisins face +x : ', MaillageComplet(i)%face(1)%voisin(:)
      ! print*, ' liste voisins face -x : ', MaillageComplet(i)%face(2)%voisin(:)
      ! print*, ' dy = ', MaillageComplet(i)%hy
      ! print*, ' liste voisins face +y : ', MaillageComplet(i)%face(3)%voisin(:)
      ! print*, ' liste voisins face -y : ', MaillageComplet(i)%face(4)%voisin(:)
      ! if (maillage_en_espace==3) then
      !   print*, ' dz = ', MaillageComplet(i)%hz
      !   print*, ' liste voisins face +z : ', MaillageComplet(i)%face(5)%voisin(:)
      !   print*, ' liste voisins face -z : ', MaillageComplet(i)%face(6)%voisin(:)
      ! end if

    end do
    print*, 'Fait.'

    print*, '-----------Generate nodes-----------'
    if (maillage_en_espace==2) then
      call generate_nodes_2d(p4est,noeuds,nb_noeuds,Lx,Ly,Lz)
    end if
    print*, 'Fait.'

  end subroutine generate_mesh


end module mesh
