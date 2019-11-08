module p4est_binding

  use ISO_C_BINDING

  interface


    subroutine p4_build_p4est(connectivity,p4est) bind(C)

      !implicit none
      import C_PTR


      ! Variables en sortie
      type(C_PTR)     :: p4est
      type(C_PTR)     :: connectivity

    end subroutine p4_build_p4est

    subroutine p4_raffiner(p4est, refine_recursive, nb_elts, ghost, mesh, noeuds, nb_noeuds) bind(C)

      !implicit none
      import C_INT, C_PTR

      ! Variables en entree
      integer(kind=C_INT),VALUE     :: refine_recursive
      type(C_PTR),VALUE             :: p4est

      !Variables en sortie
      integer(kind=C_INT),intent(out)   :: nb_elts
      type(C_PTR)                       :: ghost
      type(C_PTR)                       :: mesh
      type(C_PTR)                       :: noeuds
      integer(kind=C_INT),intent(out)   ::nb_noeuds

    end subroutine p4_raffiner


    subroutine p4_browse_mesh(p4est,ghost,mesh,compteur) bind(C)

      !implicit none
      import C_PTR, C_INT

      ! Variables en entree
      type(C_PTR),VALUE     :: p4est
      type(C_PTR),VALUE     :: ghost
      type(C_PTR),VALUE     :: mesh

      !Variables en sortie
      integer(kind=C_INT),intent(out)  :: compteur


    end subroutine p4_browse_mesh



    subroutine p4_iterate_mesh2(p4est,ghost,mesh,noeuds,num_quad,maillage_en_espace,area,length,x,y,z, &
                                  & nb_voisins,voisins_par_face,h,nb_noeuds_cel,noeuds_cel) bind(C)

      !implicit none
      import C_PTR, C_INT, C_DOUBLE

      ! Variables en entree
      type(C_PTR),VALUE             :: p4est
      type(C_PTR),VALUE             :: mesh
      type(C_PTR),VALUE             :: ghost
      type(C_PTR),VALUE             :: noeuds
      integer(kind=C_INT), VALUE    :: num_quad
      integer(kind=C_INT), VALUE    :: maillage_en_espace


      !Variables en sortie
      real(kind=C_DOUBLE),intent(out)                    :: x
      real(kind=C_DOUBLE),intent(out)                    :: y
      real(kind=C_DOUBLE),intent(out)                    :: z
      real(kind=C_DOUBLE),intent(out)                    :: area
      real(kind=C_DOUBLE),intent(out)                    :: length
      integer(kind=C_INT),dimension(0:3,0:1),intent(out) :: voisins_par_face
      integer(kind=C_INT),intent(out)                    :: nb_voisins
      real(kind=C_DOUBLE),intent(out)                    :: h
      integer(kind=C_INT),intent(out)                    :: nb_noeuds_cel
      integer(kind=C_INT),dimension(0:7),intent(out)     :: noeuds_cel

    end subroutine p4_iterate_mesh2



    subroutine p4_get_node(p4est,noeuds,num_node,x,y) bind(C)

      !implicit none
      import C_PTR, C_INT, C_DOUBLE

      ! Variables en entree
      type(C_PTR),VALUE             :: p4est
      type(C_PTR),VALUE             :: noeuds
      integer(kind=C_INT), VALUE    :: num_node


      !Variables en sortie
      real(kind=C_DOUBLE),intent(out)                    :: x
      real(kind=C_DOUBLE),intent(out)                    :: y


    end subroutine p4_get_node



    subroutine p4_save_vtk(num_test, p4est, connectivity) bind(C)

      !implicit none
      import C_PTR, C_CHAR, C_INT

      ! Variables en entree
      integer(kind=C_INT),VALUE     :: num_test
      type(C_PTR),VALUE             :: p4est
      type(C_PTR),VALUE             :: connectivity
      ! integer(kind=C_INT),VALUE     :: choix_enregistrement

    end subroutine p4_save_vtk


  end interface


end module p4est_binding
