program main

  use constantes
  use mesh
  use vtk
  use ISO_C_BINDING
  use ISO_FORTRAN_ENV

  integer(kind=C_INT) :: nb_elts, nb_MF, nb_noeuds
  real(kind=PR) :: Lx, Ly, Lz

  ! real(kind=PR), dimension(:), allocatable :: solution
  ! integer :: choix

  print*, '------ DEBUT DU PROGRAMME ------'

  print*, '--------------------------------'

  print*, '---- GENERATION DU MAILLAGE ----'

  Lx = 10.0
  Ly = 10.0
  Lz = 10.0

  call generate_mesh(maillage_en_espace,nb_elts,nb_MF,nb_noeuds,Lx,Ly,Lz)


  print*, '---------------------------------'

  print*, '---- TRANSPOSITION DE LA SOL ----'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Rajouter travail sur la sortie graphique :    !!!
  !!! 1) Lire la solution dans un fichier, l'ecrire !!!
  !!!    dans un tableau a 1 dimension et nb_elts   !!!
  !!!    elements                                   !!!
  !!! 2) Ecrire une fonction d'interpolation a      !!!
  !!!    partir de 4 valeurs (moyenne ou autre)     !!!
  !!! 3) Ecrire un fonction pour extrapoler 4       !!!
  !!!    valeurs a partir d'une seule               !!!
  !!! 4) Ecrire une routine pour transposer la sol  !!!
  !!!    sur un maillage uniforme (creer un nouveau !!!
  !!!    vecteur solution mais a plusieurs dim)     !!!
  !!! 5) Enregistrer les deux types de solutions a  !!!
  !!!    l'aide des routines du module vtk          !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print*, '---------------------------------'

  print*, '------ FIN DU PROGRAMME ------'


end program main
