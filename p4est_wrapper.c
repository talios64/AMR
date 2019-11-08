#include <p4est_connectivity.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_vtk.h>
#include <p4est_mesh.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>
#include <p4est_mesh.h>
#include <p4est_nodes.h>
#include <p4est_geometry.h>
#include <p4est.h>

#include "entrees.h"
#include "dicom.h"
#include <math.h>



/** The resolution of the image data in powers of two. */
#define P4EST_STEP1_PATTERN_LEVEL 9
/** The dimension of the image data. */
#define P4EST_STEP1_PATTERN_LENGTH (1 << P4EST_STEP1_PATTERN_LEVEL)
static const int    plv = P4EST_STEP1_PATTERN_LEVEL;    /**< Shortcut */
static const int    ple = P4EST_STEP1_PATTERN_LENGTH;   /**< Shortcut */
static const int    niveau_min = 3;
#ifdef P4_TO_P8
static const p4est_qcoord_t eighth = P4EST_QUADRANT_LEN (3);
#endif


//Pour recuperer le centre d'un quadrant
void p4_get_midpoint( p4est_t *p4est,
                      p4est_quadrant_t *q,
                      double xyz[3])
{
  p4est_qcoord_t      half_length;

  half_length = P4EST_QUADRANT_LEN (q->level) / 2;

  p4est_qcoord_to_vertex (p4est->connectivity, q->p.which_tree,
                          q->x + half_length, q->y + half_length,
                          xyz);

}



// Fonctions callback
//
// fonction d'initialisation
void p4_init_callback(p4est_t *p4est, int quad_id, p4est_quadrant_t *q)
{

}
//
// Critere raffinement : valeurs du fichier entrees.h
int p4_refine_callback_new( p4est_t *p4est,
                        p4est_topidx_t which_tree,
                        p4est_quadrant_t *q )
{
  double xyz[3];
  p4_get_midpoint(p4est,q,xyz);

  double demi_largeur_faisceau = largeur_faisceau / (2.*longueurs[face_entree]);
  double marge_normalisee = marge / longueurs[face_entree];
  double borne = demi_largeur_faisceau + marge_normalisee;

  for (int zone=0; zone<nb_zones; zone++){
    if ((xyz[1]>=xmin[zone]/longueurs[1])&&(xyz[1]<=xmax[zone]/longueurs[1])&&(xyz[0]>=ymin[zone]/longueurs[0])&&(xyz[0]<=ymax[zone]/longueurs[0])){
      if (q->level < niveau_min_p4est[zone]) {return 1;}
      if ((xyz[face_entree] >= 0.5 - borne )&&(xyz[face_entree] <= 0.5 + borne )&&(q->level < niveau_max_p4est[zone])) {return 1;}
    }
  }
  return 0;
}


// Critere raffinement : sur une image
static int refine_fn (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
  int                 tilelen;
  int                 offsi, offsj;
  int                 i, j;
  const char         *d;
  unsigned char       p[3];

  /* The connectivity chosen in main () only consists of one tree. */
  P4EST_ASSERT (which_tree == 0);

  /* We do not want to refine deeper than a given maximum level. */
  if (quadrant->level > plv) {
    return 0;
  }
  if (quadrant->level < niveau_min) {
    return 1;
  }
#ifdef P4_TO_P8
  /* In 3D we extrude the 2D image in the z direction between [3/8, 5/8]. */
  if (quadrant->level >= 3 &&
      (quadrant->z < 3 * eighth || quadrant->z >= 5 * eighth)) {
    return 0;
  }
#endif

  /* We read the image data and refine wherever the color value is dark.
   * We can then visualize the output and highlight level > PATTERN_LEVEL. */
  tilelen = 1 << (plv - quadrant->level);       /* Pixel size of quadrant */
  // offsi = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  // offsj = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  offsi = quadrant->y / P4EST_QUADRANT_LEN (plv);       /* Pixel x offset */
  offsj = quadrant->x / P4EST_QUADRANT_LEN (plv);       /* Pixel y offset */
  P4EST_ASSERT (offsi >= 0 && offsj >= 0);
  for (j = 0; j < tilelen; ++j) {
    P4EST_ASSERT (offsj + j < ple);
    for (i = 0; i < tilelen; ++i) {
      P4EST_ASSERT (offsi + i < ple);
      d =
        header_data + 4 * (ple * (ple - 1 - (offsj + j)) + (offsi + i));
      HEADER_PIXEL (d, p);
      P4EST_ASSERT (p[0] == p[1] && p[1] == p[2]);      /* Grayscale image */
      if (p[0] < 128) {
        return 1;
      }
    }
  }
  return 0;
}
// Fin des fonctions callback



// Create a forest
void p4_build_p4est ( p4est_connectivity_t **conn_out,
                      p4est_t **p4est_out)
  {
    int          mpiret;
    sc_MPI_Comm  mpicomm;

    p4est_t *p4est=NULL;
    p4est_connectivity_t *conn=NULL;

    /* Initialize MPI; see sc_mpi.h.
    * If configure --enable-mpi is given these are true MPI calls.
    * Else these are dummy functions that simulate a single-processor run. */
    mpiret = 0;
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    /* These functions are optional.  If called they store the MPI rank as a
    * static variable so subsequent global p4est log messages are only issued
    * from processor zero.  Here we turn off most of the logging; see sc.h. */
    sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init (NULL, SC_LP_PRODUCTION);

    // conn = p8est_connectivity_new_unitcube();
    conn = p4est_connectivity_new_unitsquare();

    p4est = p4est_new (mpicomm, conn, 0, NULL, NULL);

    *conn_out = conn;

    *p4est_out = p4est;

  }


//Raffinement
void p4_raffiner(p4est_t *p4est,
                 int  	      refine_recursive, //booleen pour raffiner recursivement ou non
                 int          *nb_elts,
                 p4est_ghost_t **ghost_out,
                 p4est_mesh_t **mesh_out,
                 p4est_nodes_t **nodes_out,
                 int *np_out)
{
  p4est_mesh_t *mesh=NULL;
  p4est_ghost_t *ghost=NULL;
  p4est_nodes_t *nodes;

  /////////////// JOUER SUR LA FONCTION CRITERE DE RAFFINEMENT ICI ///////////////
  int (*refine)(p4est_t *, int, p4est_quadrant_t *) = &p4_refine_callback_new;
  // int (*refine)(p4est_t *, int, p4est_quadrant_t *) = &refine_fn;
  ////////////////////////////////////////////////////////////////////////////////


  void (*init)(p4est_t *, int, p4est_quadrant_t *) = &p4_init_callback;

  p4est_refine(p4est, refine_recursive, refine , init);

  p4est_balance(p4est, P4EST_CONNECT_FACE, init); //attention a la taille entre des mailles voisines

  *nb_elts = p4est->local_num_quadrants ;

  ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
  mesh = p4est_mesh_new_ext(p4est, ghost,1,1, P4EST_CONNECT_FACE);
  nodes = p4est_nodes_new(p4est,NULL);

  *nodes_out = nodes;
  *mesh_out = mesh;
  *ghost_out = ghost;
  *np_out = nodes->indep_nodes.elem_count;

  // printf("------------------------\n" );
  // p4est_locidx_t numquads = p4est->local_num_quadrants;
  // int numquads_refine_unif = pow(4,niveau_max);
  // printf("Nombre de noeuds du domaine = %d \n", &np_out);
  // printf("Nombre de cellules (avec p4est) = %d \n", numquads);
  // printf("Nombre de cellules (sans p4est) = %d \n", numquads_refine_unif);
  // printf("------------------------\n" );

}



// Pour se balader sur la totalite du maillage
void p4_browse_mesh ( p4est_t * p4est,
                      p4est_ghost_t *ghost,
                      p4est_mesh_t *mesh,
                      int *compteur)
{
  // p4est_mesh_t *mesh;
  p4est_quadrant_t *q;
  // p4est_locidx_t q_idx;
  p4est_topidx_t which_tree;
  // p4est_ghost_t *ghost;
  p4est_locidx_t  quadrant_id, which_quad;
  p4est_mesh_face_neighbor_t mfn;
  int                 nface;
  int                 nrank;
  double xyz[3];
  int tmp=0;


  // ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
  // mesh = p4est_mesh_new_ext(p4est, ghost,1,1, P4EST_CONNECT_FACE);

  // printf("(print C) ----------------------------------------------------\n" );

  for (int qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
    which_tree = -1;
    q = p4est_mesh_quadrant_cumulative (p4est, qumid,
                                        &which_tree, &quadrant_id);
    p4est_mesh_face_neighbor_init2 (&mfn, p4est, ghost, mesh,
                                    which_tree, quadrant_id); //initialisation de l'iterateur

    p4_get_midpoint(p4est,q,xyz);

    // printf("Quadrant %d de coord x = %f et y = %f \n", quadrant_id, xyz[0], xyz[1]);
    // printf("Debut de la liste des voisins du quadrant %d : \n", qumid);
    while ((q = p4est_mesh_face_neighbor_next (&mfn, &which_tree, &which_quad,
                                               &nface, &nrank)) != NULL) //boucle sur les voisins d'un quadrant
    {
      // printf("Quadrant %d \n", which_quad);
      if (which_quad==qumid) { //identification d'une maille fictive
        tmp += 1;}
    }
    // printf("Fin de la liste des voisins du quadrant %d. \n", qumid);
    // printf("------------------------\n" );
  }
  // printf("------------------------\n" );

  // printf("Nombre de mailles fictives du maillage : %d \n", tmp);

  *compteur = tmp;

  // printf("(print C) ----------------------------------------------------\n" );

}



//Iterateur pour recuperer les info d'une maille :
// remplit un tableau avec les voisins
// ordonnes par face
void p4_iterate_mesh2( p4est_t *p4est,
                       p4est_ghost_t *ghost,
                       p4est_mesh_t *mesh,
                       p4est_nodes_t *nodes,
                      int num_quad,
                      int maillage_en_espace,
                      double *area,
                      double *l_face,
                      double *x,
                      double *y,
                      double *z,
                      int nb_voisins,
                      int voisins_par_face[4][2],
                      double *h,
                      int *Nnodes,
                      int nodes_out[8])
{
  p4est_quadrant_t *q;
  p4est_mesh_face_neighbor_t mfn;
  p4est_topidx_t which_tree;
  p4est_locidx_t  quadrant_id, which_quad;
  int                 nface1, nface2, r;
  int                 nrank;
  int qumid;
  double xyz[3], tmp1;
  int tmp2;
  int level;
  double cell_area;
  double face_length;

  // p4est_mesh_t *mesh_loc;
  // p4est_ghost_t *ghost_loc;

  // ghost_loc = p4est_ghost_new(p4est, P4EST_CONNECT_FACE);
  // mesh_loc = p4est_mesh_new_ext(p4est, ghost_loc,1,1, P4EST_CONNECT_FACE);

  qumid = num_quad;
  which_tree = -1;
  q = p4est_mesh_quadrant_cumulative (p4est, qumid,
                                      &which_tree, &quadrant_id);

  p4est_mesh_face_neighbor_init2 (&mfn, p4est, ghost, mesh,
                                    which_tree, quadrant_id); //initialisation de l'iterateur

  level = q->level;

  tmp1 = 1.0;

  for (int i=0;i<level;i++){tmp1=tmp1*2;}

  tmp1 = 1.0/tmp1;

  cell_area = 1.0;

  face_length = 1.0;

  for (int i=0;i<maillage_en_espace-1;i++){
    cell_area = cell_area*tmp1;
    face_length = face_length*tmp1;}

  cell_area = cell_area*tmp1;

  // cell_area = pow(tmp1,maillage_en_espace); //aire de la cellule
  //
  // face_length = pow(tmp1,maillage_en_espace-1); //longueur d'une de ses faces

  p4_get_midpoint(p4est,q,xyz); //recupere le centre de la maille

  *x = xyz[0];
  *y = xyz[1];
  *z = xyz[2];

  *area = cell_area;
  *l_face = face_length;
  *h = tmp1;

  // printf("(print C) ----------------------------------------------------\n" );

  // printf("Quadrant %d de coord x = %f, y = %f et z = %f \n", quadrant_id,
                                                     // xyz[0], xyz[1], xyz[2]);


  // printf("area = %f, length = %f, h = %f \n", cell_area, face_length, h);

  nb_voisins = 0;

  // printf("Debut de la liste des voisins du quadrant %d : \n", qumid);

  for (int i=0; i<4; i++){
    for (int j=0; j<2; j++){
      voisins_par_face[i][j] = -1;
    }
  }

  int k=-1;
  while ((q = p4est_mesh_face_neighbor_next (&mfn, &which_tree, &which_quad,
                                             &nface2, &nrank)) != NULL) //boucle sur les voisins d'un quadrant
  {

    nb_voisins = nb_voisins + 1;

    if (nb_voisins==1){
      nface1=nface2;
      tmp2=which_quad;}

    // printf("Quadrant %d  face %d \n", which_quad, nface2);

    if (which_quad!=qumid){

      if (nface2>=0) {
        r = nface2%8;
        if ( (k+1<2)&&(r==nface1)&&(tmp2!=qumid)){k=k+1;}
        else{k=0;}
      }

      if (nface2<0) {
        r = nface2+8;
        if ( (k+1<2)&&(r==nface1)&&(tmp2!=qumid)){k=k+1;}
        else{k=0;}
      }

    }

    if (which_quad==qumid){

      if (nface2==0){r=1;k=0;}
      if (nface2==1){r=0;k=0;}
      if (nface2==2){r=3;k=0;}
      if (nface2==3){r=2;k=0;}

    }

    // printf("Quadrant %d voisin numero %d sur la face %d \n", which_quad, k, r);
    // voisins_par_face[r][k] = which_quad;
    // printf("valeur rentree dans le tableau %d\n", voisins_par_face[r][k] = which_quad);
    voisins_par_face[r][k] = which_quad;
    // printf("valeur rentree dans le tableau %d\n", voisins_par_face[r][k]);

    nface1=r;
    tmp2=which_quad;

    // printf("------------------------\n" );

  }

  // printf("Fin de la liste des voisins du quadrant %d . \n", qumid);

  // printf("Nombre de voisins de ce quadrant : %d \n", nb_voisins);

  // printf("Liste des voisins :\n" );

  // for (int i=0;i<4;i++){
  //   printf("r = %d\n", i );
  //   for (int j=0;j<2;j++){printf(" k = %d val = %d \n", j, voisins_par_face[i][j]);}
  // }

  sc_array_t *quadrants;
  p4est_tree_t *tree;
  int i;
  p4est_quadrant_t* quad_neigh;
  int* half_neigh;
  int* cell_nodes;
  int inode;

  tree = p4est_tree_array_index(p4est->trees,0); //on recupere l'abre 0
  quadrants = &tree->quadrants;

  q = p4est_mesh_quadrant_cumulative (p4est, qumid,
                                      &which_tree, &quadrant_id);

  cell_nodes = (int*) malloc(sizeof(int)*8);
  inode = 0;
  i = num_quad;

  if (mesh->quad_to_face[4*i+2]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+2]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->y<q->y)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[1]+2]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[0]+2]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i]+1;
    inode=inode+1;
  }
  if (mesh->quad_to_face[4*i]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->x<q->x)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[0]+3]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[1]+3]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i+2]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i+2]+1;
    inode=inode+1;
  }
  if (mesh->quad_to_face[4*i+3]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+3]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->y>q->y)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[0]+1]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[1]+1]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i+3]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i+3]+1;
    inode=inode+1;
  }
  if (mesh->quad_to_face[4*i+1]<0)
  {
    half_neigh=sc_array_index(mesh->quad_to_half,mesh->quad_to_quad[4*i+1]);
    quad_neigh=sc_array_index(quadrants,half_neigh[0]);
    if (quad_neigh->x>q->x)
    {
      cell_nodes[inode]=nodes->local_nodes[4*half_neigh[1]]+1;
      cell_nodes[inode+1]=nodes->local_nodes[4*half_neigh[0]]+1;
      inode=inode+2;
    }
    else
    {
      cell_nodes[inode]=nodes->local_nodes[4*i+1]+1;
      inode=inode+1;
    }
  }
  else
  {
    cell_nodes[inode]=nodes->local_nodes[4*i+1]+1;
    inode=inode+1;
  }

  *Nnodes=inode;
  // printf("Debut de la liste des noeuds du quadrant %d : \n", i);
  for (int n=0;n<inode;n++){
    nodes_out[n]=cell_nodes[n];
    // printf("noeud %d\n", cell_nodes[n]);
  }

  // printf("(print C) ----------------------------------------------------\n" );

}



void p4_get_node(p4est_t* p4est, p4est_nodes_t* nodes, int i, double* X_out, double* Y_out)
{
  double vxyz[3];
  p4est_indep_t *indep;

  indep = sc_array_index (&nodes->indep_nodes, (size_t) i);
  p4est_qcoord_to_vertex (p4est->connectivity, 0, indep->x, indep->y, vxyz);

  *X_out=vxyz[0];
  *Y_out=vxyz[1];
}



// Enregistrement des info de la foret au format vtk
void p4_save_vtk ( int num_test,
                   p4est_t *p4est,
                   p4est_connectivity_t *conn)
{
  char filename[BUFSIZ] = "";

  snprintf(filename, BUFSIZ, P4EST_STRING "_test%02d", num_test);

  p4est_geometry_t *geom;

  geom = p4est_geometry_new_connectivity(conn);

  p4est_vtk_write_file(p4est, geom, filename);

}
