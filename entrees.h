static const double longueurs[2] = {10, 10}; // Ly, Lx (en cm)

static const int face_entree = 0; // 0=face x ; 1=face y ; 3=face z

static const double coord_centre_faisceau = 5; //coordonnees du centre du faisceau sur la face d entree : 1 coord de 2D, 2 coord en 3D

static const double largeur_faisceau = 2.; // largeur du faisceau

static const double marge = 0.5; // on rajoute une marge de mailles fines autour du faisceau

static const int nb_zones = 1; // nombre de zones definies pour changer les densites

static const double xmin[1] = {0.0}; // xmin delimite l extremite gauche de chaque zone definie

static const double xmax[1] = {10}; // xmax delimite l extremite droite de chaque zone definie

static const double ymin[1] = {0.0}; // ymin delimite l extremite en bas de chaque zone definie

static const double ymax[1] = {10}; // ymax delimite l extremite en haut de chaque zone definie

static const int niveau_max_p4est[1] = {9}; // niveau de raffinement maximum de chaque zone

static const int niveau_min_p4est[1] = {3}; // niveau de raffinement minimum de chaque zone
