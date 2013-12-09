#ifndef H_FONCTIONS_VERIFICATION
#define H_FONCTIONS_VERIFICATION

/* fonctions.h */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

/* types ======================================================================================== */

typedef enum {
        FALSE,
        TRUE
} bool_t;

/* constantes =================================================================================== */

#define PARAM "fichiers/parametres.txt"
#define THEORIQUE "fichiers/geno_et_haplo.txt"
#define RESULTATS "fichiers/resultats_haplo.txt"
#define LIRE_TAILLE 50

/* point d'entree =============================================================================== */

void lire_geno_haplos_theo(int* geno, int* haplo1, int* haplo2, FILE* theorique, int tailleGeno);
void lire_geno_haplos_res(int* geno, int* haplo1, int* haplo2, FILE* resultats, int tailleGeno);
bool_t verif_combinaison_haplo(int* geno, int* haplo1, int* haplo2, int tailleGeno);
bool_t comparaison_haplo(int* haplo1, int* haplo2, int tailleGeno);
bool_t comparaison_geno(int* geno1, int* geno2, int tailleGeno);
bool_t verif_doublon(int* seq1, int* seq2, int tailleGeno);

#endif /* H_FONCTIONS_VERIFICATION */
