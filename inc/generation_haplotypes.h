#ifndef H_GENERATION_HAPLOTYPES
#define H_GENERATION_HAPLOTYPES

/* generation_haplotypes.h */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include "../inc/liste_doublement_chainee.h"

/* variables globales =========================================================================== */

double nbIndiv;
int tailleGeno;

/* constantes =================================================================================== */

#define PARAM "./fichiers/parametres.txt"
#define TAILLE_GENO tailleGeno /* longueur du fragment d'ADN observe */
#define NB_INDIV nbIndiv       /* nombre d'individus sur lesquels les genotypes seront observes */

/* fonctions tests ============================================================================== */

void affichage_genotype(TypeGeno geno);
void affichage_haplotype(TypeHaplo haplo);
int verif_doublon(int* seq1, int* seq2);

/* point d'entree =============================================================================== */

int lire(char* chaine, int longueur, FILE* fichier);
int compare (const void *a, const void *b);
int initialisation_geno(TypeGeno* geno, int id);
void recherche_genotype_doublon(TypeGeno* geno1, TypeGeno* geno2);
void modification_nb_geno_identique(int id,TypeGeno* geno);
void recherche_haplotype_doublon(TypeGeno* geno1, TypeGeno* geno2);
int calcul_nb_haplo_non_redondant(TypeGeno* geno);
int calcul_nb_geno_non_redondant(TypeGeno* geno);
int recherche_haplo_complementaire(TypeGeno geno, int indice);
double** allouer_memoire_tableau_2d(int nb);
double recherche_frequence(int id,double** tableau);


#endif /* H_GENERATION_HAPLOTYPES */
