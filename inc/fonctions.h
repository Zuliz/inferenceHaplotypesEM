#ifndef H_FONCTIONS
#define H_FONCTIONS

/* fonctions.h */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

/* variables globales =========================================================================== */

int nbIndiv;
int tailleGeno;
int nbHaplo;
int nbLoci;

/* constantes =================================================================================== */

#define TAILLE_GENO tailleGeno /* longueur du fragment d'ADN observe */
#define NB_INDIV nbIndiv /* nombre d'individus sur lesquels les genotypes seront observe */
#define NB_HAPLO nbHaplo /* Creation d'un peu moins d'haplotype que de genotype */
#define NB_LOCI nbLoci   /* nombre maximum de loci ambigus */ 
#define PARAM "./fichiers/parametres.txt"
#define GENO_HAPLO "./fichiers/geno_et_haplo.txt"
#define GENOTYPES "./fichiers/genotypes.txt"

/* types ======================================================================================== */

typedef enum
{
	FAUX,
	VRAI
}TypeBool;

/* structures =================================================================================== */

typedef struct TypeHaplo
{
	int id;
	int* haplotype;
	TypeBool doublon;
}TypeHaplo;

typedef struct TypeGeno
{
	int id;
	TypeHaplo haplo1;
	TypeHaplo haplo2;
	int* genotype;
}TypeGeno;

/* point d'entree =============================================================================== */

void test_existence(const char* nomFichier);
int* alloue_memoire();
void initialiser_haplotypes(TypeHaplo* haplo);
/*void afficher_haplotypes(TypeHaplo haplo);*/
void afficher_genotypes(TypeGeno geno);
void afficher_haplotypes(TypeHaplo haplo);
TypeBool verification_presence_doublon(TypeHaplo haplo1, TypeHaplo haplo2);
TypeBool verification_nombre_loci(TypeHaplo haplo1, TypeHaplo haplo2);
void initialiser_genotypes(TypeGeno* adrGeno, TypeHaplo haplo[], int nbHaploNonRedondant);
void recherche_haplotype_doublon(TypeHaplo haplo[]);
int compte_nombre_doublon(TypeHaplo haplo[]);
TypeHaplo* lister_haplo_non_redondant(int compteur, TypeHaplo haplo[]);

#endif /* H_FONCTIONS */
