#ifndef H_LISTE_DOUBLEMENT_CHAINEE
#define H_LISTE_DOUBLEMENT_CHAINEE

/* liste_doublement_chainee.h */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

/* types ======================================================================================== */

typedef enum {
        FALSE,
        TRUE
} bool_t;

/* structures =================================================================================== */

/* Cellules des listes chainees contenu dans le tableau d'haplotypes non redondants 
 * Elles constituent les genotypes qui sont expliques par les haplotypes*/
typedef struct TypeGenoExplique TypeGenoExplique;
struct TypeGenoExplique
{
    int id;
    int idHaploCompl;
    int nbIndGeno;
    TypeGenoExplique* next;
    TypeGenoExplique* prev;
};

/* Cellules du tableau d'haplotypes non redondants */
typedef struct TypeHaplo TypeHaplo;
struct TypeHaplo
{
    int id;
    float freq;
    int doublon;
    int* haplotype;
    TypeGenoExplique* teteGeno;
    TypeGenoExplique* first;
    TypeGenoExplique* last;
};

/* Cellules des listes chainees contenu dans le tableau des genotypes non redondants 
 * Elles constituent les paires d'haplotypes qui expliques les genotypes */
typedef struct TypePaireHaplo TypePaireHaplo;
struct TypePaireHaplo
{
    int idHaplo1;
    int idHaplo2;
    TypePaireHaplo* next;
    TypePaireHaplo* prev;
};

/* Cellules du tableau de genotypes non redondants */
typedef struct TypeGeno TypeGeno;
struct TypeGeno
{
    int id;
    int nbLociAmbigu;
    int nbHaplo;
    int nbIdentique;
    int doublon;
    double proba;
    double probaPrec;
    int* genotype;
    TypeHaplo* matriceHaplo;
    TypePaireHaplo* tetePaireHaplo;
    TypePaireHaplo* first;
    TypePaireHaplo* last;
};

/* point d'entree =============================================================================== */

TypeHaplo* initialiser_liste_geno(int tailleGeno);
TypeGeno* initialiser_liste_haplo(int tailleGeno);
void ajout_queue_geno(TypeHaplo* liste, int id, int idHaplo);
void ajout_queue_paire_haplo(TypeGeno* liste, int id1, int id2);
void suppression_liste_geno(TypeHaplo** ptrListe);
void suppression_liste_paire_haplo(TypeGeno** ptrListe);
TypeGenoExplique* recherche_id_geno(TypeHaplo* liste, int id);
TypePaireHaplo* recherche_id_paire_haplo(TypeGeno* liste, int id);
void affichage_liste_geno(TypeHaplo* liste);
void affichage_liste_paire_haplo(TypeGeno* liste);
int* id_pos(TypeGeno* liste, int position);


#endif /* H_LISTE_DOUBLEMENT_CHAINEE */
