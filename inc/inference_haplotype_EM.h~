#ifndef H_INFERENCE_HAPLOTYPE_EM
#define H_INFERENCE_HAPLOTYPE_EM

/* inference_haplotype_EM.h */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include "../inc/generation_haplotypes.h"
#include "../inc/liste_doublement_chainee.h"

/* constantes =================================================================================== */

#define RESULTATS_FREQ "fichiers/resultats_freq.txt"
#define RESULTATS_HAPLO "fichiers/resultats_haplo.txt"

/* point d'entree =============================================================================== */

void inference_haplotype_em(double seuil, int nbGeno, int nbHaplo, int nbEtapeMax,
                             double** tabFreq, TypeGeno** tabGeno, TypeHaplo** tabHaplo);
int* recherche_meilleure_paire_haplo(TypeGeno* liste, double** tabFreqHaplo, int taille);
void ecriture_paire_haplotype(TypeGeno** tabGenoNR, double** tabFreqHaplo, int nbHaploNonRedondant, 
                              TypeHaplo** tabHaploNR, FILE* resultats_haplo, int i, int k);

                          
#endif /* H_INFERENCE_HAPLOTYPE_EM */
