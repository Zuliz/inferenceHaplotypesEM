/* fonctions_verification.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/fonctions_verification.h"

/* fonctions privees ============================================================================ */

/* Suppression du retour chariot donne par fgets */
static int lire(char* chaine, int longueur, FILE* fichier)
{
    char* positionEntree = NULL;
  
    if (fgets(chaine, longueur, fichier) != NULL)
    {
        positionEntree = strchr(chaine, '\n');
        if (positionEntree != NULL)
        {
            *positionEntree = '\0';
        }
        return 0;
    }
    else
    {
        return 1;
    }
}

/* fonctions publiques ========================================================================== */

/* Remplissage d'un tableau genotype et deux tableaux haplotypes a partir du fichier theorique */
void lire_geno_haplos_theo(int* geno, int* haplo1, int* haplo2, FILE* theorique, int tailleGeno)
{
    /* Variables Locales */
    int i = 0;
    char* chaine = NULL;
    char* sousChaine = NULL;
    
    chaine = (char*)malloc(sizeof(char) * LIRE_TAILLE);
    
    /* Recuperation du genotype */
    if(lire(chaine, LIRE_TAILLE, theorique) == 0)
    {
        sousChaine = strchr(chaine, '\t');
        sousChaine = sousChaine + 1;
        if(sousChaine != NULL)
        {
            for(i=0; i < tailleGeno; i++)
            {
                geno[i] = sousChaine[i] - '0';
            }
        }
    
    }
    
    /* Recuperation du premier haplotype */
    if(lire(chaine, LIRE_TAILLE, theorique) == 0)
    {
        sousChaine = strrchr(chaine, ' ');
        sousChaine = sousChaine + 1;
        if(sousChaine != NULL)
        {
            for(i=0; i < tailleGeno; i++)
            {
                haplo1[i] = sousChaine[i] - '0';
            }
        }
    }
    
    /* Recuperation du deuxieme haplotype */
    if(lire(chaine, LIRE_TAILLE, theorique) == 0)
    {
        sousChaine = strrchr(chaine, ' ');
        sousChaine = sousChaine + 1;
        if(sousChaine != NULL)
        {
            for(i=0; i < tailleGeno; i++)
            {
                haplo2[i] = sousChaine[i] - '0';
            }
        }
    }
    
    free(chaine);
}

/* Remplissage d'un tableau genotype et deux tableaux haplotypes a partir du fichier resultats*/
void lire_geno_haplos_res(int* geno, int* haplo1, int* haplo2, FILE* resultats, int tailleGeno)
{
    /* Variables Locales */
    int i = 0;
    char* chaine = NULL;
    char* sousChaine = NULL;
    
    chaine = (char*)malloc(sizeof(char) * LIRE_TAILLE);
    
    /* Recuperation du genotype */
    if(lire(chaine, LIRE_TAILLE, resultats) == 0)
    {
        sousChaine = strrchr(chaine, ' ');
        sousChaine = sousChaine + 1;
        if(sousChaine != NULL)
        {
            for(i=0; i < tailleGeno; i++)
            {
                geno[i] = sousChaine[i] - '0';
            }
        }
    
    }
    
    /* Recuperation du premier haplotype */
    if(lire(chaine, LIRE_TAILLE, resultats) == 0)
    {
        sousChaine = strrchr(chaine, ' ');
        sousChaine = sousChaine + 1;
        if(sousChaine != NULL)
        {
            for(i=0; i < tailleGeno; i++)
            {
                haplo1[i] = sousChaine[i] - '0';
            }
        }
    }
    
    /* Recuperation du deuxieme haplotype */
    if(lire(chaine, LIRE_TAILLE, resultats) == 0)
    {
        sousChaine = strrchr(chaine, ' ');
        sousChaine = sousChaine + 1;
        if(sousChaine != NULL)
        {
            for(i=0; i < tailleGeno; i++)
            {
                haplo2[i] = sousChaine[i] - '0';
            }
        }
    }
    
    free(chaine);
}

bool_t verif_combinaison_haplo(int* geno, int* haplo1, int* haplo2, int tailleGeno)
{
    /* Varibles locales */
    int i = 0;
    
    for(i=0; i < tailleGeno; i++)
    {
        if (geno[i] != (haplo1[i] + haplo2[i]))
            return FALSE;
    }
    
    return TRUE;    
}

/* Verifie l'egalite entre deux haplotypes */
bool_t comparaison_haplo(int* haplo1, int* haplo2, int tailleGeno)
{
    /* Varibles locales */
    int i = 0;
    
    for(i=0; i < tailleGeno; i++)
    {
        if (haplo1[i] != haplo2[i])
            return FALSE;
    }
    
    return TRUE;    
}

/* Verifie l'egalite entre deux genotypes */
bool_t comparaison_geno(int* geno1, int* geno2, int tailleGeno)
{
    /* Varibles locales */
    int i = 0;
    
    for(i=0; i < tailleGeno; i++)
    {
        if (geno1[i] != geno2[i])
            return FALSE;
    }
    
    return TRUE;    
}

/* Verifie chiffre par chiffre si les 2 haplotypes sont egaux */
bool_t verif_doublon(int* seq1, int* seq2, int tailleGeno)
{
    int verifPresenceDoublon = TRUE;
    int i = 0;
    while ((verifPresenceDoublon == TRUE)&&(i < tailleGeno))
    {
        if(seq1[i] != seq2[i])
        {
            verifPresenceDoublon = FALSE;
            return verifPresenceDoublon;
        }
        i++;
    }
    /*** TEST 3 ***/
    /*** TEST 4 ***/
    return verifPresenceDoublon;
}
