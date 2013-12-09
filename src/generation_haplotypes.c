/* generation_haplotypes.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/generation_haplotypes.h"

/* fonctions tests =========================================================================== */

/* Affichage des haplytpes */
void affichage_haplotypes(TypeGeno geno)
{
    int i;
    int j;
    for(i = 0 ; i < geno.nbHaplo ; i++ )
    {
        printf("%d : ",geno.matriceHaplo[i].id);
        for(j = 0 ; j < TAILLE_GENO ; j++)
        {
           printf("%d",geno.matriceHaplo[i].haplotype[j]); 
        }
        printf("\n");
    }
}

/* Affichage des genotypes (int *) */
void affichage_genotype(TypeGeno geno)
{
    int i;
    printf("%d : ",geno.id);
    for(i = 0 ; i < TAILLE_GENO ; i++ )
    {
        printf("%d",geno.genotype[i]);
    }
    printf("\n");
}

/* Affichage des hapltypes (int *) */
void affichage_haplotype(TypeHaplo haplo)
{
    int j;
        printf("%d : ",haplo.id);
        for(j = 0 ; j < TAILLE_GENO ; j++)
        {
           printf("%d",haplo.haplotype[j]); 
        }
        printf("\n");
}

/* fonctions privees ========================================================================= */

/* Compte le nombre de loci ambigu pour un genotype donnee */
static int compte_nombre_loci_ambigu(TypeGeno geno)
{
    int i;
    int count = 0;
    for(i = 0 ; i < TAILLE_GENO ; i++ )
    {
        if(geno.genotype[i] == 1)
        {
            count++;
        }
    }
    return count;
}

/* Rempli la matrice d'haplotypes en colonne */
static void generation_haplo_possibles(TypeHaplo* matrice, int lociMax, 
                                        int countLoci, int colonne, int n)
{
    int i;
    int modulo;
    int interval;
    if (n == 2)
    {
        for( i=0 ; i < lociMax ; i++){
            matrice[i].haplotype[colonne] = 1;
        }
    }
    else if (n == 0)
    {
        for( i=0 ; i < lociMax ; i++){
            matrice[i].haplotype[colonne] = 0;
        }
    }
    else /* Lorsque le locus est ambigu */
    {
        modulo = pow(2,countLoci);
        interval = lociMax / modulo;

        for( i=0 ; i < lociMax ; i++){
            matrice[i].haplotype[colonne] = ((i/interval)%2);
        }
    } 
}

#if 0
static void affichage_haplotypes(TypeGeno geno)
{
    int i;
    int j;
    for(i = 0 ; i < geno.nbHaplo ; i++ )
    {
        printf("%d : ",geno.matriceHaplo[i].id);
        for(j = 0 ; j < TAILLE_GENO ; j++)
        {
           printf("%d",geno.matriceHaplo[i].haplotype[j]); 
        }
        printf("\n");
    }
}
#endif



/* Verifie chiffre par chiffre si les 2 haplotypes sont egaux */
int verif_doublon(int* seq1, int* seq2)
{
    int verifPresenceDoublon = 1;
    int i = 0;
    while ((verifPresenceDoublon == 1)&&(i < TAILLE_GENO))
    {
        if(seq1[i] != seq2[i])
        {
            verifPresenceDoublon = 0;
            return verifPresenceDoublon;
        }
        i++;
    }
    /*** TEST 3 ***/
    /*** TEST 4 ***/
    return verifPresenceDoublon;
}

/* fonctions publiques ======================================================================= */

int* alloue_memoire()
{
    int* pointeur = NULL;
    pointeur = malloc(TAILLE_GENO * sizeof(int));
    if (pointeur == NULL)
    {
        fprintf(stderr,"Erreur d'allocation mémoire.\n");
        exit(1);
    }
    return pointeur;
}

/* Suppression du retour chariot donne par fgets */
int lire(char* chaine, int longueur, FILE* fichier)
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

/* Fonction de comparaison utilisee par la fonction qsort() */
int compare (const void *a, const void *b)
{
    double **pa = (double**)a;
    double **pb = (double**)b;
    double diff = (*pb)[1] - (*pa)[1];
    double ret = 0;
    if (diff > 0)
    {
            ret = 1;
    }
    else if (diff < 0)
    {
            ret = -1;
    }
    else
    {
            ret = 0;
    }
    return ret;
}

/* Initialisation des differents genotypes */
int initialisation_geno(TypeGeno* geno, int id)
{   
    int i, countLoci;           /* Compteur */
    
    /* Initialisation des parametres du genotype */
    geno->nbLociAmbigu = compte_nombre_loci_ambigu(*geno);
    geno->nbHaplo = pow(2,geno->nbLociAmbigu);
    geno->nbIdentique = 1;
    geno->doublon = 0;
    geno->proba = 0;
    countLoci = 0;

    /* Initialisation de la matriceHaplo */
    geno->matriceHaplo = malloc(sizeof(TypeHaplo) * (geno->nbHaplo));
    for(i=0 ; i < geno->nbHaplo ; i++){
        geno->matriceHaplo[i].haplotype = alloue_memoire();
        geno->matriceHaplo[i].id = id;
        geno->matriceHaplo[i].doublon = 0;
        id ++;
    }   

    /* Generation de tous les haplotypes differents possibles */
    for(i = 0 ; i < TAILLE_GENO ; i++ )
    {
        /* Permet de savoir a quelle position ambigu on se trouve */
        if(geno->genotype[i] == 1){
            countLoci ++;
        }
        generation_haplo_possibles(geno->matriceHaplo,geno->nbHaplo,countLoci,i,geno->genotype[i]);
    }
    #if 1
    printf("*********************\n");
    affichage_genotype(*geno);
    printf("---------------------\n");
    affichage_haplotypes(*geno);
    #endif
    /*** TEST 5 ***/

    return id;
}

/* Recherche de doublon(s) de genotypes */
void recherche_genotype_doublon(TypeGeno* geno1, TypeGeno* geno2)
{
    int verif = 0;
    /* Regarde parmi tous les haplotypes generes dans un second geno2 */
    verif=0;
    verif=verif_doublon(geno1->genotype,geno2->genotype);
    if (verif==1)
    {
        geno2->id = geno1->id;
        geno2->doublon = 1;
        geno1->nbIdentique = geno1->nbIdentique +1;
        /*printf("Genotype %d est redondant !\n",geno1->id);*/
        /*geno1->nbIdentique=geno1->nbIdentique+1;*/
    }
}

/* Modification du nombre de genotype identique au gentoype etudie */
void modification_nb_geno_identique(int id,TypeGeno* geno)
{
    int i;
    bool_t trouve = FALSE;
    while (trouve == FALSE)
    {
        if(id == geno[i].id)
        {
            geno[i].nbIdentique = +1;
            trouve = TRUE;
        }
        i++;
    }
}

/* Recherche de doublon(s) d'haplotypes */
void recherche_haplotype_doublon(TypeGeno* geno1, TypeGeno* geno2)
{
    /* i correspond à un haplotype de g1 */
    /* j correspond à un haplotype de g2 */
    int i, j; 
    int verif = 0;
    for(i=0 ; i<geno1->nbHaplo ; i++)
    {
        /* Regarde parmi tous les haplotypes generes dans un premier geno1 */
        for(j=0 ; j<geno2->nbHaplo ; j++)
        {
            /* Regarde parmi tous les haplotypes generes dans un second geno2 */
            verif=0;
            verif=verif_doublon(geno1->matriceHaplo[i].haplotype,geno2->matriceHaplo[j].haplotype);
            if (verif==1)
            {
                geno2->matriceHaplo[j].id = geno1->matriceHaplo[i].id;
                geno2->matriceHaplo[j].doublon=1;
            }
        }
    }
}

/* Calcul le nombre d'haplotypes non redondants */
int calcul_nb_haplo_non_redondant(TypeGeno* geno)
{
    int count = 0;
    int i,j; /*Compteur*/
    for (i=0 ; i<NB_INDIV ; i++)
    {
        for(j=0 ; j<geno[i].nbHaplo ; j++)
        {
            if (geno[i].matriceHaplo[j].doublon == 0){
                count ++;
            }
        }
    }
    return count;
}

/* Calcul le nombre d'haplotypes non redondants */
int calcul_nb_geno_non_redondant(TypeGeno* geno)
{
    int count = 0;
    int i; /*Compteur*/
    for (i=0 ; i<NB_INDIV ; i++)
    {
        if (geno[i].doublon == 0){
            count ++;
        }
    }
    return count;
}

/* Fonction de recherche d'haplotype complementaire par indice */
int recherche_haplo_complementaire(TypeGeno geno, int indice)
{
    int indiceComplementaire;
    if(geno.nbHaplo > 1)
    {
        indiceComplementaire = ((geno.nbHaplo-1) - indice); /*miroir*/ 
    }
    else
    {
        indiceComplementaire = indice;
    }
    
    return indiceComplementaire;
}

/* Allocation mémoire pour le tableau de frequence d'haplotype */
double** allouer_memoire_tableau_2d(int nb)
{
    int i;
    double** pointeur = NULL;
    pointeur = malloc(sizeof(double*) * nb);
    for(i=0 ; i<nb ; i++)
    {
        pointeur[i] = malloc(sizeof(double) * 3);
    }
    return pointeur;
}

bool_t compare_tableaux_entiers(int *tab1, int *tab2, int taille)
{
    int i = 0;
    
    for (i=0; i < taille; i++) 
    {
        if (tab1[i] != tab2[i])
        {
            return FALSE;
        }
    }
    return TRUE;
}