/* verification.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/fonctions_verification.h"


/* point d'entree =============================================================================== */

int main()
{
	/* ====== Declarations ====== */
    
    int i = 0;
    int nbErreurs = 0;      /* Nombre de mauvaises inférences d'haplotypes */
    int nbIndiv = 0;
    int tailleGeno = 0;
    int nbLoci = 0;
    int* theoGeno = NULL;   /* Contient le genotype reel */
    int* theoHaplo1 = NULL; /* Contient l'haplotype1 reel */
    int* theoHaplo2 = NULL; /* Contient l'haplotype2 reel */
    int* resGeno = NULL;    /* Contient le genotype */
    int* resHaplo1 = NULL;  /* Contient l'haplotype1 inféré */
    int* resHaplo2 = NULL;  /* Contient l'haplotype2 inféré */
    FILE* fichierParam = NULL;
    FILE* theorique = NULL;
    FILE* resultats = NULL;
    
    /* ========== Code ========== */
    
    /* Recuperation des parametres passes au logiciel de generation de genotypes */
    fichierParam = fopen(PARAM, "r");
    if (fichierParam != NULL)
    {
        fscanf(fichierParam, "%d", &nbIndiv);
        fscanf(fichierParam, "%d", &tailleGeno);
        fscanf(fichierParam, "%d", &nbLoci);
        
        /* Allocation de la memoire pour le genotype du fichier theorique et sa paire d'haplotypes */ 
        theoGeno = malloc(sizeof(int) * tailleGeno);
        theoHaplo1 = malloc(sizeof(int) * tailleGeno);
        theoHaplo2 = malloc(sizeof(int) * tailleGeno);
        
        /* Allocation de la memoire pour le genotype du fichier resultats et sa paire d'haplotypes */ 
        resGeno = malloc(sizeof(int) * tailleGeno);
        resHaplo1 = malloc(sizeof(int) * tailleGeno);
        resHaplo2 = malloc(sizeof(int) * tailleGeno);
    }
    else
    {
        fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n", PARAM);
        exit(1);
    }
    fclose(fichierParam);
    
    theorique = fopen(THEORIQUE, "r");
    if ( theorique != NULL)
    {
        resultats = fopen(RESULTATS, "r");
        if ( resultats != NULL)
        {
            for (i=0; i < nbIndiv; i++)
            {
                /* Recupere le genotype reel et ses deux haplotypes */
                lire_geno_haplos_theo(theoGeno, theoHaplo1, theoHaplo2, theorique, tailleGeno);
                if (!verif_combinaison_haplo(theoGeno, theoHaplo1, theoHaplo2, tailleGeno))
                    fprintf(stderr, "Une paire d'haplotype ne forme pas correctement son genotype.\n");
                
                /* Recupere le genotype et ses deux haplotypes inferes*/
                lire_geno_haplos_res(resGeno, resHaplo1, resHaplo2, resultats, tailleGeno);
                if (!verif_combinaison_haplo(resGeno, resHaplo1, resHaplo2, tailleGeno))
                    fprintf(stderr, "Une paire d'haplotype ne forme pas correctement son genotype.\n");
                
                /* Teste l'egalite entre les genotypes */
                if (!comparaison_geno(theoGeno, resGeno, tailleGeno))
                    fprintf(stderr, "Les genotypes ne sont pas identiques.\n");
                
                /* Verifie si les haplotypes inferes correspondent aux haplotypes reels */
                if (verif_doublon(theoHaplo1, resHaplo1, tailleGeno) || 
                    verif_doublon(theoHaplo1, resHaplo2, tailleGeno))
                {
                    printf("Indiv n°%d : OK\n",i);
                }
                else
                {
                    printf("Indiv n°%d : Erreur !\n",i);
                    nbErreurs++;
                }
            }
                printf("Nombre d'haplotypes correctement inférés: %d/%d\n", nbIndiv-nbErreurs,nbIndiv);
                printf("Pourcentage des haplotypes bien inférés : %.2f%%\n",((nbIndiv-nbErreurs)*1.0/(nbIndiv*1.0)*100));
                printf("Nombre de loci ambigus : %d\n", nbLoci);
        }
        else
        {
            fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n", RESULTATS);
            exit(1);
        }
        fclose(resultats);
    }
    else
    {
        fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n", THEORIQUE);
        exit(1);
    }
    fclose(theorique);
    
    free(theoGeno);
    free(theoHaplo1);
    free(theoHaplo2);
    free(resGeno);
    free(resHaplo1);
    free(resHaplo2);
    
    return 0;
}
