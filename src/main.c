/* main.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "../inc/fonctions.h"

#define GENO_HAPLO "./fichiers/geno_et_haplo.txt"
#define GENOTYPES "./fichiers/genotypes.txt"

/* point d'entree =============================================================================== */

int main(int argc,char* argv[])
{
    
	/* ====== Declarations ====== */
    
	int i, j = 0; /* compteurs */
    int nbNonRedondant;
	TypeHaplo* haplo = NULL;
    TypeHaplo* haploNonRedondant = NULL;
	TypeGeno* geno = NULL;
    FILE* fichier = NULL;
    
    /* ========== Code ========== */

	if (argc != 4)
    {
	    printf("Le nombre d'argument est incorrect !\n");
	    printf("la commande doit etre du type :\n");
	    printf("./main nb_individu(entier) taille_genome(entier) nb_de_loci(entier)\n");
	    exit(1);	
	}
	for (i=1 ; i < argc ; i++)
    {
	    if (atoi(argv[i]) == 0)
        {
	        printf("Impossible d'avoir un paramêtre nul ou contenant un caractère\n");
	        exit(1);
	    }
	}
    
    srand (time(NULL));
	
    /* Suppression des fichiers de genotypes si presents */
	test_existence(GENO_HAPLO);
	test_existence(GENOTYPES);
    
    /* Recuperation des parametres entres au clavier */
	nbIndiv = atoi(argv[1]);
	tailleGeno = atoi(argv[2]);
	nbLoci = atoi(argv[3]);

	/* Verifications des parametres entres au clavier */
    if (tailleGeno <= nbLoci)
    {
        printf("Le nombre de loci ambigu maximum doit etre inferieur à la taille du genotype.\n");
        exit(1);
    }

    /* 
     * Determination du nombre d'haplotypes que l'on veut generer 
     * pour creer la liste d'haplotypes qui servira a faire creer les genotypes
     * Le premier nombre peut etre modifie.
     */
	nbHaplo = 1.5*nbIndiv;
    
    /* Allocation de memoire pour les tableaux de genotypes et d'haplotypes */
    haplo = malloc(NB_HAPLO * sizeof(TypeHaplo));
    geno = malloc(NB_INDIV * sizeof(TypeGeno));

	/* Initialisation du tableau d'haplotypes */
	for (i=0 ; i < NB_HAPLO ; i++)
    {
		haplo[i].id=i;
        haplo[i].doublon=0;
        initialiser_haplotypes(&haplo[i]);
	}

	/*
     * Parcours la liste d'haplotypes pour verifier l'absence de doublon.
     * si un doublon est present, le premier croise obtient un VRAI
     * dans l'element doublon.
	 */
	recherche_haplotype_doublon(haplo);
    nbNonRedondant = compte_nombre_doublon(haplo);
    /* Creation de la liste d'haplotypes non redondants */
    haploNonRedondant = lister_haplo_non_redondant(nbNonRedondant, haplo);
    
    /* Affichage des résultats */
	printf("nb de genotypes générés : %d \n",NB_INDIV);
	printf("nb d'haplotypes générés : %d \n",NB_HAPLO);
	printf("taille du genome : %d \n",TAILLE_GENO);
	printf("nb de loci differents tolere : %d \n",NB_LOCI);
    printf("nb haplo non redondant : %d \n",nbNonRedondant);
    
    /* Stockage des parametres dans un fichier */
    fichier = fopen(PARAM, "w");
    if ( fichier != NULL)
    {
        fprintf(fichier, "%d\t%d\t%d\n", NB_INDIV, TAILLE_GENO, NB_LOCI);
    }
    else
    {
        fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n", PARAM);
        exit(1);
    }
    fclose(fichier);
    
    /* 
     * Initialisation des genotypes en tirant aleatoirement des haplotypes 
     * issues de la liste cree precedemment. 
     */
	for (i = 0 ; i < NB_INDIV ; i++)
    {
		geno[i].id=i;
		if (TAILLE_GENO < 11)
		{
			/* Tirage des 2 haplotypes dans le pool d'haplotypes */
			initialiser_genotypes_petite_taille(&geno[i],haploNonRedondant,nbNonRedondant);
		}
		else
		{
			/* Tirage d'1 haplotype dans le pool d'haplotypes et generation aleatoire du deuxieme */
			initialiser_genotypes(&geno[i],haploNonRedondant,nbNonRedondant);
		}		
	}
	
	/* Liberation de l'espace memoire alloue au tableau de genotypes */
	#if 0
	for (i=0; i < NB_INDIV; i++)
	{
		free(geno[i].haplo1.haplotype);
		free(geno[i].haplo2.haplotype);
	}
	#endif
    free(geno);
	
	/* Liberation de l'espace memoire alloue au tableau d'haplotypes */
	for (i=0; i < NB_HAPLO; i++)
	{
		for(j=0; j < TAILLE_GENO; j++)
		{
			free(haplo[i].haplotype);
			haplo[i].haplotype = NULL;
		}
	}
	free(haplo);
	
	free(haploNonRedondant);
	
    return 0;
}
