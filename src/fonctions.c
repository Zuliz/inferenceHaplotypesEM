/* fonctions.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/fonctions.h"

/* fonctions privees ============================================================================ */

/* Generation aleatoire de "0" et de "1" */
static int random_binaire()
{
	const int MIN = 0, MAX = 1;
	int nombreBinaire;
	nombreBinaire = (rand() % (MAX + 1 - MIN)) + MIN; 
	return nombreBinaire;
}

/* Fonction de tirage au sort d'haplotype */
static int tirage_au_sort(int nb)
{
	const int MIN = 0, MAX = (nb-1);
	int resultat;
	resultat = (rand() % (MAX + 1 - MIN)) + MIN; 
	return resultat;
}

/* Fonction d'ecriture d'un genotype dans un fichier */
static void ecrire_genotypes(TypeGeno geno, FILE* fic)
{
	int colonne;
	fprintf(fic,"/ind %d geno\t",geno.id);
	for (colonne = 0 ; colonne < TAILLE_GENO ; colonne++)
	{
		fprintf(fic,"%d",geno.genotype[colonne]);
	}
	fprintf(fic,"\n");
}

/* Fonction d'ecriture d'haplotype dans un fichier */
static void ecrire_haplotypes(TypeHaplo haplo, FILE* fic)
{
	int colonne;
	for (colonne = 0 ; colonne < TAILLE_GENO ; colonne++)
	{
		fprintf(fic,"%d",haplo.haplotype[colonne]);
	}
	fprintf(fic,"\n");
	/*(test)printf("\tdoublon=%d\n",haplo.doublon);*/
}

/* Fonction de creation de genotypes */
static void creation_genotypes(TypeGeno* adrGeno)
{
	int i;
	adrGeno->genotype = alloue_memoire();
	for (i = 0 ; i < TAILLE_GENO ; i++)
	{
		if (adrGeno->haplo1.haplotype[i] == 0)
		{
			if (adrGeno->haplo2.haplotype[i] == 0)
			{
				adrGeno->genotype[i]=0;
			}
			else
			{
				adrGeno->genotype[i]=1;
			}
		}
		else
		{
			if (adrGeno->haplo2.haplotype[i] == 0)
			{
				adrGeno->genotype[i]=1;
			}
			else
			{
				adrGeno->genotype[i]=2;
			}
		}
	}
}

/* Genere un fichier compose des genotypes et de leurs haplotypes reels */
static void creation_fichier_geno_haplo(TypeGeno* adrGeno, TypeHaplo haplo1, TypeHaplo haplo2)
{
	FILE* fichier=NULL;
	fichier = fopen(GENO_HAPLO,"a");

	if (fichier != NULL)
	{
		ecrire_genotypes(*adrGeno,fichier);
		fprintf(fichier,"/ind %d real haplo1 ",adrGeno->id);
		ecrire_haplotypes(haplo1,fichier);
		fprintf(fichier,"/ind %d real haplo2 ",adrGeno->id);
		ecrire_haplotypes(haplo2,fichier);
	}
	fclose(fichier);
}

/* Genere un fichier compose des genotypes */
static void creation_fichier_geno(TypeGeno* adrGeno)
{
	FILE* fichier1=NULL;
	fichier1 = fopen(GENOTYPES ,"a");

	if (fichier1 != NULL)
	{
        ecrire_genotypes(*adrGeno,fichier1);
	}
	fclose(fichier1);
}

/* fonctions publiques ========================================================================== */

/* Supprime les fichiers genotypes precedemment genere */
void test_existence(const char* nomFichier)
{
    FILE* fic = fopen(nomFichier, "r");
    int ret;

    if (fic != NULL)
    {
        fclose(fic);
        ret = remove(nomFichier);
        if (ret != 0)
        {
            printf("Erreur, un fichier n'a pu être supprimer.\n");
            exit(1);
        }
        else
        {
    		printf("Le fichier %s existe déjà. L'ancien est supprimé.\n",nomFichier);
    	}
    }
}

/* Fonction d'allocation de memoire */
int* alloue_memoire()
{
	int* pointeur = NULL;
	pointeur = malloc(TAILLE_GENO * sizeof(int));
    if (pointeur == NULL)
    {
        fprintf(stderr,"Erreur d'allocation mémoire");
        exit(1);
    }
	return pointeur;
}

/* Generation de l'haplotype */
void initialiser_haplotypes(TypeHaplo* adrHaplo){
	int colonne;
	adrHaplo->haplotype = alloue_memoire();
	for (colonne = 0 ; colonne < TAILLE_GENO ; colonne++)
	{
		adrHaplo->haplotype[colonne] = random_binaire();
	}
}

/* Fonction d'epuration des doublons */
TypeBool verification_presence_doublon(TypeHaplo haplo1, TypeHaplo haplo2){
	TypeBool verifPresenceDoublon = VRAI;
	int i = 0;
	while ((verifPresenceDoublon == VRAI)&&(i < TAILLE_GENO))
	{
		if(haplo1.haplotype[i] != haplo2.haplotype[i])
		{
			verifPresenceDoublon = FAUX;
			return verifPresenceDoublon;
		}
		i++;
	}
	/* (test)Affichage de doublon lorsqu'il y en a un de present */
	/* (test)printf("doublon-%d\n",verifPresenceDoublon);*/
	return verifPresenceDoublon;
}

/* (test)Affichage d'un tableau 2D d'haplotype */
void afficher_haplotypes(TypeHaplo haplo){
	int colonne;
	for (colonne = 0 ; colonne < TAILLE_GENO ; colonne++)
	{
		printf("%d",haplo.haplotype[colonne]);
	}
	printf("\n");
	/*(test)printf("\tdoublon=%d\n",haplo.doublon);*/
}

/* (test)Fonction d'affichage d'un genotype */
void afficher_genotypes(TypeGeno geno){
	int colonne;
	printf("/ind %d geno\t   ",geno.id);
	for (colonne = 0 ; colonne < TAILLE_GENO ; colonne++)
	{
		printf("%d",geno.genotype[colonne]);
	}
	printf("\n");
}

/*
 * Verifie si les 2 haplotypes selectionnes aleatoirement
 * sont compatibles avec les differentes regles imposees :
 *  '-> Aucun doublon
 *  '-> Il n'y a pas plus de nb_loci differents entre les 2 haplotypes. 
 */
TypeBool verification_nombre_loci(TypeHaplo haplo1, 
									TypeHaplo haplo2){
	TypeBool verifNbLoci = VRAI;
	int count = 0; /* compte le nombre de differences observees */
	int i = 0;
	while ((verifNbLoci == VRAI) && (i < TAILLE_GENO))
	{
		if(haplo1.haplotype[i] != haplo2.haplotype[i])
		{
			count ++;
		}
		/* si le comptage est superieur au nb de loci autorise,
           les haplotypes ne sont pas selectionnes */
		if(count > NB_LOCI)
		{ 
			verifNbLoci = FAUX;
		}
		i++;
	}
	
	#if 0
	afficher_haplotypes(haplo1);
	afficher_haplotypes(haplo2);
	printf("NB LOCI : %d\n",NB_LOCI);
	printf("*** %d ***\n",verifNbLoci);
	#endif
	return verifNbLoci;
}

/*
 * Fait la selection des 2 haplotypes repondant aux criteres, 
 * initialise ensuite haplo1 et 2 du genotype
 * et enfin initialise le genotype a partir des 2 haplotypes.
 */
void initialiser_genotypes(TypeGeno* adrGeno, TypeHaplo haplo[], int nbHaploNonRedondant)
{
	int h1, h2; /* nb tire aleatoirement parmi 0 et nb_haplo-1 */
	TypeBool nbLociDiffCorrect;
	h1 = tirage_au_sort(nbHaploNonRedondant);
	h2 = tirage_au_sort(nbHaploNonRedondant);
	nbLociDiffCorrect = verification_nombre_loci(haplo[h1],haplo[h2]);
	while (nbLociDiffCorrect == FAUX)
	{
		h2 = tirage_au_sort(nbHaploNonRedondant);
		nbLociDiffCorrect = verification_nombre_loci(haplo[h1],haplo[h2]);
	}
	adrGeno->haplo1.haplotype = alloue_memoire();
	adrGeno->haplo1.haplotype = haplo[h1].haplotype;
	adrGeno->haplo2.haplotype = alloue_memoire();
	adrGeno->haplo2.haplotype = haplo[h2].haplotype;
	creation_genotypes(adrGeno);
	creation_fichier_geno_haplo(adrGeno, adrGeno->haplo1, adrGeno->haplo2);
	creation_fichier_geno(adrGeno);
    free(adrGeno->genotype);
}

/* Recherche des haplotypes redondant et les marques */
void recherche_haplotype_doublon(TypeHaplo haplo[])
{
	int i;
	int j;
	TypeBool doublon;
	for (i=0 ; i < NB_HAPLO ; i++)
	{
		for (j=i+1 ; j < NB_HAPLO ; j++)
   	 	{
			doublon = 0;
            doublon = verification_presence_doublon(haplo[i],haplo[j]);
			if (doublon==VRAI)
            {
                haplo[j].id = haplo[i].id;
				haplo[j].doublon = VRAI;
			}
    	}
	}
}

/* Compte le nombre d'haplotypes non redondant */
int compte_nombre_doublon(TypeHaplo haplo[])
{
	int i;
	int compteur = 0;
	for(i=0 ; i<NB_HAPLO ; i++)
	{
		if(haplo[i].doublon == 0)
		{

			compteur ++;
		}
	}
	return compteur;
}

/* Repertorie les haplotypes non redondant */
TypeHaplo* lister_haplo_non_redondant(int compteur, TypeHaplo haplo[])
{
	int i;
	int compte = 0;
	TypeHaplo* haploNonRedondant = NULL;
	haploNonRedondant = malloc(compteur * sizeof(TypeHaplo));
	for(i=0 ; i<NB_HAPLO ; i++)
	{
		if(haplo[i].doublon == FAUX)
		{
			haploNonRedondant[compte] = haplo[i];
			compte ++;
		}
	}
	return haploNonRedondant;
}
