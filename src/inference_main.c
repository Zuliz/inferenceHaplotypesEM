/* inference_main.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/generation_haplotypes.h"
#include "../inc/liste_doublement_chainee.h"
#include "../inc/inference_haplotype_EM.h"

/* point d'entree =============================================================================== */

int main(int argc, char* argv[])
{  
    /* ====== Declarations ====== */

    int c, i, j, k = 0; /* compteurs */
    int id = 0;
    int idGenoMax = 0;
    int nbLoci = 0;
    int nbEtapeMax = 1000;
    int nbHaploNonRedondant = 0;
    int nbGenoNonRedondant = 0;
    int lireTaille = 0;         /* taille memoire pour la lecture de ligne */
    /*int* idsPaireHaplo = NULL;   tableau contenant les ids de la paire d'haplos la plus vrai semblable */
    char* nomFichier;           /* fichier contenant les genotypes */
    char* chaine = NULL;        /* chaine qui contiendra les lignes du fichier geno */
    char* sousChaine = NULL;    /* pointeur sur la partie code du genome */
    double** tabFreqHaplo = NULL;
    double seuil = 0.000000000000000000000000001;
    FILE* fichier = NULL;       
    FILE* fichierParam = NULL;
    FILE* resultats_freq = NULL;
    FILE* resultats_haplo = NULL;
    TypeGeno* geno = NULL;
    TypeHaplo** tabHaploNR = NULL; /* tableau d'haplotypes non redondant */
    TypeGeno** tabGenoNR = NULL;   /* tableau de genotypes non redondant */
    
     /* ========== Code ========== */
    
    if (argc != 2)
    {
        printf("Le nombre d'argument est incorrect !\n");
        printf("la commande doit etre du type :\n");
        printf("./__FILE__ fichier_genotype(chaine) \n");
        exit(1);    
    }
    
    nomFichier = argv[1];
    
    fichier = fopen(nomFichier, "r");
    if (fichier != NULL)
    {
        fichierParam = fopen(PARAM, "r");
        if (fichierParam != NULL)
        {
            fscanf(fichierParam, "%lf", &nbIndiv);
            fscanf(fichierParam, "%d", &tailleGeno);
            fscanf(fichierParam, "%d", &nbLoci);
        }
        else
        {
            fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n", PARAM);
            exit(1);
        }
        fclose(fichierParam);
        
        /* Determination de la taille d'une ligne */
        while((c != '\n') && (c != EOF))
        {
            c = fgetc(fichier);
            lireTaille++;
        }
        /* taille de la ligne + '\0' + augmentation de caracteres par nbIndiv: 10->100->1000.. */
        lireTaille = lireTaille + 1 + (int)((log10(NB_INDIV) + 1)); 
        
        fseek(fichier, 0, SEEK_SET);
        chaine = (char*)malloc(sizeof(char) * lireTaille);
        
        /* Stockage des genotypes du fichier dans un tableau de structures */
        c = 0;
        geno = malloc(sizeof(TypeGeno) * NB_INDIV);
        for(i=0; i < NB_INDIV; i++)
        {
            geno[i].genotype = (int*)malloc(sizeof(int) * TAILLE_GENO);
        }
        while((lire(chaine, lireTaille, fichier) == 0))
        {
            geno[c].id = c;
            sousChaine = strchr(chaine, '\t');
            sousChaine = sousChaine + 1;
            if(sousChaine != NULL)
            {
                for(i=0; i < TAILLE_GENO; i++)
                {
                    geno[c].genotype[i] = sousChaine[i] - '0';
                    /*printf("%d",geno[c].genotype[i]);*/
                }
                c++;
                /*printf("\n");*/
            }
        }
        free(chaine);
    }
    else
    {
        fprintf(stderr, 
                "Une erreur s'est produite durant l'ouverture du fichier des génotypes.\n");
        exit(1);
    }
    fclose(fichier);
    
    /* Generation des haplotypes possibles */
    for(i=0 ; i<NB_INDIV ; i++)
    {
        /* Recuperation de l'id pour continuer la numerotation 
         * des haplotypes du genotype suivant */
        id = initialisation_geno(&geno[i],id); 
    }

    /* Indique pour chaque genotype et haplotype si celui-ci est un doublon ou non */
    for (i=0 ; i < NB_INDIV ; i++)
    {
        j = i+1;
        while (j < NB_INDIV)
        {
            /* Verification des genotypes identiques et remplacement de l'id 
             * en cas de similarite */ 
            recherche_genotype_doublon(&geno[i], &geno[j]);
            recherche_haplotype_doublon(&geno[i], &geno[j]);
            j++;
        }
        /*afficher_haplotypes(haplo[i],TAILLE_GENO);*/
        /*** TEST 1 ***/
    }

    #if 0
    for (i=0 ; i < NB_INDIV ; i++)
    {
        printf("*********************\n");
        affichage_genotype(geno[i]);
        printf("---------------------\n");
        affichage_haplotypes(geno[i]);
    }
    #endif

    #if 0
    for (i=0 ; i < NB_INDIV ; i++)
    {
        printf("G%d (d=%d) + nbIdentique = %d : ",geno[i].id,geno[i].doublon,geno[i].nbIdentique);
        for(j=0 ; j<geno[i].nbHaplo ; j++)
        {
            printf("H%d - ",geno[i].matriceHaplo[j].id);
        }
        printf("\n");
    }
    #endif

    /* Calcul du nombre total et du nombre d'haplotypes non redondant */
    nbHaploNonRedondant = calcul_nb_haplo_non_redondant(geno);
    nbGenoNonRedondant = calcul_nb_geno_non_redondant(geno);
    #if 0
    printf("Nb de Haplo non redondant : %d\n", nbHaploNonRedondant);
    printf("Nb de Geno non redondant : %d\n", nbGenoNonRedondant);
    #endif

    /* Allocation du tableau d'haplotypes non redondant */
    tabHaploNR = (TypeHaplo**)malloc(sizeof(TypeHaplo*) * nbHaploNonRedondant);
    if (tabHaploNR == NULL)
    {
        fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
        exit(1);
    }

	/* Initialisation de chaque liste du tableau d'haplotypes */
    for (i = 0; i < nbHaploNonRedondant; i++)
    {
        tabHaploNR[i] = initialiser_liste_geno(TAILLE_GENO);
    }
    
    /* Affection des haplotypes dans leur tableau */
    c = 0;
    for (i=0 ; i < NB_INDIV ; i++)
    {
        for (j=0; j < geno[i].nbHaplo; j++)
        {
            /*printf("G%d - H%d\n",geno[i].id,geno[i].matriceHaplo[j].id);*/
            if (geno[i].matriceHaplo[j].doublon == 0)
            {
                tabHaploNR[c]->id = geno[i].matriceHaplo[j].id;
                for (k=0; k < TAILLE_GENO; k++)
                {
                    tabHaploNR[c]->haplotype[k] = geno[i].matriceHaplo[j].haplotype[k];
                }
                c++;
            }
        }
    }
    #if 0 /* Test qui regarde les haplotypes dans chaque cellule du tableau */
    c=0;
    for (i=0; i < nbHaploNonRedondant; i++)
    {
        for (k=0; k < TAILLE_GENO; k++)
        {
            printf("%d", tabHaploNR[c]->haplotype[k]);
        }
        printf("\n");
        c++;
    }
    /*** TEST 2 ***/
	#endif
    
	/* Ajout des genotypes expliques par l'haplotype dans la liste correspondante */
    for (c=0; c < nbHaploNonRedondant; c++)
    {
        idGenoMax = 0;
        for (i=0; i < NB_INDIV; i++)
        {
            /* Genotype actuel non redondant donc introduit dans le tableau d'haplotypes au niveau
             * des listes chainees */
            if(idGenoMax <= geno[i].id)
            {
                for (j=0; j < geno[i].nbHaplo; j++)
                {
                    if (tabHaploNR[c]->id == geno[i].matriceHaplo[j].id)
                    {
                        ajout_queue_geno(tabHaploNR[c], geno[i].id,geno[i].matriceHaplo[recherche_haplo_complementaire(geno[i],j)].id);
                    }
                }
                idGenoMax = geno[i].id;
            }
        }        
    }
    #if 0 /* Affiche le contenu de toutes les listes chainees du tableau d'haplotypes NR */
    printf("-----test-----\n");
    for (i=0; i < nbHaploNonRedondant; i++)
    {
        printf("H%d: ", tabHaploNR[i]->id);
        affichage_liste_geno(tabHaploNR[i]);
    }
    #endif
    
	/* Initialisation du tableau de frequences d'haplotypes */
    tabFreqHaplo = allouer_memoire_tableau_2d(nbHaploNonRedondant); 
    for (i=0 ; i<nbHaploNonRedondant ; i++)
    {
        tabFreqHaplo[i][0] = tabHaploNR[i]->id;
        tabFreqHaplo[i][1] = 1.00/nbHaploNonRedondant*100;
        tabFreqHaplo[i][2] = 1.00/nbHaploNonRedondant*100;
    }
    #if 0
    for(i=0 ; i<nbHaploNonRedondant ; i++)
    {
        printf("%f %f %f\n",tabFreqHaplo[i][0],tabFreqHaplo[i][1],tabFreqHaplo[i][2]);
    }
    #endif

    /* Allocation du tableau de genotypes */
    tabGenoNR = (TypeGeno**)malloc(sizeof(TypeGeno*) * nbGenoNonRedondant);
    if (tabGenoNR == NULL)
    {
        fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
        exit(1);
    }

    /* Initialisation de chaque liste du tableau d'haplotypes */
    for (i = 0; i < nbGenoNonRedondant; i++)
    {
        tabGenoNR[i] = initialiser_liste_haplo(TAILLE_GENO);
    }
    
    /* Affection des genotypes dans leur tableau */
    c = 0;
    for (i=0; i < NB_INDIV; i++)
    {
        if (geno[i].doublon == 0)
        {
            tabGenoNR[c]->id = geno[i].id;
            tabGenoNR[c]->nbIdentique = geno[i].nbIdentique;
            tabGenoNR[c]->nbHaplo = geno[i].nbHaplo;

            for (j=0; j < TAILLE_GENO; j++)
            {
                tabGenoNR[c]->genotype[j] = geno[i].genotype[j];
            }
            c++;
        }
    }
    #if 0 /* Regarde l'id et le genotype de chaque genotype du tableau */
    for (i=0; i < nbGenoNonRedondant; i++)
    {
        printf("ID:%d -->", tabGenoNR[i]->id);
        for (j=0; j < TAILLE_GENO; j++)
        {
            printf("%d", tabGenoNR[i]->genotype[j]);
        }
        printf("\n");
    }
    #endif
    
    /* Ajout des paires d'haplotypes a chaque liste chainee du tableau de genotypes */
    c = 0; /* Compteur permettant d'avancer dans le tableau de genotypes lorsque il y a un 
            * nouveau genotype non redondant de trouve */
    for (i=0; i < NB_INDIV; i++)
    {
        if (geno[i].doublon == 0)
        {
            /* Si le genotype ne possede pas de loci ambigus (homozygotie) */
            if (geno[i].nbHaplo == 1)
            {
                ajout_queue_paire_haplo(tabGenoNR[c], geno[i].matriceHaplo[0].id, 
                                        geno[i].matriceHaplo[0].id);
            }
            else /* Heterozygotie */

            {
                for (j=0; j < (geno[i].nbHaplo / 2); j++)
                {
                    ajout_queue_paire_haplo(tabGenoNR[c], geno[i].matriceHaplo[j].id, 
                        geno[i].matriceHaplo[recherche_haplo_complementaire(geno[i], j)].id);
                }
            }
            c++; /* Incremente lorsque le genotype observe n'est pas un doublon */

        }
    }
    
    #if 1 /* Affiche le contenu de toutes les listes chainees du tableau de genotypes NR */
    printf("-----test-----\n");
    for (i=0; i < nbGenoNonRedondant; i++)
    {
        printf("G%d (nb = %d): ", tabGenoNR[i]->id, tabGenoNR[i]->nbIdentique);
        affichage_liste_paire_haplo(tabGenoNR[i]);
    }
    #endif

    /* Inference d'haplotypes EM */
    if (NB_INDIV < 100)
    {
        inference_haplotype_em_seuil(seuil, nbGenoNonRedondant, nbHaploNonRedondant, nbEtapeMax,
                           tabFreqHaplo, tabGenoNR, tabHaploNR);
    }
    else
    {
       inference_haplotype_em(nbGenoNonRedondant, nbHaploNonRedondant, nbEtapeMax,
                           tabFreqHaplo, tabGenoNR, tabHaploNR);
    }
        
    /* Tri du tableau de frequence d'haplotype */ 
    qsort (tabFreqHaplo, nbHaploNonRedondant , sizeof*tabFreqHaplo, compare);
    
    #if 0
    for(i=0 ; i<nbHaploNonRedondant ; i++)
    {
        printf("%f %f %f\n",tabFreqHaplo[i][0],tabFreqHaplo[i][1],tabFreqHaplo[i][2]);
    }
    #endif

    /* Creation du fichier des haplotypes avec leur fréquences triés par ordre decroissant */
    resultats_freq = fopen(RESULTATS_FREQ, "w");
    if ( fichier != NULL)
    {
        for (i=0; i < nbHaploNonRedondant ;i++)
        {
            for (j=0; j < nbHaploNonRedondant; j++)
            {
                if ((tabHaploNR[j]->id == tabFreqHaplo[i][0]) && (tabFreqHaplo[i][1] > 0.01))
                {
                    for (c=0; c < TAILLE_GENO; c++)
                    {
                        fprintf(resultats_freq, "%d", tabHaploNR[i]->haplotype[c]);
                    }
                    fprintf(resultats_freq, "\t%.2f\n",tabFreqHaplo[i][1]);
                }
            }
        }
    }
    else
    {
        fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n",
                RESULTATS_FREQ);
        exit(1);
    }
    fclose(resultats_freq);
    
    /* Creation du fichier de genotypes par individu avec leur paires d'haplotypes la plus 
     * vraisemblable */ 
    resultats_haplo = fopen(RESULTATS_HAPLO, "w");
    if ( fichier != NULL)
    {
        for (i=0; i < NB_INDIV; i++)
        {
             /* Ecriture de la ligne geno */
            fprintf(resultats_haplo, "/ind %d geno   ", i);
            for (j=0; j < TAILLE_GENO; j++)
            {
                fprintf(resultats_haplo, "%d", geno[i].genotype[j]);
            }
            fprintf(resultats_haplo, "\n");
            id = 0;
            /* Regarde à quel indice dans le tableau de genotypes non redondant se trouve
             * le genotype de l'individu */
            for (c=0; c < nbGenoNonRedondant; c++)
            {
                    if (tabGenoNR[c]->id == geno[i].id)
                        id = c;
            }
            /* Ecriture de la paire d'haplo explicative */
            ecriture_paire_haplotype(tabGenoNR, tabFreqHaplo, nbHaploNonRedondant, 
                                     tabHaploNR, resultats_haplo, id, i);
        }
    }
    else
    {
        fprintf(stderr, "Une erreur s'est produite dans la création du fichier %s\n", 
                RESULTATS_HAPLO);
        exit(1);
    }
    fclose(resultats_haplo);

    /* Liberation de la memoire alloue au tableau de frequences des haplotypes */
    for (i=0; i < nbHaploNonRedondant; i++)
    {
        free(tabFreqHaplo[i]);
        tabFreqHaplo[i] = NULL;
    }
    free(tabFreqHaplo);
    
    /* Liberation de la memoire allouee au tableau de genotypes contenant les donnees du fichier 
     * de genotypes */
    for(i=0; i < NB_INDIV; i++)
    {   
        for (j=0; j < TAILLE_GENO; j++)
        {
            free(geno[i].genotype);
            geno[i].genotype = NULL;
        }
        
        for (j=0; j < geno[i].nbHaplo; j++)
        {
            free(geno[i].matriceHaplo[j].haplotype);
        }
        
        free(geno[i].matriceHaplo);
        geno[i].matriceHaplo = NULL;
    }
    free(geno);
    
    /* Liberation de la memoire de chaque element du tableau d'haplotypes puis du tableau */
    for (i = 0; i < nbHaploNonRedondant; i++)
    {
        /* Liberation de l'espace alloue pour les tableaux d'entiers representant les haplotypes */
        for (j=0; j < TAILLE_GENO; j++)
        {
            free(tabHaploNR[i]->haplotype);
            tabHaploNR[i]->haplotype = NULL;
        }
        /* Suppression de chaque liste chainee */
        suppression_liste_geno(&tabHaploNR[i]);
        
        /* Liberation de la memoire allouer pour chaque cellule du tableau d'haplotypes */
        free(tabHaploNR[i]);
        tabHaploNR[i] = NULL;
    }
    /* Liberation du pointeur sur le tableau d'haplotypes */
    free(tabHaploNR);
    
     /* Liberation de la memoire de chaque element du tableau de genotypes puis du tableau */
    for (i = 0; i < nbGenoNonRedondant; i++)
    {
        /* Liberation de l'espace alloue pour les tableaux d'entiers representant les genotypes */
        for (j=0; j < TAILLE_GENO; j++)
        {
            free(tabGenoNR[i]->genotype);
            tabGenoNR[i]->genotype = NULL;
        }
        /* Suppression de chaque liste chainee */
        suppression_liste_paire_haplo(&tabGenoNR[i]);
        
        /* Liberation de la memoire allouer pour chaque cellule du tableau de genotypes */
        free(tabGenoNR[i]);
        tabGenoNR[i] = NULL;
    }
    /* Liberation du pointeur sur le tableau de genotypes */
    free(tabGenoNR);

    return 0;
}
