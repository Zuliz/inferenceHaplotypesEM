/* inference_haplotype_EM.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/inference_haplotype_EM.h"
#include "../inc/liste_doublement_chainee.h"

/* fonctions privees ============================================================================ */

 /* Recherche et retourne la frequence associee a l'id de l'haplotype voulu */
static double recherche_frequence_precedente(int id, double** tableauFrequence)
{
    int i = 0;
    bool_t trouve = FALSE;
    double frequencePrecedente;
    while (trouve == FALSE)
    {
        if(id == tableauFrequence[i][0])
        {
            trouve = TRUE;
            frequencePrecedente = tableauFrequence[i][2];
        }
        i++;
    }
    return frequencePrecedente;
}

/* Calcul des probabilites initiales pour chaque genotypes */
static void calculer_proba_geno(int nbGeno, TypeGeno** geno, double** tabFreq)
{
    int i;
    double freqH1, freqH2, proba, ppart;
    TypePaireHaplo* listePaireHaplo;
    for(i=0 ; i<nbGeno ; i++)
    {
        proba = 0.0;
        listePaireHaplo = geno[i]->first;
        while(listePaireHaplo != NULL)
        {
            freqH1 = recherche_frequence_precedente(listePaireHaplo->idHaplo1,tabFreq);
            freqH2 = recherche_frequence_precedente(listePaireHaplo->idHaplo2,tabFreq);
            #if 0
            printf("____________________ G%d ______________________\n",geno[i]->id);
            printf("freqH1 : %f - freqH2 : %f \n ",freqH1,freqH2);
            #endif
            if(listePaireHaplo->idHaplo1 == listePaireHaplo->idHaplo2)
            {
                ppart = freqH1 * freqH1;
            }
            else
            {
                ppart = 2 * freqH1 * freqH2;
            }
            proba = proba + ppart;
            #if 0
            printf("Proba : %f\n",proba);
            #endif
            listePaireHaplo = listePaireHaplo->next;
        }
        geno[i]->probaPrec = proba;
        #if 0
        printf("G%d - proba=%.10f\n",geno[i]->id,proba);
        #endif
    }
}




/* Recherche d'un genotype complet a partir d'un id */
static TypeGeno* recherche_genotype(int id, TypeGeno** geno)
{
    int i = 0;
    bool_t trouve = FALSE;
    TypeGeno* genoTrouve;
    while (trouve == FALSE)
    {
        if(id == geno[i]->id)
        {
            genoTrouve = geno[i];
            trouve = TRUE;
        }
        i++;
    }
    return genoTrouve;
}

/* Mise a jour du tableau de frequences d'haplotypes */
static void modification_frequence_precedente(int id, double freq, double** tabFreq)
{
    int i=0;
    bool_t trouve = FALSE;
    while (trouve == FALSE)
    {
        if(id == tabFreq[i][0])
        {
            tabFreq[i][2] = freq; /* Rempli la colonne des frequences precedentes */
            trouve = TRUE;
        }
        i++;
    }
}

/* Calcul du log de la vraisemblance et mise à jour des probas pour chaque genotype */
static double estimation_esperance(int nbGeno, TypeGeno** geno, double** tabFreq)
{
    int i;
    double freqH1, freqH2, proba, ppart;
    double logLikelihood = 0;
    TypePaireHaplo* listePaireHaplo;
    for(i=0 ; i<nbGeno ; i++)
    {
        proba = 0;
        listePaireHaplo = geno[i]->first;
        while(listePaireHaplo != NULL)
        {
            freqH1 = recherche_frequence_precedente(listePaireHaplo->idHaplo1,tabFreq);
            freqH2 = recherche_frequence_precedente(listePaireHaplo->idHaplo2,tabFreq);
            if(listePaireHaplo->idHaplo1 == listePaireHaplo->idHaplo2)
            {
                ppart = freqH1 * freqH1;
            }
            else
            {
                ppart = 2 * freqH1 * freqH2;
            }
            proba = proba + ppart;
            listePaireHaplo = listePaireHaplo->next;
        }
        geno[i]->probaPrec = proba;
        logLikelihood = logLikelihood + geno[i]->nbIdentique * log10(proba);
        #if 0
        printf("proba = %f ; Log10 = %f ; NGeno = %d ; vraisemblance = %f\n",proba, 
            log10(proba),geno[i]->nbIdentique,logLikelihood);
        #endif
    }
    return logLikelihood;
}







/* Calcul des frequences d'haplotypes */
static void maximisation(int nbHaplo, TypeHaplo** haplo, TypeGeno** geno, double** tabFreq)
{
    int i;
    double freqPrec1, freqPrec2, freq, contribution;
    TypeGenoExplique* listeGenoExp;
    TypeGeno* genoCourant;
    for(i=0 ; i<nbHaplo ; i++)
    {
        #if 0
        printf("========= H%d ===============\n",haplo[i]->id);
        #endif
        freqPrec1 = recherche_frequence_precedente(haplo[i]->id,tabFreq);
        #if 0
        printf("freqPrec1 : %f\n",freqPrec1);
        #endif
        freq = 0;
        /* Parcours de chaque liste chainee de genotypes */
        listeGenoExp = haplo[i]->first;
        while(listeGenoExp != NULL)
        {
            /* Recuperation des informations du genotype observe */
            genoCourant = recherche_genotype(listeGenoExp->id, geno);
            #if 0
            printf("id : %d - id2 : %d ",haplo[i]->id,listeGenoExp->idHaploCompl);
            #endif
            /* Homozygotie */
            if(haplo[i]->id == listeGenoExp->idHaploCompl)
            {
                contribution = (2.0*((freqPrec1 * freqPrec1)/genoCourant->probaPrec)*
                    (genoCourant->nbIdentique/NB_INDIV))*100;
            }
            else /* Heterozygotie */
            {
                /* printf("id pour freqH2 : %d\n",listeGenoExp->idHaploCompl); */
                freqPrec2 = recherche_frequence_precedente(listeGenoExp->idHaploCompl,tabFreq);
                contribution = (2.0*((freqPrec1 * freqPrec2)/genoCourant->probaPrec)*
                    (genoCourant->nbIdentique/NB_INDIV))*100;
                #if 0
                printf("Ngeno = %d - Ngeno/Nindiv = %f - ",genoCourant->nbIdentique,
                    genoCourant->nbIdentique/NB_INDIV);
                printf("freqH2 : %.5f - ",freqPrec2);
                printf("Contribution : %.5f\n",contribution);
                #endif
            }
            freq = freq + contribution;
            listeGenoExp = listeGenoExp->next;
        }
        freq = freq/2.0;
        #if 0
        printf("freq = %.10f\n",freq);
        #endif
        modification_frequence_precedente(haplo[i]->id,freq,tabFreq);
    }
}

/* Mise à jour des probabilites de genotypes */
static void mise_a_jour_proba(int nbGeno, TypeGeno** tabGeno)
{
    int i;
    for(i=0 ; i < nbGeno ; i++)
    {
        tabGeno[i]->proba = tabGeno[i]->probaPrec;
    }
}

/* Mise à jour des frequences d'haplotypes */
static void mise_a_jour_freq(int nbHaplo, double** tabFreq)
{
    int i;
    for(i=0 ; i < nbHaplo ; i++)
    {
        tabFreq[i][1] = tabFreq[i][2];
    }
}



/* fonctions publiques ========================================================================== */

/* Inference d'haplotypes EM avec la "methode" du seuil => pour un faible nb d'individus (< 100) */
void inference_haplotype_em_seuil(double seuil, int nbGeno, int nbHaplo, int nbEtapeMax,
                            double** tabFreq, TypeGeno** tabGeno, TypeHaplo** tabHaplo)
{
    /* Variables locales */
    bool_t convergence = FALSE;
    int nbEtape = 0;
    double vraisemblance;
    double vraisemblancePrec = 0.0;

    /* Debut */
    calculer_proba_geno(nbGeno, tabGeno, tabFreq);
    while((convergence == FALSE) && (nbEtape < nbEtapeMax))
    {
        nbEtape++;
        #if 0
        printf("*********** ITERATION %d *****************\n",nbEtape);
        #endif
        maximisation(nbHaplo, tabHaplo, tabGeno, tabFreq);
        vraisemblance = estimation_esperance(nbGeno, tabGeno, tabFreq);
        #if 0
        printf("V = %f\n",vraisemblance);
        #endif
        convergence = (fabs(vraisemblance-vraisemblancePrec)/vraisemblancePrec) <= seuil;
        #if 0
        printf("etape : %d => V : %f - VPrec: %f - ",nbEtape,vraisemblance,vraisemblancePrec);
        printf("calcul convergence : %f\n",fabs(vraisemblance-vraisemblancePrec)/vraisemblancePrec);
        #endif

        if(convergence == FALSE)
        {
            vraisemblancePrec = vraisemblance;
            mise_a_jour_proba(nbGeno, tabGeno);
            mise_a_jour_freq(nbHaplo, tabFreq);
        }
        else
        {
            printf("Nb iterations : %d\n",nbEtape);
        }
    }
    /*FIN*/
}
































/* 
 * Inference d'haplotypes EM avec la "methode" de la convergence 
 * => pour un nb d'individus > 100 : une valeur de convergence precedente et actuelle
 *  sont comparees, si ces valeurs sont identiques, la convergence est atteinte 
 */
void inference_haplotype_em(int nbGeno, int nbHaplo, int nbEtapeMax,
                            double** tabFreq, TypeGeno** tabGeno, TypeHaplo** tabHaplo)
{
    /* Variables locales */
    bool_t convergence = FALSE;
    int nbEtape = 0;
    double valConvergence=0.0;
    double valConvergencePrec=0.0;
    double vraisemblance;
    double vraisemblancePrec = 0.0;

    /* Debut */
    calculer_proba_geno(nbGeno, tabGeno, tabFreq);
    while((convergence == FALSE) && (nbEtape < nbEtapeMax))
    {
        nbEtape++;
        #if 0
        printf("*********** ITERATION %d *****************\n",nbEtape);
        #endif
        maximisation(nbHaplo, tabHaplo, tabGeno, tabFreq);
        vraisemblance = estimation_esperance(nbGeno, tabGeno, tabFreq);
        #if 0
        printf("V = %f\n",vraisemblance);
        #endif
        valConvergencePrec = (fabs(vraisemblance-vraisemblancePrec)/vraisemblancePrec);
        convergence = valConvergencePrec == valConvergence;
        #if 0
        printf("etape : %d => V : %f - VPrec: %f - ",nbEtape,vraisemblance,vraisemblancePrec);
        printf("calcul convergence : %f\n",fabs(vraisemblance-vraisemblancePrec)/vraisemblancePrec);
        #endif

        if(convergence == FALSE)
        {
            valConvergence = valConvergencePrec;
            vraisemblancePrec = vraisemblance;
            mise_a_jour_proba(nbGeno, tabGeno);
            mise_a_jour_freq(nbHaplo, tabFreq);
        }
        else
        {
            printf("Nb iterations : %d\n",nbEtape);
        }
    }
    /*FIN*/
}


























/* Recherche la paire d'haplotype la plus vraisemblable */
int* recherche_meilleure_paire_haplo(TypeGeno* liste, double** tabFreqHaplo, int taille)
{
    float freq1, freq2, freqTotal = 0.00;   
    int i, j, k = 0;            /* Compteurs */
    static int idsHaplo[2];     /* Ids de la paire d'haplotypes ayant la meilleure frequence */
    int* tabIds = NULL;         /* Recoit la paire d'id d'haplotypes a tester */
    static int* ids = idsHaplo; /* Pointeur sur la paire d'ids qui sera retourné */
    
    /* Parcours de chaque paire d'haplo */                      
    for (i=1; i <= liste->nbHaplo; i++)
    {
        tabIds = id_pos(liste, i);
        if (tabIds[0] == -1) /* Il n'y plus de paire d'haplotypes dans la liste chainee */
            break;
            
        /* Parcours du tableau des frequences pour comparer les valeurs de frequences */
        for (j=0; j < taille; j++)
        {   
            if (tabIds[0] == (int)tabFreqHaplo[j][0])
            {
                freq1 = tabFreqHaplo[j][1];
            }
            for (k=0; k < taille; k++)
            {
                 if (tabIds[1] == (int)tabFreqHaplo[k][0])
                 {
                    freq2 = tabFreqHaplo[k][1];
                 }
            }
            if ( freqTotal < freq1 + freq2)
            {
                idsHaplo[0] = tabIds[0];
                idsHaplo[1] = tabIds[1];
                freqTotal = freq1 + freq2;
            }
        }
    }

    return ids; /* Retour des ids de la paire d'haplotypes ayant la plus grande frequence */
}








/* Fonction permettant l'ecriture de la paire la plus vraisemblable pour un genotype donne */
void ecriture_paire_haplotype(TypeGeno** tabGenoNR, double** tabFreqHaplo, int nbHaploNonRedondant, 
                              TypeHaplo** tabHaploNR, FILE* resultats_haplo, int i, int k)
{
    /* Variables locales */
    int j, c = 0;
    int* idsPaireHaplo = NULL;  /* tableau contenant les ids de la paire d'haplos 
                                 * la plus vraisemblable */
    
    /* Ecriture du premier haplo de la paire */
    fprintf(resultats_haplo, "/ind %d haplo1 ", k);
    idsPaireHaplo = recherche_meilleure_paire_haplo(tabGenoNR[i], tabFreqHaplo, nbHaploNonRedondant);
    /* Recherche du 1er haplotype de la paire dans le tableau d'haplotypes non redondant */
    for (c=0; c < nbHaploNonRedondant; c++)
    {
        if (idsPaireHaplo[0] == tabHaploNR[c]->id)
        {
            for (j=0; j < TAILLE_GENO; j++)
            {
                fprintf(resultats_haplo, "%d", tabHaploNR[c]->haplotype[j]);
            }
            fprintf(resultats_haplo, "\n");
        }
    }
    
    /* Ecriture du deuxieme haplo de la paire */
    fprintf(resultats_haplo, "/ind %d haplo2 ", k);
    /* Recherche du 2eme haplotype de la paire dans le tableau d'haplotypes non redondant */
    for (c=0; c < nbHaploNonRedondant; c++)
    {
        if (idsPaireHaplo[1] == tabHaploNR[c]->id)
        {
            for (j=0; j < TAILLE_GENO; j++)
            {
                fprintf(resultats_haplo, "%d", tabHaploNR[c]->haplotype[j]);
            }
            fprintf(resultats_haplo, "\n");
        }
    } 
}