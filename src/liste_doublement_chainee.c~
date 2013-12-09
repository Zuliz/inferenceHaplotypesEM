/* liste_doublement_chainee.c */

/*
 * Auteurs: Julie Pelletier
 *          Nicolas Belouin
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/liste_doublement_chainee.h"

#define NB_HAPLO 5
#define TAILLE_GENO 10

/* fonctions publiques ========================================================================== */

/* Fonction d'initilisation d'une liste chainee de genotype */
TypeHaplo* initialiser_liste_geno(int tailleGeno)
{
    /* Allocation de l'haplotype */
    TypeHaplo* liste = malloc(sizeof(TypeHaplo));
    if (liste == NULL)
    {
        fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
        exit(1);
    }
    
    liste->haplotype = (int*)malloc(sizeof(int) * tailleGeno);
    liste->id = 0;
    liste->freq = 0.00;
    liste->taille = 0;
    liste->teteGeno = NULL;
    liste->first = NULL;
    liste->last = NULL;
    
    return liste;
}

/* Focntion d'initilisation d'une liste chainee de paires d'haplotypes */
TypeGeno* initialiser_liste_haplo(int tailleGeno)
{
    /* Allocation de l'haplotype */
    TypeGeno* liste = malloc(sizeof(TypeGeno));
    if (liste == NULL)
    {
        fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
        exit(1);
    }
    
    liste->id = 0;
    liste->taille = 0;
    liste->nbHaplo = 0;
    liste->nbIdentique = 1;
    liste->probaPrec = 0;
    liste->proba = 0;
    liste->genotype = (int*)malloc(sizeof(int) * tailleGeno);
    liste->tetePaireHaplo = NULL;
    liste->first = NULL;
    liste->last = NULL;
    
    return liste;
}

/* Fonction d'ajout d'un element en debut de liste */
void ajout_tete(TypeHaplo* liste, int id)
{
    if (liste != NULL)
    {
        /* Creation du nouvel element de la liste */
        TypeGenoExplique* nouvelleCase = malloc(sizeof(TypeGenoExplique));
        if (nouvelleCase == NULL)
        {
            fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
            exit(1);
        }
        nouvelleCase->id =id;
        nouvelleCase->prev = NULL;
        
        /* Cas ou la liste est vide */
        if (liste->first == NULL)
        {
            liste->first = nouvelleCase;
            liste->last = nouvelleCase;
            nouvelleCase->next = NULL;
        }
        else
        {
            nouvelleCase->next = liste->first;
            liste->first->prev = nouvelleCase;
            liste->first = nouvelleCase;
        }
        liste->taille++; /* Incrementation de la taille de la liste */
    }  
}

/* Fonction d'ajout d'un genotype en fin de liste */
void ajout_queue_geno(TypeHaplo* liste, int id, int idHaplo)
{
    if (liste != NULL)
    {
        /* Creation du nouvel element de la liste */
        TypeGenoExplique* nouvelleCase = malloc(sizeof(TypeGenoExplique));
        if (nouvelleCase == NULL)
        {
            fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
            exit(1);
        }
        nouvelleCase->id =id;
        nouvelleCase->idHaploCompl = idHaplo;
        nouvelleCase->next = NULL;
        
        /* Cas ou la liste est vide */
        if (liste->last == NULL)
        {
            nouvelleCase->prev = NULL;
            liste->first = nouvelleCase;
            liste->last = nouvelleCase;
        }
        else
        {
            nouvelleCase->prev = liste->last;
            liste->last->next = nouvelleCase;
            liste->last = nouvelleCase;
        }
        liste->taille++; /* Incrementation de la taille de la liste */
    }
}

/* Fonction d'ajout d'une paire d'haplotypes en fin de liste */
void ajout_queue_paire_haplo(TypeGeno* liste, int id1, int id2)
{
    if (liste != NULL)
    {
        /* Creation du nouvel element de la liste */
        TypePaireHaplo* nouvelleCase = malloc(sizeof(TypeGenoExplique));
        if (nouvelleCase == NULL)
        {
            fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
            exit(1);
        }
        nouvelleCase->idHaplo1 = id1;
        nouvelleCase->idHaplo2 = id2;
        nouvelleCase->next = NULL;
        
        /* Cas ou la liste est vide */
        if (liste->last == NULL)
        {
            nouvelleCase->prev = NULL;
            liste->first = nouvelleCase;
            liste->last = nouvelleCase;
        }
        else
        {
            nouvelleCase->prev = liste->last;
            liste->last->next = nouvelleCase;
            liste->last = nouvelleCase;
        }
        liste->taille++; /* Incrementation de la taille de la liste */
    }
}

#if 0
/* Fonction permettant l'ajout d'un element a la position souhaitee */
void ajout_pos(TypeHaplo* liste, int id, int position)
{
    int i = 1;          /* Compteur */
    TypeGenoExplique* ptr = NULL;   /* Pointeur de parcours de liste */
    
    
    if (liste != NULL)
    {
        ptr = liste->first;
        
        while ((i <= position) && (ptr != NULL))
        {
            if (position == i)
            {
                /* Cas ou position est en fin de liste */
                if (ptr->next == NULL)
                {
                    ajout_queue_geno(liste, id);
                }
                /* Cas ou position est en tete de liste */
                else if (ptr->prev == NULL)
                {
                    ajout_tete(liste, id);
                }
                else
                {
                    /* Creation du nouvel element de la liste */
                    TypeGenoExplique* nouvelleCase = malloc(sizeof(TypeGenoExplique));
                    if (nouvelleCase == NULL)
                    {
                        fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
                        exit(1);
                    }
                    nouvelleCase->id =id;
                    ptr->next->prev = nouvelleCase;
                    ptr->prev->next = nouvelleCase;
                    nouvelleCase->prev = ptr->prev;
                    nouvelleCase->next = ptr;
                    liste->taille++;
                }
            }
            else
            {
                ptr = ptr->next;
            }
            i++;
        }
    }
}
#endif

/* Fonction supprimant un element en fonction de sa position */
void sup_pos(TypeHaplo* liste, int position)
{
    int i = 1;                   /* Compteur */
    TypeGenoExplique* ptr = liste->first; /* Pointeur de parcours de liste */
    
    if (liste != NULL)
    {
        while ((ptr != NULL) && (i <= position))
        {
            if (position == i)
            {
                /* Cas ou l'on se trouve en fin de liste */
                if (ptr->next == NULL)
                {
                    liste->last = ptr->prev;
                    liste->last->next = NULL;
                }
                /* Cas ou l'on se trouve en début de liste */
                else if (ptr->prev == NULL)
                {
                    liste->first = ptr->next;
                    liste->first->prev = NULL;
                }
                else
                {
                    ptr->next->prev = ptr->prev;
                    ptr->prev->next = ptr->next;
                }
                free(ptr);
                liste->taille--;
            }
            else
            {
                ptr = ptr->next;
            }
            i++;
        }
    }
}

/* Fonction supprimant la case contenant l'id passe en parametre */
void sup_case_id(TypeHaplo* liste, int id)
{
    bool_t trouve = 0;
    TypeGenoExplique* ptr = liste->first; /* Pointeur de parcours de liste */
    
    if (liste != NULL)
    {
        while ((ptr != NULL) && (!trouve))
        {
            if (ptr->id == id)
            {
                /* Cas ou l'on se trouve en fin de liste */
                if (ptr->next == NULL)
                {
                    liste->last = ptr->prev;
                    liste->last->next = NULL;
                }
                /* Cas ou l'on se trouve ne debut de liste */
                else if (ptr->prev == NULL)
                {
                    liste->first = ptr->next;
                    liste->first->prev = NULL;
                }
                else
                {
                    ptr->next->prev = ptr->prev;
                    ptr->prev->next = ptr->next;
                }
                free(ptr);
                liste->taille--;
                trouve = 1;
            }
            else
            {
                ptr=ptr->next;
            }
        }
        
    }
}

/* Fonction supprimant tous les elements ayant le meme id */
void sup_ids(TypeHaplo* liste, int id)
{
    TypeGenoExplique* ptr = liste->first; /* Pointeur de parcours de liste */
    
    if (liste != NULL)
    {
        while (ptr != NULL)
        {
            if (ptr->id == id)
            {
                /* Cas ou l'on se trouve en fin de liste */
                if (ptr->next == NULL)
                {
                    liste->last = ptr->prev;
                    liste->last->next = NULL;
                }
                /* Cas ou l'on se trouve ne debut de liste */
                else if (ptr->prev == NULL)
                {
                    liste->first = ptr->next;
                    liste->first->prev = NULL;
                }
                else
                {
                    ptr->next->prev = ptr->prev;
                    ptr->prev->next = ptr->next;
                }
                free(ptr);
                liste->taille--;
            }
            else
            {
                ptr=ptr->next;
            }
        }
        
    }
}

/* Fonction supprimant entierement la liste chainee de genotypes */
void suppression_liste_geno(TypeHaplo** ptrListe)
{
    TypeGenoExplique* ptr = (*ptrListe)->first; /* pointeur de parcours de liste */
    TypeGenoExplique* ptrDel = NULL;            /* pointeur de suppression d'element */
    
    if (*ptrListe != NULL)
    {
        while (ptr != NULL)
        {
            ptrDel = ptr;
            ptr = ptr->next;
            free(ptrDel);
        }
        free(*ptrListe);
        *ptrListe = NULL;
    }
}

/* Fonction supprimant entierement la liste chainee de paires d'haplotypes */
void suppression_liste_paire_haplo(TypeGeno** ptrListe)
{
    TypePaireHaplo* ptr = (*ptrListe)->first; /* pointeur de parcours de liste */
    TypePaireHaplo* ptrDel = NULL;            /* pointeur de suppression d'element */
    
    if (*ptrListe != NULL)
    {
        while (ptr != NULL)
        {
            ptrDel = ptr;
            ptr = ptr->next;
            free(ptrDel);
        }
        free(*ptrListe);
        *ptrListe = NULL;
    }
}

/* Fonction de recherche en fonction de l'id dans une liste chainee de genotype */
TypeGenoExplique* recherche_id_geno(TypeHaplo* liste, int id)
{
    TypeGenoExplique* ptr = liste->first; /* Pointeur de parcours de liste */
    
    if (liste != NULL)
    {
        while (ptr != NULL)
        {
            if (ptr->id == id)
            {
                return ptr;
            }
            else
            {
                ptr=ptr->next;
            }
        }
        
    }
    
    return ptr;
}

/* Fonction de recherche en fonction de l'id dans une liste chainee de paire d'haplotypes */
TypePaireHaplo* recherche_id_paire_haplo(TypeGeno* liste, int id)
{
    TypePaireHaplo* ptr = liste->first; /* Pointeur de parcours de liste */
    
    if (liste != NULL)
    {
        while (ptr != NULL)
        {
            if (ptr->idHaplo1 == id)
            {
                return ptr;
            }
            else
            {
                ptr=ptr->next;
            }
        }
        
    }
    
    return ptr;
}

/* Fonction modifiant les informations d'un element de la liste */
void modif_liste(TypeHaplo* liste, int id, int new_id)
{
    TypeGenoExplique* ptr = recherche_id_geno(liste, id);
    if (ptr == NULL)
    {
        fprintf(stderr, "La modification de liste n'a pu se faire.\n");
        return;
    }
    else
    {
        ptr->id = new_id;
    }
}

/* Fonction affichant chaque genotype de la liste */
void affichage_liste_geno(TypeHaplo* liste)
{
    TypeGenoExplique* ptr = liste->first; /* Poitneur de parcours de liste */
    
    if (liste != NULL)
    {
        while (ptr != NULL)
        {
            printf("%d/%d > ", ptr->id,ptr->idHaploCompl);
            ptr = ptr->next;
        }
        printf("NULL\n");
    }
    else
    {
        printf("La liste est vide.\n");
    }
}

/* Fonction affichant chaque paire d'haplotypes de la liste */
void affichage_liste_paire_haplo(TypeGeno* liste)
{
    TypePaireHaplo* ptr = liste->first; /* Poitneur de parcours de liste */
    
    if (liste != NULL)
    {
        while (ptr != NULL)
        {
            printf("%d/%d > ", ptr->idHaplo1, ptr->idHaplo2);
            ptr = ptr->next;
        }
        printf("NULL\n");
    }
    else
    {
        printf("La liste est vide.\n");
    }
}

/* Fonction retournant la taille d'une liste */
int taille_liste(TypeHaplo* liste)
{
    int taille = 0;
    
    if (liste != NULL)
    {
        taille = liste->taille;
    }
    
    return taille;
}

/* Fonction donnant l'idHaplo1 d'une paire d'haplotypes à telle position */
int* id_pos(TypeGeno* liste, int position)
{
    int i = 1; /* Compteur */
    int tabIds[2];
    int* ids = tabIds;
    
    TypePaireHaplo* ptr = liste->first; /* Pointeur de parcours de liste */
    
    if (liste != NULL)
    {
        while ((ptr != NULL) && (i <= position))
        {
            if (position == i)
            {
                /*return ptr->idHaplo1;*/
                tabIds[0] = ptr->idHaplo1;
                tabIds[1] = ptr->idHaplo2;
                return ids;
            }
            else
            {
                ptr = ptr->next;
            }
            i++;
        }
    }
    ids[0] = -1;
    ids[1] = -1;
    return ids;
}


/* point d'entree =============================================================================== */
#if 0
int main()
{
    /* ====== Declarations ====== */

    int i = 0;
    
    /* Allocation du tableau d'haplotypes */
    TypeHaplo** tabListe = (TypeHaplo**)malloc(sizeof(TypeHaplo*) * 100);
    if (tabListe == NULL)
    {
        fprintf(stderr, "Un problème d'allocation mémoire est survenu.\n");
        exit(1);
    }
    
    /* Initialisation de chaque liste du tableau d'haplotypes */
    for (i = 0; i < 100; i++)
    {
        tabListe[i] = initialiser_liste(10);
    }

    /* ========== Code ========== */
    
    ajout_queue_geno(tabListe[0], 1);
    affichage_liste(tabListe[0]);
    printf("%f\n", tabListe[0]->freq);
    
    /* Liberation de la mémoire de chaque element du tableau puis du tableau */
    for (i = 0; i < 100; i++)
    {
        suppression_liste_geno(&tabHaploNR[i]);
        
        free(tabListe[i]);
        tabListe[i] = NULL;
    }
    free(tabListe);
    
    return 0;
}
#endif
