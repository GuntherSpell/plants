#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <vector>
#include <array>
#include <fstream>

#include "patch.h"


/**
 * @file
 */

 /**
  * @brief
  * Énumération qui permet de choisir
  * le type de distribution que l'on
  * veut pour les mutations.
  */
typedef enum _mut_
{
    gaussian = 0,
    uniform = 1,
} distrMut;

/**
 * @brief
 * Contient les caractéristiques d'un monde.
 * Cette classe contient également les méthodes
 * permettant l'exécution du modèle.
 */

class World
{
public:

    /**
     * @brief
     * Constructeur d'un monde
     *
     * @param idWorld   Permet de générer des rapports différents pour chaque monde.
     * @param Kmin      La capacité d'accueil minimale.
     * @param Kmax      La capacité d'accueil maximale.
     * @param sigmaK    Le degré de varition de K dans l'espace.
     * @param Pmin      La probabilité maximale qu'un patch soit pollinisé.
     * @param Pmin      La probabilité minimale qu'un patch soit pollinisé.
     * @param sigmaP    Le degré de varition de P dans l'espace.
     */
    World(int idWorld, int NPatch, double delta, double c, bool relationshipIsManaged, int typeMut, double mu, double sigmaZ, int Kmin, int Kmax, int sigmaK,
          double Pmin, double Pmax, double sigmaP, double sInit, double dInit, int NGen, int genReport);

    /**
     * @brief
     * Méthode qui lance la simulation
     *
     * @param idWorld   Permet d'identifier le monde sur l'écran de progression
     */
    void run(int idWorld);

private:

    int NPatch; /**< @brief Nombre de patchs du monde */

    std::vector<Patch> patches; /**< @brief Vecteur qui contient tous les patchs du monde */

    double delta; /**< @brief La dépression de consanguinité */
    double c; /**< @brief Le coût de dispersion */

    bool relationshipIsManaged; /**< @brief Indique si on doit gérer l'apparentement */
    int Ktot; /**< @brief La somme de la capacité d'accueil de tous les patchs. Utile uniquement si on gère l'apparentement */

    distrMut typeMut; /**< @brief La distribution de l'ampleur mutation */
    double mu; /**< @brief La probabilité de mutation */
    double sigmaZ; /**< @brief L'ampleur de la mutation */

    int NGen; /**< @brief Le nombre de générations à créer */
    int genReport; /**< @brief Le nombre de générations entre chaque rapport .txt */

    std::mt19937_64 generator; /**< @brief Générateur de nombre aléatoire */


    /**
     * @brief Deux vecteurs qui contiennent temporairement la nouvelle génération d'un patch
     *
     * Pour plus de détails, voir la méthode createNextGen
     */
    std::array<std::vector<Individual>,2> juveniles;

    /**
     * @brief
     * Vecteur de vecteurs (demi-matrice) qui contient
     * tous les apparentements entre tous les individus.
     *
     * La diagonale n'est pas remplie car, en cas d'autof, on
     * utilisera le coefficient de consanguinité de la mère.
     * On remplit la moitié supérieure de la matrice, c'est-à-dire
     * que la première ligne est pleine et la dernière est vide.
     */
    std::array<std::vector<std::vector<double>>> relationship;

    /**
     * @brief
     * vecteur qui contient tous les pères choisis pour pouvoir récréer
     * la matrice d'apparentement. NB: en cas d'autof, père = mère.
     */
    std::vector<int> fathers;

    /**
     * @brief
     * vecteur qui contient toutes les mères choisies pour pouvoir récréer
     * la matrice d'apparentement.
     */
    std::vector<int> mothers;

    /**
     * @brief
     * Méthode qui permet de retourner une valeur pour un patch selon sa position.
     *
     * @param minVal    Valeur maximale
     * @param maxVal    Valeur minimale
     * @param sigma     Degré de variation
     * @param i         La position du patch
     *
     * @return          La valeur pour le patch donné.
     */
     double distr(double minVal, double maxVal, double sigma, int posPatch);

    /**
     * @brief
     * Méthode qui crée la nouvelle génération d'un patch en argument
     *
     * Afin d'éviter de garder en mémoire des données inutiles
     * c'est seulement à cette étape que sont calculées les pressions
     * en propagules. De plus, on ne calcule que les pressions utiles
     * pour la création de la nouvelle génération de ce patch.
     * La nouvelle génération est stockée temporairement puis elle
     * remplace l'ancienne dès que cette dernière n'est plus utile
     * pour les calculs (i.e. quand on a créé la génération à gauche
     * du patch, puisqu'on progresse de gauche à droite).
     *
     * @param idPatch   L'identifiant du patch dont on souhaite créer la nouvelle génération.
     */
    void createNextGen(int idPatch);

    /**
     * @brief
     * Méthode qui crée un nouvel individu selon la propagule choisie
     *
     * @param whr           indique dans quel vecteur temporaire il faut stocker la génération.
     * @param patchMother   identifiant du patch de la mère
     * @param mother        identifiant de la mère
     * @param autof         si la graine est issue d'autof ou non
     */
    void newInd(int whr, int patchMother, int mother, bool autof);


    /**
     * @brief
     * Méthode qui va renvoyer la position absolue d'un individu dans le monde
     * (e.g. les indivdus du patch 0 sont les premiers).
     *
     * Cette méthode ne sert que quand on gère l'apparentemment.
     */
    int getAbsolutePos (int posInPatch, int idPatch);

    /**
     * @brief
     * Méthode qui applique une mutation aléatoire
     * sur un des traits de l'individu
     *
     * @param IndToMutate   L'individu à muter
     */
    void mutation(Individual& IndToMutate);

    /**
      * @brief
      * Crée une mutation selon une loi uniforme.
      *
      * @param t        Valeur du trait à muter.
      *
      * @return         La valeur du trait après mutation.
      */
    double unifMutation (double t);

    /**
      * @brief
      * Crée une mutation selon une loi normale.
      *
      * @param t        Valeur du trait à muter.
      *
      * @return         La valeur du trait après mutation.
      */
    double gaussMutation (double t);

    /**
     * @brief
     * Méthode qui cherche un père aléatoirement dans le patch de la mère
     * si le mode de reproduction est l'allofécondation.
     *
     * Le père sera nécessairement différent de la mère.
     * Note importante: si la mère est seule dans son patch,
     * le programme rencontrera une boucle infinie.
     *
     * @param idPatch       patch de la mère
     * @param mother        identifiant de la mère
     *
     * @return              L'identifiant du père
     */
    int getFather(int patchMother, int mother);

    void calcNewRelationships(int whr);

    /**
     * @brief
     * Méthode qui affiche la progression à l'écran du terminal.
     *
     * @param progress  état de la progression
     */
     void printProgress(int progress);

    std::ofstream report; /**< @brief Variable permettant d'écrire le rapport */

    /**
     * @brief
     * Méthode qui écrit l'entête dans le rapport
     *
     * @param Kmin      La capacité d'accueil minimale.
     * @param Kmax      La capacité d'accueil maximale.
     * @param sigmaK    Le degré de varition de K dans l'espace.
     * @param Pmin      La probabilité maximale qu'un patch soit pollinisé.
     * @param Pmin      La probabilité minimale qu'un patch soit pollinisé.
     * @param sigmaP    Le degré de varition de P dans l'espace.
     */
    void writeHeader(int Kmin, int Kmax, int sigmaK, double Pmin, double Pmax, double sigmaP);

    /**
     * @brief
     * Méthode qui écrit un rapport pour un génération donnée
     *
     * @param gen   la génération pour laquelle il faut écrire le rapport
     */
    void writeReport(int gen);

    /**
     * @brief
     * Méthode qui permet de vider un vecteur et de libérer entièrement la mémoire
     */
    void clear_and_freeVector(std::vector<double>& toClear);
};

#endif // WORLD_H_INCLUDED
