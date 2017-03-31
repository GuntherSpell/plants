#ifndef WORLD_H_INCLUDED
#define WORLD_H_INCLUDED

#include <vector>
#include <array>
#include <fstream>

#include "patch.h"


/**
 * @file
 */

 /** @brief Énumération qui permet de choisir le type de distribution que l'on veut pour les mutations. */
typedef enum _mut_
{
    gaussian = 0,
    uniform = 1,
} distrMut;


/** @brief Structure qui permet de stocker la position relative (dans un patch) pour un individu. */
typedef struct _IndividualPosition_
{
    int posInPatch;
    int patch;
} IndividualPosition;

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
    World(int idWorld, int NPatch, double delta, double c, bool relationshipIsManaged,
          int typeMut, double mu, double sigmaZ, double d_s_relativeMutation, int Kmin, int Kmax, int sigmaK,
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

    /**
     * @brief Vecteur qui contient, pour chaque individu, le numéro de son patch et sa postion dans celui-ci.
     *
     * Ce vecteur permet de connaitre la position relative (dans un patch) d'un individu
     * à partir de sa posititon absolue (dans le monde).
     * Ce vecteur est très utile pour récupérer les infos de la mère après l'avoir tirée au sort.
     * Il sert également pour gérer l'apparentement entre individus.
     */
    std::vector<IndividualPosition> globalPop;

    double delta; /**< @brief La dépression de consanguinité */
    double c; /**< @brief Le coût de dispersion */

    bool relationshipIsManaged; /**< @brief Indique si on doit gérer l'apparentement */

    distrMut typeMut; /**< @brief La distribution de l'ampleur de mutation */
    double mu; /**< @brief La probabilité de mutation */
    double sigmaZ; /**< @brief L'ampleur de la mutation */
    double d_s_relativeMutation; /** @brief Mutation relative de d et s. Si égale à 1, seul d mute. */

    int NGen; /**< @brief Le nombre de générations à créer */
    int genReport; /**< @brief Le nombre de générations entre chaque rapport .txt */
    int genCount; /**< @brief Compteur de générations */

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
     * Pour les générations paires, les parents sont dans la case 0.
     * Pour les générations impaires, c'est l'inverse.
     */
    std::array<std::vector<std::vector<double>>, 2> relationship;

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

    std::ofstream report; /**< @brief Variable permettant d'écrire le rapport */

    /**
     * @brief
     * Méthode qui permet de retourner une valeur pour un patch selon sa position.
     *
     * @param minVal    Valeur maximale
     * @param maxVal    Valeur minimale
     * @param sigma     Degré de variation
     * @param i         La position du patch
     *
     * @return          La valeur pour le patch donné
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
     * @param mother        identifiant globale de la mère
     * @param autof         si la graine est issue d'autof ou non
     */
    void newInd(int whr, int mother, bool autof);

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
    double unifMutation(double t);

    /**
      * @brief
      * Crée une mutation selon une loi normale.
      *
      * @param t        Valeur du trait à muter.
      *
      * @return         La valeur du trait après mutation.
      */
    double gaussMutation(double t);

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

    /**
     * @brief
     * Méthode qui redéfinit l'état de pollinisation d'un patch
     * à partir de sa probabilité d'être pollinisé.
     *
     * @param p     La probabilité d'être pollinisé.
     *
     * @return      Si le patch est pollinisé ou non.
     */
    bool redefinePollination(double p);

    /** @brief Méthode qui va recalculer les apparentements entre tous les individus. */
    void calcNewRelationships(void);

    /**
     * @brief
     * Méthode qui affiche la progression à l'écran du terminal.
     *
     * @param progress  état de la progression
     */
     void printProgress(int progress);

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

    /** @brief Méthode qui écrit un rapport pour un génération donnée */
    void writeReport(void);

    /**
     * @brief
     * Méthode qui permet de vider un vecteur et de libérer entièrement la mémoire
     *
     * @param toClear   le vecteur à vider, qui doit contenir des doubles
     */
    void clear_and_freeVector(std::vector<double>& toClear);
};

#endif // WORLD_H_INCLUDED
