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
  * Structure qui permet de choisir
  * le type de distribution que l'on
  * veut pour les mutations.
  */

/**
 * @brief
 * Contient les caractéristiques d'un monde.
 * Cette classe contient également les méthodes
 * permettant l'exécution du modèle.
 */

class World
{
public:

    World(int id, int NPatch, double delta, double c, int typeMut, double mu, double sigmaZ, int Kmin, int Kmax, int sigmaK,
    double Pmin, double Pmax, double sigmaP, double sInit, double dInit, int NGen, int genReport);

    void run(void); /**< @brief Méthode qui lance la simulation */

private:

    /**
     * @brief Identifiant du monde.
     *
     * Utile uniquement pour générer des rapports
     * différents s'il y a plusieurs mondes.
     */
    int id;

    int NPatch; /**< @brief Nombre de patchs du monde */

    std::vector<Patch> patches; /**< @brief Vecteur qui contient tous les patchs du monde */

    double delta; /**< @brief La dépression de consanguinité */
    double c; /**< @brief Le coût de dispersion */

    distrMut typeMut; /**< @brief La distribution de l'ampleur mutation */
    double mu; /**< @brief La probabilité de mutation */
    double sigmaZ; /**< @brief L'ampleur de la mutation */

    int Kmin; /**< @brief La capacité d'accueil minimale */
    int Kmax; /**< @brief La capacité d'accueil maximale */
    int sigmaK; /**< @brief Le degré de varition de K dans l'espace */

    double Pmin;/**< @brief La probabilité maximale qu'un patch soit pollinisé */
    double Pmax; /**< @brief La probabilité minimale qu'un patch soit pollinisé */
    double sigmaP; /**< @brief Le degré de varition de P dans l'espace */

    double sInit; /**< @brief La valeur initiale du taux d'autofécondation */
    double dInit; /**< @brief La valeur initiale du taux de disperion */

    int NGen; /**< @brief Le nombre de générations à créer */
    int genReport; /**< @brief Le nombre de générations entre chaque rapport .txt */

    std::default_random_engine generator; /**< @brief Générateur de nombre aléatoire */


    /**
     * @brief Deux vecteurs qui contiennent temporairement la nouvelle génération d'un patch
     *
     * Pour plus de détails, voir la méthode createNextGen
     */
    std::array<std::vector<Individual>,2> tmp;

    /**
     * @brief
     * Méthode qui permet de retourner une valeur pour un patch selon sa position.
     *
     * @param minVal    Valeur maximale
     * @param maxVal    Valeur minimale
     * @param sigma     Degré de variation
     * @param i         La position du patch
     */
     double distr(double minVal, double maxVal, double sigma, int i);

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
     * @param id    L'identifiant du patch dont on souhaite créer
     *              la nouvelle génération.
     */
    void createNextGen(int id);

    /**
     * @brief
     * Méthode qui crée un nouvel individu selon la propagule choisie
     *
     * @param whr           indique dans quel vecteur temporaire il faut
     *                      stocker la génération
     * @param patchMother   identifiant du patch de la mère
     * @param mother        identifiant de la mère
     * @param autof         si la graine est issue d'autof ou non
     */
    void newInd(int whr, int patchMother, int mother, bool autof);

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
     * @param fatherTraits  vecteur qui stocke les traits du père
     */
    void getFather(int idPatch, int mother, std::array<double,2>& fatherTraits);

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
     * @param report    le fichier où écrire le rapport
     */
    void writeHeader(void);

    /**
     * @brief
     * Méthode qui écrit un rapport pour un génération donnée
     *
     * @param gen       la génération pour laquelle il faut écrire le rapport
     * @param report    le fichier où écrire le rapport
     */
    void writeReport(int gen);
};

#endif // WORLD_H_INCLUDED
