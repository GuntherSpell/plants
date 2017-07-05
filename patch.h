#ifndef PATCH_H_INCLUDED
#define PATCH_H_INCLUDED

#include <vector>
#include <array>

#include "individual.h"

/**
 * @file
 */


/**
 * @brief
 * Contient les caractéristiques d'un patch et les
 * méthodes nécessaires au fonctionnement du modèle.
 */

class Patch
{
public:

    Patch(double p, int K, double sInit, double dInit, int pos_of_first_ind);

    int K; /**< @brief La capacité d'accueil du patch */

    /** @brief Un vecteur qui contient tous les individus du patch */
    std::vector<Individual> population;

    double p; /**< @brief La probabilité d'être pollinisé */

    bool pollenized; /**< @brief L'état de pollinisation */

    /** @brief Si le taux de dispersion des individus du patch est stable dans le temps. */
    bool d_hasConverged;
    /** @brief Si le taux d'autofécondation des individus du patch est stable dans le temps. */
    bool s_hasConverged;

    std::vector<double> previous_d_means; /**< @brief n dernières valeurs de la moyenne de d */
    std::vector<double> previous_s_means; /**< @brief n dernières valeurs de la moyenne de s */

    /**
     * @brief
     * Matrice qui indique si les états précédents de la moyenne de s du patch sont similaires ou non.
     */
    std::vector<std::vector<bool>> prev_gens_s_similarity_matrix;

    /**
     * @brief
     * Matrice qui indique si les états précédents de la moyenne de d du patch sont similaires ou non.
     */
    std::vector<std::vector<bool>> prev_gens_d_similarity_matrix;

    /**
    * @brief La position absolue (dans le monde entier) du premier individu du patch.
    *
    * C'est une info utile pour retrouver la position absolue de n'importe quel individu du patch.
    */
    int pos_of_first_ind;

    /**
     * @brief
     * Vecteur qui contient les pressions en graines dispersantes.
     * Cela permet de ne pas devoir les calculer 2 fois quand
     * on est au patch à gauche puis à droite.
     */
    std::vector<double> dispSeeds;

    /**
     * @brief
     * Méthode qui détermine si le patch est pollinisé
     * en fonction de la probabilité qu'il le soit
     */
    void isPollenized(void);

    /**
     * @brief
     * Méthode qui rassemble toutes les pression de
     * propagules dispersantes de tous les individus du patch.
     *
     * @param delta         La dépression de consanguinité
     * @param c             Le coût de dispersion
     * @param press         Le vecteur de pression à remplir
     */
    void getDispPress(double delta, double c, std::vector<double>& press);

    /**
     * @brief
     * Méthode qui rassemble toutes les pression de
     * propagules résidentes de tous les individus du patch.
     *
     * @param delta         La dépression de consanguinité
     * @param press         Le vecteur de pression à remplir
     */
    void getResidPress(double delta, std::vector<double>& press);

    /**
     * @brief
     * Méthode qui vérifie l'état de convergence du patch par rapport aux deux états précédents
     *
     * @param checkCount            Le nombre de fois que la convergence a été vérifiée
     * @param NGenToConverge        Le nombre de vérifications consécutives identiques pour considérer que le patch a convergé
     * @param n_choose_2            Le nombre de combinaisons de 2 parmi les générations à vérifier pour la convergence
     * @param relativeConvergence   Le critère de variation relative pour juger de l'état de convergence
     * @param absoluteConvergence   Le critère de variation absoule pour juger de l'état de convergence
     *
     * @return                      1 si le patch a convergé, 0 sinon. Utile pour facilement connaitre le nbr de patchs convergés.
     */
    int check_convergence(int checkCount, int NGenToConverge, int n_choose_2, double relativeConvergence, double absoluteConvergence);

    /**
     * @brief
     * Méthode qui compare deux moyennes et juge si elles sont suffisamment similaires selon deux critères.
     *
     * @param first_mean            La première valeur à comparer
     * @param second_mean           La seconde valeur à comparer
     * @param relativeConvergence   Le critère de variation relative pour juger si les moyennes sont similaires
     * @param absoluteConvergence   Le critère de variation absoule pour juger si les moyennes sont similaires
     *
     * @return                      Vrai si les moyennes sont suffisamment similaires, faux sinon
     */
    bool check_stats(double first_mean, double second_mean, double relativeConvergence, double absoluteConvergence);

    /**
     * @brief
     * Méthode qui calcule la moyenne d'un ensemble de valeurs numériques.
     *
     * @param values    Vecteur qui contient les valeurs
     *
     * @return          La valeur de la moyenne
     */
    double calc_mean(std::vector<double> values);
};

#endif // PATCH_H_INCLUDED
