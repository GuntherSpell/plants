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

    std::array<double, 2> previous_d_means; /**< @brief 2 dernières valeurs de la moyenne de d. */
    std::array<double, 2> previous_d_vars; /**< @brief 2 dernières valeurs de la variance de d. */
    std::array<double, 2> previous_s_means; /**< @brief 2 dernières valeurs de la moyenne de s. */
    std::array<double, 2> previous_s_vars; /**< @brief 2 dernières valeurs de la variance de d. */

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

    int check_convergence(int reportCount);

    bool check_stats(std::array<double, 2> previous_means, std::array<double, 2> previous_vars, std::array<double, 2> new_vals);

    void calc_mean_var(std::vector<double> values, std::array<double, 2>& mean_var);
};

#endif // PATCH_H_INCLUDED
