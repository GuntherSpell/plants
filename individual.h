#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <random>

/**
 * @file
 */

/**
 * @brief
 * Contient les caractéristiques d'un individu
 * et la manière dont il se reproduit
 */

class Individual
{
public:

    Individual(double s, double d, double f);

    double s; /**< @brief Le taux d'autofécondation */
    double d; /**< @brief Le taux de dispersion */

    double f; /**< @brief Le taux de consanguinité, utilisé si on gère l'apparentement */

    /**
     * @brief
     * Méthode qui calcule les pressions en propagules dipersantes générées par cet individu.
     *
     * @param delta         La dépression de consanguinité
     * @param c             Le coût de la dispersion
     * @param pollenized    L'état de pollinisation du patch de l'individu
     * @param press         Le vecteur de pression à remplir
     */
    void calcDispPress(double delta, double c, bool pollenized, std::vector<double>& press);

    /**
     * @brief
     * Méthode qui calcule les pressions en propagules résidentes générées par cet individu.
     *
     * @param delta         La dépression de consanguinité
     * @param c             Le coût de la dispersion
     * @param pollenized    L'état de pollinisation du patch de l'individu
     * @param press         Le vecteur de pression à remplir
     */
    void calcResidPress(double delta, bool pollenized, std::vector<double>& press);

private:

    /**
     * @brief
     * Méthode qui calcule la dépression de consanguinité
     * subie par l'individu à partir de son taux de consanguinité.
     *
     * @param delta La dépression de consanguinité maximale (si le taux de consanguinité est de 1).
     * @param f     Le taux de consanguinité de l'individu.
     *
     * @return      La dépression de consanguinité subie par l'individu.
     */
    double f_to_delta(double delta, double f);
};


#endif // INDIVIDUAL_H_INCLUDED

