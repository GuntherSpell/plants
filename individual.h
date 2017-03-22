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
     * Méthode qui calcule les pressions en
     * propagule générées par cet individu.
     *
     * Cette méthode ne génère que les pressions dispersantes
     * ou résidentes, selon les besoins du modèle.
     *
     * @param delta         La dépression de consanguinité
     * @param c             Le coût de la dispersion
     * @param pollenized    L'état de pollinisation du patch de l'individu
     * @param dispNeeded    Si on souhaite les pressions dispersantes ou résidentes
     * @param press         Le vecteur de pression à remplir
     */
    void calcPress(double delta, double c, bool pollenized, bool dispNeeded, std::vector<double>& press);
};


#endif // INDIVIDUAL_H_INCLUDED

