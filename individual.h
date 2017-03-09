#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <vector>

/**
 * @file
 */


typedef enum _mut_
{
    gaussian = 0,
    uniform = 1,
} distrMut;

/**
 * @brief
 * Contient les caractéristiques d'un individu,
 * la manière dont il se reproduit
 * et les modalités de la mutation
 */

class Individual
{
public:

    Individual(double s, double d);

    double s; /**< @brief Le taux d'autofécondation */
    double d; /**< @brief Le taux de dispersion */

    /**
     * @brief
     * Méthode qui applique une mutation aléatoire
     * sur un des traits de l'individu
     *
     * Actuellement, l'ampleur de la mutation suit
     * une loi normale centrée.
     *
     * @param mu        La probabilité de mutation
     * @param sigmaZ    L'ampleur de la mutation.
     */
    void mutation(double mu, double sigmaZ, distrMut typeMut);

    /**
     * @brief
     * Méthode qui calcule les pressions en
     * propagule générées par cet individu.
     *
     * Cette méthode ne génère que les pressions
     * dispersantes ou résidentes, selon
     * les besoins du modèle
     *
     * @param delta         La dépression de consanguinité
     * @param c             Le coût de la dispersion
     * @param pollenized    L'état de pollinisation du patch de l'individu
     * @param disp          Si on souhaite les pressions dispersantes ou résidentes
     * @param press         Le vecteur de pression à remplir
     */
    void calcPress(double delta, double c, bool pollenized, bool disp, std::vector<double>& press);

private:

     /**
      * @brief
      * Crée une mutation selon une loi uniforme.
      *
      * @param sigmaZ   L'ampleur de la mutation.
      * @param t        Valeur du trait à muter.
      */
    double unifMutation (double sigmaZ, double t);

    /**
      * @brief
      * Crée une mutation selon une loi normale.
      *
      * @param sigmaZ   L'ampleur de la mutation.
      * @param t        Valeur du trait à muter.
      */
    double gaussMutation (double sigmaZ, double t);
};


#endif // INDIVIDUAL_H_INCLUDED

