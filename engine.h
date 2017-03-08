#ifndef ENGINE_H_INCLUDED
#define ENGINE_H_INCLUDED

#include <vector>

/**
 * @file
 */

/**
 * @brief
 * Méthode que prend les paramètres du fichier texte config.txt
 *
 * Pour plus d'informations, regarder le fichier config.txt
 *
 * @param param vecteur où mettre les paramètres du (des) monde(s)
 */
int getParamsFromTxt(std::vector<double>& params);

/**
 * @brief
 * Méthode qui lance la simulation pour tous les mondes
 */
void runSimu(void);

#endif // ENGINE_H_INCLUDED
