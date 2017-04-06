#ifndef ENGINE_H_INCLUDED
#define ENGINE_H_INCLUDED

#include <vector>
#include <array>

/**
 * @file
 */

/**
 * @brief
 * Fonction qui prend les paramètres du fichier texte config.txt
 *
 * Pour plus d'informations, regarder le fichier config.txt
 *
 * @param param         Vecteur où mettre les paramètres du (des) monde(s)
 * @param NWorld        Valeur du nombre de mondes à modifier
 * @param NReplicats    Valeur du nombre de réplicats à modifier
 */
void getParamsFromTxt(std::vector<std::array<double, 23>>& params, int& NWorld, int& NReplicats);

/** @brief Fonction qui lance la simulation pour tous les mondes */
void runSimu(void);

#endif // ENGINE_H_INCLUDED
