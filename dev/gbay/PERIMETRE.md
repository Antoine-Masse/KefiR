# Périmètre de g.bay()

Établi le 11 novembre 2026 à partir de `DescTools::GTest` 0.99.60 (code lu) et de la voie D17 de m.bay. Branche `claude/feat-g-bay-1f2987`, partie de `feat/m.bay` (c8e8976).

## Entrées

| Entrée | GTest | g.bay | Cas |
|---|---|---|---|
| `x` vecteur numérique d'effectifs, `y` absent | oui | oui | ajustement à `p` |
| `x` table à une dimension, ou matrice 1 × k / k × 1 | oui (ramenée à un vecteur) | oui, avec un « Changement de stratégie » | ajustement |
| `x` facteur ou texte, `y` absent | non (erreur) | oui : `table(x)` | ajustement |
| `x` matrice, table 2D, data.frame d'effectifs | oui | oui | indépendance, plan joint |
| `x` et `y` deux vecteurs de même longueur | oui | oui (manquants retirés, dit dans le récit) | indépendance, plan joint |
| `x` et `y` noms de colonnes avec `data` | non | oui | indépendance, plan joint |
| formule `~ A` avec `data` | non | oui | ajustement |
| formule `~ A + B` | non | oui | indépendance, plan joint |
| formule `n ~ A + B`, `n ~ A` (colonne d'effectifs, une ligne par combinaison) | non | oui | comme ci-dessus, pondéré |
| formule `y ~ g` (réponse qualitative) | non | oui | indépendance, plan indépendant, marge de g fixée |
| `p`, `rescale.p` | oui | oui, mêmes règles | ajustement |
| `correct` (`"williams"`, `"yates"`) | oui | accepté, ignoré, dit dans le récit (G1) | — |
| table à 3 dimensions ou plus | erreur | refus en clair | — |

## Réglages propres à g.bay

`sampling` (`"auto"`, `"jointMulti"`, `"indepMulti"`, `"poisson"`, `"hypergeom"`), `fixed` (`"rows"`, `"cols"`), `prior_h1`, `sensitivity`, `draws`, `seed`, `verbose`, `return`, `detail`. Aucun ne s'ajuste sur les données.

## Sorties

| GTest | g.bay |
|---|---|
| statistique G, ddl, p-value | facteur de Bayes, probabilité a posteriori de chaque hypothèse, qualificatif ; jamais de p-value |
| `expected` | `expected` (mêmes valeurs) |
| — | rapport à l'indépendance par case (ou à la proportion théorique par modalité) avec intervalle de crédibilité (G6) |
| — | V de Cramér a posteriori (G7), w de Cohen pour l'ajustement |
| — | rapport des cotes et intervalle de crédibilité en 2 × 2 (G4) |
| — | sensibilité aux a priori et au plan d'échantillonnage |
| — | récit en sept étapes, bilan, phrase de compte rendu, en français et en anglais |

## Hors périmètre (refus ou renvoi)

- Tables à trois dimensions et plus (Cochran-Mantel-Haenszel, log-linéaire) : refus.
- Données appariées (McNemar) : non détectables depuis une table ; la page d'aide le signale.
- Quels groupes diffèrent deux à deux (post-hoc par regroupements) : c'est `m.bay(y ~ g)`, que la page d'aide indique.
- Ordre des modalités : ignoré, comme dans `GTest` (une réponse ordonnée relève de `m.bay`).
- Cases structurellement nulles (quasi-indépendance) : non traitées ; une ligne ou colonne entièrement vide est retirée (G3).
