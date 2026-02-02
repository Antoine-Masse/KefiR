# Guide d'installation - Intégration medpb2

**Date:** 2025-01-04
**Package:** KefiR v02 bêta

---

## PROBLÈME RENCONTRÉ

Vous avez testé `m.test()` et obtenu :
```
c) Posthoc - Comparaisons robustes de médianes [medpb2() de {WRS2}].
    Attention ! Test medpb2 par paires échoué. Utilisation des tests précédents.
```

**Cause:** Le fichier `R/sys_pairwise_medpb2.R` n'est pas chargé dans votre session R.

---

## SOLUTION 1 : Recharger le package (RECOMMANDÉ)

### Option A : avec devtools (recommandé)

```r
# Installer devtools si nécessaire
if (!require(devtools)) install.packages("devtools")

# Se placer dans le répertoire du package
setwd("/chemin/vers/KefiR")

# Recharger TOUT le package
devtools::load_all(".")

# Tester
x <- rexp(90, rate = c(1, 0.5, 0.3))
g <- factor(rep(c("A", "B", "C"), each = 30))
m.test(x ~ g, verbose = TRUE, plot = FALSE)
```

### Option B : Réinstaller le package

```r
# Depuis le répertoire parent de KefiR
setwd("/chemin/vers/parent")

# Désinstaller l'ancienne version
remove.packages("KefiR")

# Réinstaller
install.packages("KefiR", repos = NULL, type = "source")

# Recharger
library(KefiR)
```

---

## SOLUTION 2 : Chargement manuel (temporaire)

Si vous ne pouvez pas recharger le package immédiatement :

```r
# Charger manuellement le nouveau fichier
source("R/sys_pairwise_medpb2.R")

# Vérifier que la fonction est bien chargée
exists(".pairwise_medpb2")  # Doit retourner TRUE

# Tester
x <- rexp(90, rate = c(1, 0.5, 0.3))
g <- factor(rep(c("A", "B", "C"), each = 30))
m.test(x ~ g, verbose = TRUE, plot = FALSE)
```

---

## VÉRIFICATION APRÈS INSTALLATION

### Test 1 : Vérifier que la fonction existe

```r
exists(".pairwise_medpb2", mode = "function")
# Doit retourner: TRUE
```

### Test 2 : Tester la fonction directement

```r
library(WRS2)

set.seed(123)
x <- c(rnorm(20, 0), rnorm(20, 1), rnorm(20, 2))
g <- factor(rep(c("A", "B", "C"), each = 20))

# Appel direct
result <- .pairwise_medpb2(x, g, alpha = 0.05, nboot = 100, debug = FALSE)
print(result$groups)

# Doit afficher:
#   categories medpb2_Holm
# 1          A           a
# 2          B           b
# 3          C           c
```

### Test 3 : Tester via m.test()

```r
# Données avec forte asymétrie (pour déclencher med1way)
set.seed(456)
x <- rexp(90, rate = c(1, 0.5, 0.3))
g <- factor(rep(c("A", "B", "C"), each = 30))

# Lancer m.test
result <- m.test(x ~ g, verbose = TRUE, plot = FALSE, return = TRUE)

# Vérifier les colonnes
print(colnames(result$posthoc$groups))

# Doit contenir: "categories" "Bootstrap" "Wilcoxon_Holm" "medpb2_Holm"
```

---

## SI LE PROBLÈME PERSISTE

### Diagnostic approfondi

```r
# 1. Vérifier les fichiers R/
list.files("R", pattern = "medpb2")
# Doit retourner: "sys_pairwise_medpb2.R"

# 2. Vérifier que WRS2 est installé
if (!require(WRS2)) install.packages("WRS2")

# 3. Tester medpb2 directement (2 groupes)
library(WRS2)
x <- c(rnorm(20), rnorm(20, 1))
g <- factor(rep(c("A", "B"), each = 20))
medpb2(x ~ g, nboot = 100)
# Doit fonctionner sans erreur

# 4. Sourcer manuellement et tester
source("R/sys_pairwise_medpb2.R")
source("R/sys_dbg.R")

x <- c(rnorm(20, 0), rnorm(20, 1), rnorm(20, 2))
g <- factor(rep(c("A", "B", "C"), each = 20))
result <- .pairwise_medpb2(x, g, nboot = 100, debug = TRUE)
print(result)
```

### Message d'erreur attendu (si fonction non chargée)

Si vous relancez `m.test()` AVEC le message d'erreur détaillé ajouté :

```
ERROR in .pairwise_medpb2(): .pairwise_medpb2() function not found.
Did you source R/sys_pairwise_medpb2.R?
```

→ Cela confirme que vous devez recharger le package (Solution 1).

---

## RÉSULTAT ATTENDU APRÈS CORRECTION

Quand tout fonctionne correctement, vous devriez voir :

```
10) Tests post-hocs non paramétriques de comparaison des groupes.
    a) Posthoc - Analyse des différences de médianes par bootstrap [pairwise.boot() de {KefiR}].
    b) Posthoc - Test de Wilcoxon-Mann-Whitney [pairwise.wilcox.test()].
    c) Posthoc - Comparaisons robustes de médianes [medpb2() de {WRS2}, correction de Holm].
        Méthode bootstrap percentile (100 itérations).
        Note : medpb2 est la méthode recommandée pour comparaisons
        de médianes avec valeurs ex-aequo et forte asymétrie (Wilcox, 2017).
```

**SANS le message d'erreur "Attention ! Test medpb2 par paires échoué"**

---

## FICHIERS MODIFIÉS À RECHARGER

Pour que medpb2 fonctionne, ces fichiers doivent être chargés :

1. ✅ `R/m.test.R` (ligne 175) - import WRS2::medpb2
2. ✅ `R/sys_posthoc.R` (lignes 1750-1817) - intégration medpb2
3. ✅ `R/sys_pairwise_medpb2.R` (NOUVEAU) - wrapper medpb2
4. ✅ `R/sys_dbg.R` (existant) - fonction debug

---

## COMMANDE RAPIDE (tout-en-un)

```r
# Nettoyage complet et rechargement
rm(list = ls())
if (require(devtools)) {
  setwd("/chemin/vers/KefiR")
  devtools::load_all(".")
} else {
  source("R/sys_pairwise_medpb2.R")
  source("R/sys_posthoc.R")
  source("R/sys_dbg.R")
  library(KefiR)
}

# Test immédiat
set.seed(123)
x <- rexp(90, rate = c(1, 0.5, 0.3))
g <- factor(rep(c("A", "B", "C"), each = 30))
m.test(x ~ g, verbose = TRUE, plot = FALSE)
```

---

## CONTACT

Si le problème persiste après avoir suivi ce guide :

**Email:** antoine.masse@u-bordeaux.fr

**Inclure dans le message:**
1. Résultat de `exists(".pairwise_medpb2")`
2. Résultat de `list.files("R", pattern = "medpb2")`
3. Méthode utilisée (devtools::load_all vs source manuel)
4. Message d'erreur complet
5. Version de R : `R.version.string`
6. Packages installés : `sessionInfo()`

---

**FIN DU GUIDE**
