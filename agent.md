# Agent KefiR — Rôle et périmètre d’intervention

Le package **KefiR** doit être publié sur GitHub (nouvelle mise à jour à venir).

---

## 1. Principes généraux

- Modifier **le moins possible** les fichiers existants.
- Créer **le moins possible** de nouveaux fichiers.
- Conserver **le style, la structure, la documentation et la logique actuelle**.
- Les interventions peuvent concerner :
  - la documentation roxygen2,
  - DESCRIPTION / NAMESPACE,
  - métadonnées du package,
  - organisation interne,
  - scripts de test,
  - **sans forcément modifier un script R**.

---

## 2. Règles de développement — bp.log (AUTORITÉ)

Les règles officielles du projet sont définies dans :
dev/bp.log à lire nécessairement et régulièrement dans un contexte de modfication de scripts R.

Tu as le droit d'utiliser R pour tester les commandes. Si tu génères des scripts de tests tu les sauvegardes dans KefiR/dev et le nom du fichier doit commencer par "agent_"
Dans dev, tu crée un sous-dossier save où tu sauvegardes les différentes versions des scripts R modifier avec un historique de leur modification optimisé pour un retour en arrière.

Tu dois aussi créer un dossier test dans dev où tu crées des scripts R pour tester les fonctions que tu modifies. Ces scripts doivent être nommés de façon à ce qu'ils soient faciles à comprendre et à utiliser. Ils doivent aussi être commentés de façon à ce qu'ils soient faciles à comprendre et à utiliser.

Dans dev, load_all_kefir.R permet de tout charger comme si KefiR était installé afin de le tester.

.testeur.R est un script qui permet de tester les fonctions de KefiR. Il est commenté de façon à ce qu'il soit facile à comprendre et à utiliser et l'ensemble de ses tests (sauf ceux qui testent les messages error doivent donner des résultats satisfaisants

## 3. Règles de préparation du package : autorisation à utilise

Autorisation pour exécuter la commande devtools::check() sous R.
