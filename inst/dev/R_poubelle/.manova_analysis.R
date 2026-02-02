#* Note perso :
#* Ce qu'il faut faire sur cette fonction...
#* 1) Lire tout ce qu'il a dans le présent fichier et construire le pipeline approprié
#* ... pour traiter les MANOVA dans tous les scénarios (tailles, effet aléatoire, situations paramétriques et non paramétriques).
#* 2) S'inspirer de .multi_factor.analysis() et .one_factor_analysis() pour utiliser le même style, mettre en place
#* ... les bons contrôles des assomptions et les verbose .vbse()
#* 3) A chaque fois, vérifier que les contrôles sont acadédmiques - mettre en commentaire des sources vérifiables (idée, APA, lien)
#* 4) Rédiger l'entête de la fonction et disperser cette nopte dans les commentaires.
#* 5) Utiliser le code .testeur.R quitte à le mettre à jour.







.manova_analysis <- function(x=NULL,g=NULL,formula=NULL,data=NULL,
	paired = FALSE, id = NULL,alpha = 0.05,
	k=NULL,code=FALSE,debug=FALSE,verbose=FALSE,boot=TRUE, silent=TRUE, return=TRUE) {

	###################
	#	Un code à conserver ?
	###################
	if (!is.null(formula)) {
	  print("DEBUG: Début du traitement")

	  clean_formula <- function(expr) {
		if (is.call(expr)) {
		  func_name <- as.character(expr[[1]])
		  transform_funcs <- c("I", "log", "sqrt", "exp", "scale", "as.numeric", "as.factor")
		  preserve_funcs <- c("cbind", "interaction", "poly", "bs", "ns")

		  print(paste("DEBUG: Traitement fonction:", func_name))

		  if (func_name %in% transform_funcs) {
			return(clean_formula(expr[[2]]))
		  } else if (func_name %in% preserve_funcs) {
			return(expr)
		  } else if (func_name == ":") {
			return(expr)
		  }
		  return(as.call(lapply(expr, clean_formula)))
		}
		return(expr)
	  }

	  print("DEBUG: Nettoyage formule")
	  cleaned_formula <- clean_formula(formula)
	  vars_needed <- all.vars(cleaned_formula)
	  print(paste("DEBUG: Variables nécessaires:", paste(vars_needed, collapse=", ")))

	  response_expr <- formula[[2]]
	  original_response_expr <- response_expr
	  cleaned_response_expr <- clean_formula(response_expr)

	  print(paste("DEBUG: Expression réponse originale:", deparse(original_response_expr)))

	  response_vars <- if (is.call(response_expr)) {
		func_name <- as.character(response_expr[[1]])
		if (func_name == "cbind") {
		  unique(unlist(lapply(response_expr[-1], all.vars)))
		} else {
		  all.vars(cleaned_response_expr)
		}
	  } else {
		all.vars(cleaned_response_expr)
	  }

	  print(paste("DEBUG: Variables réponse:", paste(response_vars, collapse=", ")))

	  response <- if (all(response_vars %in% colnames(data))) {
		print("DEBUG: Réponse dans data")
		if (length(response_vars) == 1) {
		  if (is.call(original_response_expr)) {
			print("DEBUG: Évaluation transformation sur colonne data")
			eval(original_response_expr, envir = data, enclos = parent.frame())
		  } else {
			data[[response_vars]]
		  }
		} else {
		  if (is.call(original_response_expr) && as.character(original_response_expr[[1]]) == "cbind") {
			eval(original_response_expr, envir = data, enclos = parent.frame())
		  } else {
			data[, response_vars, drop = FALSE]
		  }
		}
	  } else {
		print("DEBUG: Réponse externe à data")
		tryCatch({
		  if (is.call(original_response_expr)) {
			print("DEBUG: Évaluation expression externe")
			eval(original_response_expr, envir = list2env(as.list(data), parent = parent.frame()))
		  } else {
			get(response_vars, envir = parent.frame())
		  }
		}, error = function(e) {
		  print(paste("DEBUG: Erreur lors de l'évaluation:", e$message))
		  stop(e)
		})
	  }

	  # Suite du code identique...
	  g_vars <- setdiff(vars_needed, response_vars)
	  print(paste("DEBUG: Variables explicatives:", paste(g_vars, collapse=", ")))

	  g_list <- lapply(g_vars, function(var) {
		if (var %in% colnames(data)) {
		  data[[var]]
		} else {
		  get(var, envir = parent.frame())
		}
	  })

	  g <- if (length(g_list) > 1) {
		as.data.frame(setNames(g_list, g_vars))
	  } else if (length(g_list) == 1) {
		g_list[[1]]
	  } else {
		NULL
	  }
	}

}

############################
# Une piste :
1.	MANOVA classique
r
Copy
# Pour variables dépendantes Y1, Y2, Y3 et facteurs A, B
fit <- manova(cbind(Y1, Y2, Y3) ~ A * B)
summary(fit, test="Wilks") # Autres tests: "Pillai", "Hotelling-Lawley", "Roy"
2.	MANOVA avec effet aléatoire
r
Copy
library(nlme)
# random = ~1|Subject pour effet aléatoire sur Subject
fit <- lme(cbind(Y1, Y2, Y3) ~ A * B, random = ~1|Subject)
anova(fit)
3.	MANOVA mesures répétées
r
Copy
library(car)
# Pour données en format long
fit <- Anova(lm(DV ~ IV * Time + Error(Subject/Time)),
             idata = data.frame(Time), 
             idesign = ~Time)
4.	MANOVA robuste
r
Copy
library(rrcov)
fit <- Wilks.test(Y ~ group, data = data, method = "MCD")
# Ou avec CovRob
library(CovRob)
fit <- robustMANOVA(Y ~ group, data = data)




1.	Normalité Multivariée La normalité multivariée est une condition essentielle pour la MANOVA. Elle stipule que le vecteur des variables dépendantes doit suivre une distribution normale multivariée  => BHEP ? Faut-il faire le BHEP ?
1
. Cette condition est souvent vérifiée à l'aide du test de Shapiro-Wilk pour chaque variable-réponse à chaque niveau de la variable de groupement. Si les données ne suivent pas une distribution normale, cela peut affecter la validité des résultats de la MANOVA.  cf. pages envoyés par mails.
2.	Homogénéité des Matrices de Covariance Les matrices de variance-covariance des variables dépendantes doivent être égales entre les groupes 
2
. Cette condition est vérifiée à l'aide du Test M de Box, qui est très sensible et dont la significativité est généralement déterminée à un niveau alpha de 0,001. Si cette condition est violée, cela peut indiquer que les groupes ont des variances différentes, ce qui peut biaiser les résultats de l'analyse.
3.	Indépendance des Observations Chaque observation doit être indépendante des autres 
3
. Cela signifie que chaque sujet ne doit appartenir qu'à un seul groupe, et il ne doit pas y avoir de lien entre les observations de chaque groupe. Cette condition est cruciale pour éviter les biais dans les résultats.
4.	Absence de Multicollinéarité Les variables dépendantes ne doivent pas être trop corrélées entre elles 
4
. Une corrélation trop élevée (généralement supérieure à r = 0,90) peut indiquer une multicollinéarité, ce qui peut rendre difficile l'interprétation des résultats de la MANOVA.
5.	Linéarité Il doit y avoir une relation linéaire entre toutes les variables-réponses pour chaque groupe 
5
. Cette condition assure que les relations entre les variables sont correctement modélisées par la MANOVA.
6.	Absence de Valeurs Aberrantes Les valeurs aberrantes, qu'elles soient univariées ou multivariées, peuvent influencer de manière significative les résultats de la MANOVA 
6
. Les valeurs aberrantes univariées peuvent être identifiées à l'aide de boxplots, tandis que les valeurs aberrantes multivariées peuvent être détectées à l'aide de la distance de Mahalanobis.
7.	Échantillonnage Aléatoire Les données doivent être collectées par un processus d'échantillonnage aléatoire 
7
. Cela garantit que l'échantillon est représentatif de la population, permettant ainsi de généraliser les résultats de l'échantillon à la population.
Pour vérifier ces conditions, plusieurs tests statistiques sont couramment utilisés :
•	Le test de Shapiro-Wilk pour la normalité univariée
•	Le test M de Box pour l'homogénéité des matrices de covariance
•	La distance de Mahalanobis pour détecter les valeurs aberrantes multivariées
•	Le test de Levene pour l'homogénéité des variances 
2
•	Les graphiques Q-Q pour évaluer la normalité multivariée
•	Les matrices de nuages de points pour vérifier la linéarité
Il est important de noter que si ces conditions ne sont pas remplies, des ajustements ou des alternatives peuvent être envisagés. Par exemple, on peut utiliser des transformations de données, des méthodes robustes moins sensibles aux violations des hypothèses, ou encore des approches non paramétriques comme le test de Kruskal-Wallis multivarié.En conclusion, la vérification minutieuse de ces conditions est essentielle pour garantir la validité et la fiabilité des résultats de la MANOVA. Une attention particulière doit être portée à chaque condition pour s'assurer que les conclusions tirées de l'analyse sont précises et interprétables correctement.


La normalité multivariée implique un contrôle de l'ensemble des données des variables dépendantes, plutôt qu'un simple test de normalité par variable dépendante. Cependant, il est important de noter que la vérification de la normalité multivariée est plus complexe et implique des considérations supplémentaires par rapport à la normalité univariée. Voici une explication détaillée :
1.	Définition de la normalité multivariée La normalité multivariée est une extension de la loi normale unidimensionnelle à plusieurs dimensions. Elle concerne la distribution conjointe de plusieurs variables aléatoires 
1
. Cela signifie que non seulement chaque variable individuelle doit suivre une distribution normale, mais aussi que toute combinaison linéaire de ces variables doit également être normalement distribuée.
2.	Contrôle de l'ensemble des données La vérification de la normalité multivariée nécessite un examen de l'ensemble des données des variables dépendantes simultanément. Cela est dû au fait que la normalité multivariée prend en compte les relations et les interdépendances entre les variables, ce qui ne peut pas être capturé par des tests de normalité univariés individuels 
2
.
3.	Insuffisance des tests univariés Bien que les tests de normalité univariés (comme le test de Shapiro-Wilk) puissent être appliqués à chaque variable dépendante individuellement, cela n'est pas suffisant pour établir la normalité multivariée. Il est possible que chaque variable soit normalement distribuée individuellement, mais que leur distribution conjointe ne soit pas normale 
3
.
4.	Méthodes de vérification spécifiques Pour vérifier la normalité multivariée, des méthodes et tests spécifiques sont utilisés :a. Test de Mardia : Ce test évalue la normalité multivariée en examinant la skewness et la kurtosis multivariées 
4
.b. Test de Royston : Une extension du test de Shapiro-Wilk pour les données multivariées 
5
.c. Test de Henze-Zirkler : Un test basé sur une mesure de distance entre la fonction caractéristique empirique et la fonction caractéristique d'une distribution normale multivariée 
6
.d. Méthodes graphiques : Les parcelles Q-Q multivariées et les matrices de nuages de points peuvent être utilisées pour une évaluation visuelle de la normalité multivariée 
4
.
5.	Implications pratiques La vérification de la normalité multivariée est cruciale pour de nombreuses analyses statistiques multivariées, telles que MANOVA, PCA, et l'analyse canonique. Elle affecte directement la validité des tests statistiques, la précision des interprétations et le contrôle du taux d'erreur de type I 
7
.
6.	Complexité accrue La vérification de la normalité multivariée est intrinsèquement plus complexe que la vérification de la normalité univariée. Elle nécessite des techniques statistiques plus avancées et une interprétation plus nuancée des résultats 
8
.==> Penser à la fin à une analyse discriminante. Et se renseigner sur les tests post-hocs multivariés. A priori, on ne fait pas de standardisation sauf si grand effet d’une variable.
En conclusion, bien que la vérification de la normalité pour chaque variable dépendante individuellement soit une étape utile, elle n'est pas suffisante pour établir la normalité multivariée. Un contrôle de l'ensemble des données des variables dépendantes est nécessaire pour évaluer correctement la normalité multivariée, en utilisant des méthodes spécifiquement conçues pour cette tâche.
Pour effectuer une analyse de normalité multivariée sous R, vous pouvez utiliser le package MVN, qui offre une gamme complète d'outils pour tester et visualiser la normalité multivariée. Voici un guide étape par étape pour réaliser cette analyse :
1.	Installation et chargement du package
Commencez par installer et charger le package MVN :
R
# Installation du package (à faire une seule fois)
install.packages("MVN")

# Chargement du package
library(MVN)
2.	Préparation des données
Assurez-vous que vos données sont dans un format approprié. Par exemple, vous pouvez utiliser le jeu de données iris pour la démonstration :
R
# Chargement du jeu de données iris
data(iris)

# Sélection des variables numériques
iris_data <- iris[, 1:4]
3.	Test de normalité multivariée
Vous pouvez effectuer plusieurs tests de normalité multivariée. Voici comment réaliser les tests de Mardia, Henze-Zirkler et Royston :
R
# Test de Mardia
result_mardia <- mvn(data = iris_data, mvnTest = "mardia")
print(result_mardia$multivariateNormality)

# Test de Henze-Zirkler
result_hz <- mvn(data = iris_data, mvnTest = "hz")
print(result_hz$multivariateNormality)

# Test de Royston
result_royston <- mvn(data = iris_data, mvnTest = "royston")
print(result_royston$multivariateNormality)
4.	Interprétation des résultats
Pour chaque test, examinez la valeur p (p-value). Si la valeur p est inférieure à 0,05, cela suggère que les données ne suivent pas une distribution normale multivariée 
1
 
2
.
5.	Visualisation graphique
La méthode graphique la plus courante pour évaluer la normalité multivariée est le QQ-plot multivarié :
R
# Création d'un QQ-plot multivarié
mvn(data = iris_data, multivariatePlot = "qq")
Interprétez le QQ-plot : si les points s'écartent significativement de la ligne y = x, cela indique des écarts par rapport à la normalité multivariée 
3
.
6.	Analyse complète
Pour une analyse plus complète, vous pouvez combiner les tests statistiques et les méthodes graphiques en une seule commande :
R
result_complet <- mvn(data = iris_data, mvnTest = "mardia", multivariatePlot = "qq")
Cette commande effectuera le test de Mardia et générera un QQ-plot multivarié 
4
.En suivant ces étapes, vous pouvez effectuer une analyse approfondie de la normalité multivariée de vos données sous R. N'oubliez pas que la normalité multivariée est une condition plus stricte que la normalité univariée, et il est important de considérer à la fois les résultats des tests statistiques et les représentations graphiques pour une évaluation complète.


