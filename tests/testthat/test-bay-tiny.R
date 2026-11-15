# m.bay : groupe d'une seule valeur (I3, passe 3 de l'expert, D32), peu de valeurs (3 + 3),
# seuils du soutien alignes sur l'echelle des facteurs de Bayes imprimee dans le bilan.

test_that("un groupe d'une valeur : pas d'erreur brute, lois a variance commune, soutien plafonne", {
	set.seed(32)
	d <- data.frame(y = c(rnorm(20, 10, 1), 16), g = factor(c(rep("a", 20), "b")))
	msg <- capture.output(out <- capture.output(r <- m.bay(y ~ g, data = d, draws = 1000)), type = "message")
	all_txt <- paste(c(out, msg), collapse = "\n")
	expect_false(grepl("not enough observations", all_txt))
	expect_false(grepl("Error in", all_txt))
	expect_true(isTRUE(r$tiny_group))
	expect_true(grepl("une seule valeur|one value", all_txt))
	expect_true(all(r$laws$law %in% c("normal_eq", "t_eq")))
	expect_false(grepl("soutien fort|soutien très fort|strong support", all_txt))
})

test_that("3 + 3 valeurs : la forme du bruit n'est pas annoncee quand aucune ne se detache", {
	set.seed(33)
	d <- data.frame(y = c(rnorm(3, 10), rnorm(3, 11)), g = factor(rep(c("a", "b"), each = 3)))
	out <- gsub("[[:space:]]+", " ", paste(capture.output(r <- m.bay(y ~ g, data = d, draws = 1000)), collapse = " "))
	if (max(r$laws$posterior) < 0.5) {
		expect_true(any(grepl("trop peu de valeurs pour reconnaître la forme du bruit|too few values to recognise the shape of the noise", out)))
	} else {
		expect_false(any(grepl("trop peu de valeurs pour reconnaître|too few values to recognise", out)))
	}
})

test_that("soutien : 20 et 150 bornent le soutien fort, comme l'echelle imprimee", {
	expect_equal(.bay_support(19 / 20)$level, 3L)
	expect_equal(.bay_support(0.96)$level, 4L)
	expect_equal(.bay_support(149 / 150)$level, 4L)
	expect_equal(.bay_support(0.995)$level, 5L)
})
