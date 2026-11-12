# Turns every non-ASCII character of an R file into a \uXXXX escape (they all sit in strings).
f <- commandArgs(TRUE)
for (p in f) {
	x <- readLines(p, encoding = "UTF-8", warn = FALSE)
	y <- vapply(x, function(s) {
		cp <- utf8ToInt(enc2utf8(s))
		paste(ifelse(cp > 127, paste0(intToUtf8(92L), "u", sprintf("%04x", cp)), intToUtf8(cp, multiple = TRUE)), collapse = "")
	}, character(1), USE.NAMES = FALSE)
	writeLines(y, p, useBytes = TRUE)
	cat(p, sum(x != y), "lines escaped\n")
}
