# Fonction interne : retourne un message en fonction de la langue système
.msg <- function(en, fr) {
    # Obtenir la locale système pour LC_CTYPE
    sys_locale <- Sys.getlocale("LC_CTYPE")
    
    # Extraire la langue principale (avant le premier "_" ou ".")
    lang <- sub("_.*", "", sys_locale)
    lang <- sub("\\..*", "", lang)

    # Retourner le message en fonction de la langue détectée
    if (grepl("French", lang)) {
        return(fr)
    } else {
        return(en)
    }
}

# Message anglais absent : NULL, vide ou NA. On affiche alors le francais.
# (L'ancien test length(ang) == 1 remplacait TOUT message anglais d'une chaine.)
.ang_missing <- function(ang) {
  is.null(ang) || length(ang) == 0 || all(is.na(ang)) || !any(nzchar(ang))
}
