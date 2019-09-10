library(tidyverse)

# fonction de projection de population utilisant un modèle par cohortes et composantes (j'ai écrit ça rapidement
# et ce n'est sans doute pas la façon la plus élégante de faire ça, mais ça ira pour les besoins de ce billet)
projection_population <- function(population,
                                  fécondité,
                                  mortalité,
                                  solde_migratoire,
                                  ratio_sexe_naissance,
                                  période) {
  # nombre de classes d'âge dans les vecteurs de population fournis en argument
  nb_classes_âge <- length(population$hommes_européens)
  
  # vecteur de population au format matriciel pour les hommes d'origine étrangère européenne et natifs
  population_hommes_européens <- matrix(population$hommes_européens,
                                        nrow = nb_classes_âge,
                                        ncol = 1)
  # vecteur de population au format matriciel pour les femmes d'origine étrangère européenne et natives
  population_femmes_européennes <- matrix(population$femmes_européennes,
                                          nrow = nb_classes_âge,
                                          ncol = 1)
  
  # vecteur de population au format matriciel pour les hommes immigrés non-européens
  population_hommes_immigrés_non_européens <- matrix(population$hommes_immigrés_non_européens,
                                                     nrow = nb_classes_âge,
                                                     ncol = 1)
  
  # vecteur de population au format matriciel pour les femmes immigrées non-européennes
  population_femmes_immigrées_non_européennes <- matrix(population$femmes_immigrées_non_européennes,
                                                        nrow = nb_classes_âge,
                                                        ncol = 1)
  
  # vecteur de population au format matriciel pour les hommes natifs d'origine non-européenne
  population_hommes_natifs_non_européens <- matrix(population$hommes_natifs_non_européens,
                                                   nrow = nb_classes_âge,
                                                   ncol = 1)
  
  # vecteur de population au format matriciel pour les femmes natives d'origine non-européenne
  population_femmes_natives_non_européennes <- matrix(population$femmes_natives_non_européennes,
                                                      nrow = nb_classes_âge,
                                                      ncol = 1)
  
  # vecteur de population pour l'ensemble des groupes
  population_totale <- population_hommes_européens +
    population_femmes_européennes +
    population_hommes_immigrés_non_européens +
    population_femmes_immigrées_non_européennes +
    population_hommes_natifs_non_européens +
    population_femmes_natives_non_européennes
  
  # matrice de taux de survie pour les hommes
  taux_survie_hommes <- cbind(diag(1 - mortalité$hommes), rep(0, nb_classes_âge - 1))
  
  # matrice de taux de survie pour les femmes
  taux_survie_femmes <- cbind(diag(1 - mortalité$femmes), rep(0, nb_classes_âge - 1))
  
  # matrice de Leslie pour les hommes européens et/ou natifs
  leslie_hommes_européens <- rbind(matrix(rep(0, nb_classes_âge),
                                          nrow = 1,
                                          ncol = nb_classes_âge),
                                   taux_survie_hommes)
  
  # matrice de Leslie pour les femmes européennes et/ou natives
  leslie_femmes_européennes <- rbind(matrix(fécondité$européennes,
                                            nrow = 1,
                                            ncol = nb_classes_âge),
                                     taux_survie_femmes)
  
  # matrice de Leslie pour les hommes immigrés non-européens
  leslie_hommes_immigrés_non_européens <- rbind(matrix(rep(0, nb_classes_âge),
                                                       nrow = 1,
                                                       ncol = nb_classes_âge),
                                                taux_survie_hommes)
  
  # matrice de Leslie pour les femmes immigrées non-européennes
  leslie_femmes_immigrées_non_européennes <- rbind(matrix(fécondité$immigrées_non_européennes,
                                                          nrow = 1,
                                                          ncol = nb_classes_âge),
                                                   taux_survie_femmes)
  
  # matrice de Leslie pour les hommes natifs d'origine non-européene
  leslie_hommes_natifs_non_européens <- rbind(matrix(rep(0, nb_classes_âge),
                                                     nrow = 1,
                                                     ncol = nb_classes_âge),
                                              taux_survie_hommes)
  
  # matrice de Leslie pour les femmes natives d'origine non-européenne
  leslie_femmes_natives_non_européennes <- rbind(matrix(fécondité$natives_non_européennes,
                                                        nrow = 1,
                                                        ncol = nb_classes_âge),
                                                 taux_survie_femmes)
  
  # vecteur au format matriciel qui servira à calculer à chaque itération le solde migratoire par
  # âge des hommes non-européens à partir du vecteur de population pour l'ensemble des groupes
  solde_hommes <- matrix(solde_migratoire$total / sum(population_totale) * solde_migratoire$distribution_âge_hommes,
                         nrow = nb_classes_âge,
                         ncol = 1)
  
  # vecteur au format matriciel qui servira à calculer à chaque itération le solde migratoire par
  # âge des femmes non-européennes à partir du vecteur de population pour l'ensemble des groupes
  solde_femmes <- matrix(solde_migratoire$total / sum(population_totale) * solde_migratoire$distribution_âge_femmes,
                         nrow = nb_classes_âge,
                         ncol = 1)
  
  # crée la structure dans laquelle seront enregistrés les résultats
  résultats <- tibble(année = période,
                      proportion_immigrés_non_européens = numeric(length(période)),
                      proportion_natifs_non_européens = numeric(length(période)),
                      proportion_total_non_européens = numeric(length(période)),
                      proportion_nouveau_nés_non_européens = numeric(length(période)))
  
  # fait tourner le modèle
  for (i in seq_along(période)) {
    # enregistre les résultats de l'itération précédente ou, dans le cas où i = 1,
    # les conditions initiales dans la structure prévue à cet effet
    résultats$proportion_immigrés_non_européens[i] <- sum(population_hommes_immigrés_non_européens + population_femmes_immigrées_non_européennes) / sum(population_totale)
    résultats$proportion_natifs_non_européens[i] <- sum(population_hommes_natifs_non_européens + population_femmes_natives_non_européennes) / sum(population_totale)
    résultats$proportion_total_non_européens[i] <- sum(population_hommes_immigrés_non_européens + population_femmes_immigrées_non_européennes + population_hommes_natifs_non_européens + population_femmes_natives_non_européennes)/sum(population_totale)
    résultats$proportion_nouveau_nés_non_européens[i] <- (population_hommes_immigrés_non_européens[1, 1] + population_femmes_immigrées_non_européennes[1, 1] + population_hommes_natifs_non_européens[1, 1] + population_femmes_natives_non_européennes[1, 1]) / population_totale[1, 1]
    
    # calcule la population immigrée par âge et sexe
    population_hommes_immigrés_non_européens <- round(leslie_hommes_immigrés_non_européens %*% population_hommes_immigrés_non_européens +
      sum(population_totale) * solde_hommes)
    population_femmes_immigrées_non_européennes <- round(leslie_femmes_immigrées_non_européennes %*% population_femmes_immigrées_non_européennes +
      sum(population_totale) * solde_femmes)
    
    # calcule la population européenne par âge et sexe
    population_hommes_européens <- round(leslie_hommes_européens %*% population_hommes_européens)
    population_femmes_européennes <- round(leslie_femmes_européennes %*% population_femmes_européennes)
    population_hommes_européens[1, 1] <- round(population_femmes_européennes[1, 1] * ratio_sexe_naissance / (1 + ratio_sexe_naissance))
    population_femmes_européennes[1, 1] <- round(population_femmes_européennes[1, 1] * (1 - ratio_sexe_naissance / (1 + ratio_sexe_naissance)))
    
    # calcule la population des natifs d'origine non-européenne par âge et sexe
    population_hommes_natifs_non_européens <- round(leslie_hommes_natifs_non_européens %*% population_hommes_natifs_non_européens)
    population_femmes_natives_non_européennes <- round(leslie_femmes_natives_non_européennes %*% population_femmes_natives_non_européennes)
    population_hommes_natifs_non_européens[1, 1] <- round(population_femmes_natives_non_européennes[1, 1] * ratio_sexe_naissance / (1 + ratio_sexe_naissance) +
      population_femmes_immigrées_non_européennes[1, 1] * ratio_sexe_naissance / (1 + ratio_sexe_naissance))
    population_femmes_natives_non_européennes[1, 1] <- round(population_femmes_natives_non_européennes[1, 1] * (1 - ratio_sexe_naissance / (1 + ratio_sexe_naissance)) +
      population_femmes_immigrées_non_européennes[1, 1] * (1 - ratio_sexe_naissance / (1 + ratio_sexe_naissance)))
    
    # enlève de la population immigrée tous les individus nés en France ajoutés plus haut
    # pour ne garder que ceux issus de l'immigration
    population_hommes_immigrés_non_européens[1, 1] <- round((population_totale * solde_hommes)[1, 1])
    population_femmes_immigrées_non_européennes[1, 1] <- round((population_totale * solde_femmes)[1, 1])
    
    # calcule la population totale en additionnant les différentes composantes
    population_totale <- population_hommes_européens +
      population_femmes_européennes +
      population_hommes_immigrés_non_européens +
      population_femmes_immigrées_non_européennes +
      population_hommes_natifs_non_européens +
      population_femmes_natives_non_européennes
  }
  
  return (résultats)
}

# source : https://www.insee.fr/fr/statistiques/1892259?sommaire=1912926#titre-bloc-5 (fécondité par âge en 2015)
fécondité_insee <- read_csv2("fécondité2015.csv") %>%
  mutate(TAUX = TAUX / 10000)

# crée les vecteurs de taux de fécondité par âge pour les Européens nés en France et les immigrés
# européens, les individus d'origine non-européenne nés en France et les immigrés non-européens
# note : 1) je fais l'hypothèse que les descendants d'immigrés non-européens ont les mêmes taux de
# fécondité par âge que les individus d'origine européenne et que ces derniers ont également tous les
# mêmes taux qu'ils soient immigrés ou descendants d'immigrés et 2) les taux sont ajustés de manière
# approximative (cette méthode suppose notamment que la seule différence entre les Européennes et les
# immigrés non-européennes est le nombre d'enfants qu'elles ont et pas l'âge auquel elles les ont, ce
# qui est sans doute faux, mais ne change rien à long terme) sur la base des chiffres de l'INSEE qui
# sont publiés sur cette page web : https://www.insee.fr/fr/statistiques/3675496
fécondité <- tibble(européennes = c(rep(0, 14), fécondité_insee$TAUX, rep(0, 51)) * 1.8 / 1.92,
                    natives_non_européennes = c(rep(0, 14), fécondité_insee$TAUX, rep(0, 51)) * 1.8 / 1.92,
                    immigrées_non_européennes = c(rep(0, 14), fécondité_insee$TAUX, rep(0, 51)) * 3 / 1.92)

# source : https://www.insee.fr/fr/statistiques/3311422?sommaire=3311425 (mortalité par âge en 2016)
mortalité_insee <- read_csv2("mortalité2016.csv") %>%
  filter(AGE < 100) %>%
  mutate(TAUX = TAUX/100000)

# crée les vecteurs de mortalité
# note : je fais l'hypothèse dans le modèle que tous les groupes d'origine ont les mêmes taux de mortalité,
# ce qui est faux à n'en pas douter, mais je doute que ça soit de nature à faire une grosse différence
mortalité <- tibble(hommes = mortalité_insee$TAUX[mortalité_insee$SEXE == 1],
                    femmes = mortalité_insee$TAUX[mortalité_insee$SEXE == 2])

# source : https://www.insee.fr/fr/statistiques/3565914?sommaire=3558417 (recensement de 2015 localisé à la
# région, fichier détail don't j'ai conservé seulement les variables qui m'intéressaient)
recensement_insee <- readRDS("recensement2015.rds")

# crée une variable indiquant si un individu est un immigré non-européen
recensement_insee <- recensement_insee %>%
  mutate(NON_EURO = ifelse(PNAI12 > 6 & INAT != 11, TRUE, FALSE))

# crée une structure avec la population par âge, sexe et origine
population_insee <- recensement_insee %>%
  group_by(AGEREV, SEXE, NON_EURO) %>%
  summarize(N = round(sum(IPONDI)))

# exclut tous les individus de plus de 100 ans, âge au-delà duquel on considère dans
# la projection que tous les gens sont morts
population_insee <- population_insee %>%
  filter(AGEREV <= 100)

# crée les vecteurs de population
population <- tibble(hommes_européens = population_insee$N[population_insee$SEXE == 1 & population_insee$NON_EURO == FALSE],
                     femmes_européennes = population_insee$N[population_insee$SEXE == 2 & population_insee$NON_EURO == FALSE],
                     hommes_immigrés_non_européens = population_insee$N[population_insee$SEXE == 1 & population_insee$NON_EURO == TRUE],
                     femmes_immigrées_non_européennes = population_insee$N[population_insee$SEXE == 2 & population_insee$NON_EURO == TRUE],
                     hommes_natifs_non_européens = rep(0, 101),
                     femmes_natives_non_européennes = rep(0, 101))

# calcule le nombre d'immigrés arrivés depuis moins de 5 ans par âge
distribution_sexe_âge_solde <- recensement_insee %>%
  filter(ANARR == 0 & NON_EURO == TRUE & AGEREV <= 100) %>%
  group_by(AGEREV, SEXE) %>%
  summarize(N = round(sum(IPONDI))) %>%
  ungroup() %>%
  complete(AGEREV, SEXE, fill = list(N = 0)) %>%
  mutate(proportion = N / sum(N))

# crée la structure avec les informations sur le solde migratoire pour les immigrés
# non-européens, sur la base du chiffre donné par Hervé Le Bras au Nouvel Observateur
# source : https://www.nouvelobs.com/societe/20180109.OBS0357/immigration-que-repondre-a-votre-beau-frere-qui-croit-au-grand-remplacement.html
solde_migratoire_non_européen <- list(total = 130000,
                                      distribution_âge_hommes = distribution_sexe_âge_solde$proportion[distribution_sexe_âge_solde$SEXE == 1],
                                      distribution_âge_femmes = distribution_sexe_âge_solde$proportion[distribution_sexe_âge_solde$SEXE == 2])

# fait tourner le modèle avec le scénario dans lequel le solde migratoire reste au
# niveau actuel en proportion de la population totale durant toute la période
résultats1 <- projection_population(population,
                                   fécondité,
                                   mortalité,
                                   solde_migratoire_non_européen,
                                   1.05,
                                   2015:2100)

# prépare la structure qui contient les résultats pour ggplot
résultats1 <- résultats1 %>%
  gather(key = groupe,
         value = proportion,
         proportion_immigrés_non_européens,
         proportion_natifs_non_européens,
         proportion_total_non_européens,
         proportion_nouveau_nés_non_européens) %>%
  filter(groupe != "proportion_nouveau_nés_non_européens")

# crée un graphique qui montre l'évolution de la part des différentes composantes d'origine
# non-européenne de la population entre 2015 et 2100 d'après la projection
ggplot(résultats1, aes(x = année, y = proportion, group = groupe, color = groupe)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100\n(scénario avec immigration au niveau actuel)") +
  xlab("Année") +
  ylab("Proportion") +
  scale_color_discrete(name = "Groupe", labels = c("Immigrés d'origine\nnon-européenne", "Natifs d'origine\nnon-européenne", "Ensemble des individus\nd'origine non-européenne")) +
  scale_x_continuous(breaks = seq(2015, 2100, 5)) +
  scale_y_continuous(breaks = seq(0, plyr::round_any(max(résultats1$proportion), 0.05), 0.05), labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100 (scénario avec immigration au niveau actuel).png", width = 12, height = 6)

# fait tourner le modèle avec le scénario dans lequel le solde migratoire double par rapport
# au niveau actuel en proportion de la population totale durant toute la période
solde_migratoire_non_européen$total <- 260000
résultats2 <- projection_population(population,
                                   fécondité,
                                   mortalité,
                                   solde_migratoire_non_européen,
                                   1.05,
                                   2015:2100)

# prépare la structure qui contient les résultats pour ggplot
résultats2 <- résultats2 %>%
  gather(key = groupe,
         value = proportion,
         proportion_immigrés_non_européens,
         proportion_natifs_non_européens,
         proportion_total_non_européens,
         proportion_nouveau_nés_non_européens) %>%
  filter(groupe != "proportion_nouveau_nés_non_européens")

# crée un graphique qui montre l'évolution de la part des différentes composantes d'origine
# non-européenne de la population entre 2015 et 2100 d'après la projection
ggplot(résultats2, aes(x = année, y = proportion, group = groupe, color = groupe)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100\n(scénario avec une immigration 2 fois plus importante qu'actuellement)") +
  xlab("Année") +
  ylab("Proportion") +
  scale_color_discrete(name = "Groupe", labels = c("Immigrés d'origine\nnon-européenne", "Natifs d'origine\nnon-européenne", "Ensemble des individus\nd'origine non-européenne")) +
  scale_x_continuous(breaks = seq(2015, 2100, 5)) +
  scale_y_continuous(breaks = seq(0, plyr::round_any(max(résultats2$proportion), 0.05), 0.05), labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100 (scénario avec une immigration 2 fois plus importante qu'actuellement).png", width = 12, height = 6)

# fait tourner le modèle avec le scénario dans lequel le solde migratoire reste au niveau actuel en proportion
# de la population totale durant toute la période, mais la fécondité des femmes immigrées non-européennes est
# identique à celle des autres femmes
solde_migratoire_non_européen$total <- 130000
fécondité$immigrées_non_européennes <- c(rep(0, 14), fécondité_insee$TAUX, rep(0, 51)) * 1.8 / 1.92
résultats3 <- projection_population(population,
                                   fécondité,
                                   mortalité,
                                   solde_migratoire_non_européen,
                                   1.05,
                                   2015:2100)

# prépare la structure qui contient les résultats pour ggplot
résultats3 <- résultats3 %>%
  gather(key = groupe,
         value = proportion,
         proportion_immigrés_non_européens,
         proportion_natifs_non_européens,
         proportion_total_non_européens,
         proportion_nouveau_nés_non_européens) %>%
  filter(groupe != "proportion_nouveau_nés_non_européens")

# crée un graphique qui montre l'évolution de la part des différentes composantes d'origine
# non-européenne de la population entre 2015 et 2100 d'après la projection
ggplot(résultats3, aes(x = année, y = proportion, group = groupe, color = groupe)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100\n(scénario avec immigration au niveau actuel et fécondité des immigrées non-européennes identique à celle des autres femmes)") +
  xlab("Année") +
  ylab("Proportion") +
  scale_color_discrete(name = "Groupe", labels = c("Immigrés d'origine\nnon-européenne", "Natifs d'origine\nnon-européenne", "Ensemble des individus\nd'origine non-européenne")) +
  scale_x_continuous(breaks = seq(2015, 2100, 5)) +
  scale_y_continuous(breaks = seq(0, max(résultats3$proportion), 0.05), labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100 (scénario avec immigration au niveau actuel et fécondité des immigrées non-européennes identique à celle des autres femmes).png", width = 12, height = 6)

population$hommes_natifs_non_européens[0:30] <- population$hommes_européens[0:30] * 0.10
population$femmes_natives_non_européennes[0:30] <- population$femmes_européennes[0:30] * 0.10

population$hommes_européens <- population$hommes_européens - population$hommes_natifs_non_européens
population$femmes_européennes <- population$femmes_européennes - population$femmes_natives_non_européennes

fécondité$immigrées_non_européennes <- c(rep(0, 14), fécondité_insee$TAUX, rep(0, 51)) * 3 / 1.92

résultats4 <- projection_population(population,
                                    fécondité,
                                    mortalité,
                                    solde_migratoire_non_européen,
                                    1.05,
                                    2015:2100)

# prépare la structure qui contient les résultats pour ggplot
résultats4 <- résultats4 %>%
  gather(key = groupe,
         value = proportion,
         proportion_immigrés_non_européens,
         proportion_natifs_non_européens,
         proportion_total_non_européens,
         proportion_nouveau_nés_non_européens) %>%
  filter(groupe != "proportion_nouveau_nés_non_européens")

# crée un graphique qui montre l'évolution de la part des différentes composantes d'origine
# non-européenne de la population entre 2015 et 2100 d'après la projection
ggplot(résultats4, aes(x = année, y = proportion, group = groupe, color = groupe)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100\n(scénario avec immigration au niveau actuel et présence de descendants d'immigrés non-européens)") +
  xlab("Année") +
  ylab("Proportion") +
  scale_color_discrete(name = "Groupe", labels = c("Immigrés d'origine\nnon-européenne", "Natifs d'origine\nnon-européenne", "Ensemble des individus\nd'origine non-européenne")) +
  scale_x_continuous(breaks = seq(2015, 2100, 5)) +
  scale_y_continuous(breaks = seq(0, max(résultats4$proportion), 0.05), labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave("Projection de la part des individus d'origine non-européenne dans la population entre 2015 et 2100 (scénario avec immigration au niveau actuel et présence de descendants d'immigrés non-européens).png", width = 12, height = 6)
