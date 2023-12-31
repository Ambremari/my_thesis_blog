---
title: "Analyse spectrale"
author: "Mathilde Couteyen Carpaye"
date: 2022-03-01
categories: ["Notes"]
tags: ["Time serie", "Spectral analysis"]
math: true
---



# Exploration des séries temporelles

Nous disposons des séries temporelles d'abondance (en cellules/mL)  pour plusieurs groupes d'organismes de 2012 à 2021. Ces séries sont celles obtenues après un lissage par régression polynomiale locale sur les données échantillonnées en moyenne deux fois par mois. 

La Figure 1 présente la série temporelle pour de l'abondance de Synechococcus à Banyuls. 

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-1-1.png" alt="Série temporelle de l'abondance de Synechococcus à Banyuls de 2012 à 2021." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-1"></span>Figure 1: Série temporelle de l'abondance de Synechococcus à Banyuls de 2012 à 2021.</p>
</div>

La série semble dirigée par des composants périodiques. Ce pourquoi on choisit de réaliser une analyse dans le domaine de fréquence. Ainsi, on étudie  la série comme variation périodiques des phénomènes sous-jacents. Cette approche apporte l'avantage de pouvoir évaluer un comportement périodique commun à deux séries. 

## Périodogramme

Pour tout échantillon de séries temporelles $\(x_1, ...x_n\)$ avec $\(n\)$ impaire on peut écrire : 

$x_t=a_0 + \sum_{j=1}^{(n-1)/2}[a_jcos(2\pi t j/n)+b_jsin(2\pi t j/n)]$

Le périodogramme échelonné peut être définit de la manière suivante : 

$P(j/n)=a^2_j+b^2_j$

Le périodogramme indique quelles composantes de fréquences ont une grande magnitude et lesquelles ont une faible magnitude. De grandes valeurs de $\(P(j/n)\)$ indique donc quelles fréquences $\(w_j=j/n\)$ sont prédominantes dans la série. Les petites valeurs de $\(P(j/n)\)$ sont plus suceptibles d'être associées au bruit. 

Les valeurs de $\(a_j\)$ et $\(b_j\)$ peuvent être obtenues avec la *fast Fourrier transform* (FFT) : 

$|d(j/n)|^2=\frac{1}{n}\Big(\sum_{t=1}^nx_tcos(2\pi tj/n\Big)^2+\frac{1}{n}\Big(\sum_{t=1}^nx_tsin(2\pi tj/n\Big)^2$

On peut obtenir le périodogramme de la manière suivante : 

$P(j/n)=\frac{4}{n}|d(j/n)|^2$

La première étape consiste à s'affranchir de la tendance présente dans notre série. En effet, les tendances vont donner de l'importance à des composantes de faible fréquence qui vont obscurcir les plus grandes fréquences. En général, on réalise un centrage des données de la forme $\(x_t-\bar{x}\)$ en amont d'une analyse spectrale. 

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" alt="Périodogramme correspondant à l'abondance de Synechococcus à Banyuls de 2012 à 2021 avec échelle de fréquence annuelle." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-2"></span>Figure 2: Périodogramme correspondant à l'abondance de Synechococcus à Banyuls de 2012 à 2021 avec échelle de fréquence annuelle.</p>
</div>

## Densité spectrale 

La densité spectrale est analogue à une densité de probabilité. Elle donne la même information que la fonction d'autocovariance. Cependant au lieu d'exprimer la covariance en terme de lags, elle l'exprime en terme de cycle. 

Soit la fonction d'autocovariance d'un processus stationnaire $\(\gamma(h)\)$, la densité spectrale s'exprime comme la transformation inverse de l'autocovariance : 

$f(\omega)=\sum_{h=-\infty}^\infty \gamma(h)e^{-2\pi i\omega h} \quad -1/2\leq \omega \leq 1/2$

# Estimation spectrale non-paramétrique

La bande de fréquence $\(\mathcal B\)$ de $\(L<<n\)$ fréquences fondamentales contigues centrée en $\(\omega_j\)$, s'exprime de la manière suivante : 

$\mathcal{B}=\Big\{\omega^*:\omega_j-\frac{m}{n}\leq\omega^*\leq\omega_j+\frac{m}{n}\Big\}$

avec $\(L=2m+1\)$ impair, choisi de manière à ce que les valeurs spectrales dans l'intervalle $\(\mathcal B\)$ : 

$f(\omega_j+k/n), \quad k=-m,...,0,...m$

soit approximativemment égales à $\(f(\omega)\)$. Les valeurs du spectre à l'intérieure de cette bande doivent être rlativement constante pour que le spectre lissé soit un bon estimateur. 

Le périodogramme moyen (ou lissé) et défini comme la moyenne des valeurs du périodogramme sur la bande $\(\mathcal B\)$.

$\bar{f}(\omega)=\frac{1}{L}\sum_{k=-m}^mI(\omega_j+k/n)$

Pour de grande valeurs de $\(n\)$, les périodogrammes moyens sont distribués comme des variables aléatoires indépendantes $\(f(\omega)\chi_2^2/2\)$ pour $\(0<\omega<1/2\)$, ainsi

$\frac{2L\bar{f}(\omega)}{f(\omega)}\sim \chi_{2L}^2$

On parle alors de fenêtre pour faire référence à la largeur de l'intervalle de fréquence définit par $\(B=L/n\)$. On peut ainsi obtenir une approximation d'un intervalle de confiance à $\(100(1-\alpha)\%\)$ de la forme : 

$\frac{2L\bar{f}(\omega)}{\chi_{2L}^2(1-\alpha/2)}\leq f(\omega) \leq \frac{2L\bar{f}(\omega)}{\chi_{2L}^2(\alpha/2)}$

La FFT est redondante lorsque $\(n\)$ est un multiple de plusieurs facteurs de 2,3 ou 5. Pour profiter de cette propriété, il est commun de rajouter des 0 pour atteindre l'entier composite le plus proche $\(n'\)$. Dans ce cas, il faut remplacer $\(2L\)$ par $\(2L/n'\)$. On définit les degrés de liberté ajustés : 

$df=\frac{2Ln}{n'}$

et l'intervalle de confiance devient :

$\frac{df\bar{f}(\omega)}{\chi_{df}^2(1-\alpha/2)}\leq f(\omega) \leq \frac{df\bar{f}(\omega)}{\chi_{df}^2(\alpha/2)}$

Après des tests avec plusieurs valeurs de fenêtre, celle qui semble la plus compatbile avec le premier spectre (Figure 2) est celle obtenue avec $\(L=3\)$ avec un noyau de Daniell (Figure 3). 

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" alt="Périodogramme moyen correspondant à l'abondance de Synechococcus à Banyuls de 2012 à 2021 avec échelle de fréquence en cycle annuel. $\(n=\)$ 3243 , $\(n'=\)$ 3375 , $\(L=\)$ 3 , $\(df=\)$ 5.765" width="672" />
<p class="caption"><span id="fig:unnamed-chunk-3"></span>Figure 3: Périodogramme moyen correspondant à l'abondance de Synechococcus à Banyuls de 2012 à 2021 avec échelle de fréquence en cycle annuel. $\(n=\)$ 3243 , $\(n'=\)$ 3375 , $\(L=\)$ 3 , $\(df=\)$ 5.765</p>
</div>

Il faut noter que les fréquences qui apparaissent sur la Figure 3 correspondent à des cycles annuels. On retrouve donc bien notre cycle annuel (1/1) et notre cycle semestriel (1/2 an = 6 mois) ainsi que des pics mineurs au niveau des harmoniques. Les intrevalles de confiance associés sont présentés dans le Tableau 1. Dans les deux cas, la borne inférieure de l'intrvalle de confiance a une puissance supérieure au reste du spectre sans ces deux pics. 



\begin{table}[h]
\caption{Interval de confiance à $\(95\%\)$ pour le spectre de l'abondance de Synechococcus à Banyuls de 2012 à 2021.}
\centering
\begin{tabular}{c c c c c}
\hline
$\(\omega\)$ & Période & Puissance & Inférieur & Supérieur\\
\hline
1/365 & 1 an & 1.2\times 10^{8} & 5.1\times 10^{7} & 6.3\times 10^{8}\\
1/180 & 6 mois & 4.8\times 10^{7} & 2\times 10^{7} & 2.5\times 10^{8} \\
\hline
\end{tabular}
\end{table}

## Harmoniques

Comme c'est le cas pour le spectre de l'abondance de Synechococcus, il est commun d'observer des pics mineurs au niveau des harmoniques. Ici, on a un pic important à $\(\omega = 1\Delta\)$ un cycle d'un an, un pic moins important à $\(\omega=2\Delta\)$ et des pics mineurs aux autres harmoniques $\(\omega=k\Delta\)$ pour $\(k=3,4\)$ ... Ce cas est courant car la plupart des signaux ne sont pas des sinusoïdes parfait. Les harmoniques permettent alors de capturer le comportement non-sinyusoïdal du signal. La question est de décider si un pic est significatif ou non. On s'intéresse alors aux affirmations simultanées sur une gamme de fréquences. 

Soit $\(K\)$ affirmations $\(S_1, S_2,..., S_k\)$ réalisées à un niveau de confiance $\(\alpha\)$, $\(P\{S_k\}=1-\alpha\)$ alors la probabilité que toutes les affirmations soient vraies satisfait l'inéquation de Bonferroni : 

$P\{\text{tous }S_K \text{ vrais}\}\geq 1-K\alpha$

Il est préferrable de fixer le niveau de confiance pour les tests de chaque fréquence à $\(\alpha/K\)$ s'il y a $\(K\)$ potentielles fréquences d'intérêt. 

Il est possible de renforcer l'estimation en utilisant une moyenne pondérée : 

$\hat f (\omega) = \sum_{k=-m}^mh_kI(\omega_j+k/n) $

avec les poids $\(h_k>0\)$ qui satisfont : 

$\sum_{k=-m}^mh_k=1$

On fixe $\(h_k=L^{-1}\)$ pour tous $\(k\)$, dans ce cas $\(\sum_{k=-m}^mh_k^2=L^{-1}\)$. On définit ainsi :

$L_h=\bigg(\sum_{k=-m}^mh_k^2\bigg)^{-1}$

et réalise l'approximation : 

$\frac{2L_h\hat f(\omega)}{f(\omega)}\sim\chi_{2L_h}^2$

La fenêtre est alors $\(B=L_h/n\)$ et l'intervalle de confiance est de la forme : 

$\frac{2L_h\hat{f}(\omega)}{\chi_{2L_h}^2(1-\alpha/2)}\leq f(\omega) \leq \frac{2L_h\hat{f}(\omega)}{\chi_{2L_h}^2(\alpha/2)}$

On choisis d'appliquer un noyau de Daniell modifié à $\(m=1\)$ deux fois. On obtient ainsi $\(L_h=$r Lh\)$ proche de notre $\(L=3\)$ précédent. 




## Tapering

Lors de la création de fenêtre, les fenêtre sans *taper* peuvent produire des ondulations sur la bande voulu et à l'extérieur de la bande. Ces dernières sont appelées *sidelobes*. Ainsi, il est possible que certaines fréquences en dehors de l'intervalle voulu contaminent la bande. Pour supprimer cet effet de fuite, il est possible d'utiliser des *taper*. 

Les *taper* ont généralement une forme qui amplifie le centre des données par rapport aux extrémités, comme une cloche cosinus de la forme : 

$h_t=0.5\bigg[1 + cos\bigg(\frac{2\pi(t- \bar t)}{n}\bigg)\bigg]$

On peut remplacer une série orginiale par une sérier *tapered* : 

$y_t = h_tx_t$

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" alt="Périodogramme moyen de l'abondance de Synechococcus à banyuls de 2012 à 2021 avec noyau de Daniell modifié sans taper (ligne pointillée) et avec tapering complet par une cloche cosinus (ligne pleine)." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-6"></span>Figure 4: Périodogramme moyen de l'abondance de Synechococcus à banyuls de 2012 à 2021 avec noyau de Daniell modifié sans taper (ligne pointillée) et avec tapering complet par une cloche cosinus (ligne pleine).</p>
</div>


# Filtres linéaires

Les filtres linéaires peuvent être utilisés pour extraire des signaux des séries temporelles. 

Si notre série filtrée $\(y_t\)$ peut s'écrire de la manière suivante :

$y_t=\sum_{j=-\infty}^\infty a_jx_{t-j} \quad \sum_{j=-\infty}^\infty|a_j|<\infty$

alors $\(x_t\)$ a un spectre $\(f_{xx}(\omega)\)$ et $\(y_t\)$ a un spectre 

$f_{yy}(\omega)=|A_{yx}(\omega)|^2f_{xx}(\omega)$

Avec $\(A_{yx}(\omega)\)$ la fonction de fréquence de réponse : 

$A_{yx}(\omega)=\sum_{j=-\infty}^\infty a_je^{-2\pi i\omega j}$

Autrement dit, l'effet du filtrage peut être caractérisé par le produit de la fréquence et de la magnitude de la fonction de fréquence réponse au carré. 

On présente deux filtres, la première différence : 

$y_t=\nabla x_t= x_t-x_{t-1}$

et le filtre par moyenne annuelle symétrique mobile, avec $\(t\)$ exprimé en mois ou noyau de Daniell modifié avec $\(m=6\)$ : 

$y_t=\frac{1}{24}(x_{t-6}+x_{t+6}) + \frac{1}{12}\sum_{r=-5}^5x_{t-r}$

Comme on peut le voir sur la Figure 3, l'effet de la différenciation est de *roughen* la série en retenant les plus grandes fréquences, c'est un exemple de filtre *high-pass*. A l'inverse, l'effet de la moyenne mobile centrée est de lisser la série en gardant les basses fréquences et en atténuant les hautes fréquences. C'est un exemple de filtre *low-pass*. 

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" alt="Comparaison de la série temporelle d'abondance de Synechococcus à banyuls (haut), avec la série différenciée (milieu) et la moyenne mobile centrée à 12 mois." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-7"></span>Figure 5: Comparaison de la série temporelle d'abondance de Synechococcus à banyuls (haut), avec la série différenciée (milieu) et la moyenne mobile centrée à 12 mois.</p>
</div>

Pour comprendre de quelle manière les filtres affectent le signal original, on étudie le carré de la fréquence de réponse $\(|A_{yx}(\omega)|^2\)$ (Figure 4). Pour la différenciation on a :

$A_{yx}(\omega)=1-e^{-2\pi i \omega}$

Pour moyenne mobile centrée à 12 mois on a : 

$A_{yx}(\omega)=\frac{1}{12}\Big[1+cos(12\pi\omega) + 2\sum_{k=1}^5cos(2\pi\omega k)\Big]$

Le filtrage par différenciation va renforcer les hautes fréquences notamment à partir 3 ce qui correspond à une période de 120 jours (Figure 6). A l'inverse, le filtrage par moyenne mobile va fortement atténuer l'effet des fréquences supérieures à 0.8. Ainsi, on voit disparaître les cycles annuels et semestriels. 

Dans le deux cas on a éliminé le cycle annuel. Dans le premier cas les variations à petites échelles sont mises en avant. Dans le second cas on a une mise en évidence de la tendance sur la période globale de la série. 

On va tenter d'améliorre le flitrage du signal pour ne garder que les fréquences inférieures à 1 cycle annuel et ne pas introduire des fréquences supérieures à 1 comme le réalise légèrement le filtre par moyenne mobile. 


<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" alt="Carré de la fréquence de réponse pour les filtres par différentation (haut) et par moyenne mobile sur 12 mois (bas)." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-8"></span>Figure 6: Carré de la fréquence de réponse pour les filtres par différentation (haut) et par moyenne mobile sur 12 mois (bas).</p>
</div>


# Extraction du signal et filtrage optimal

Dans le cas d'un modèle de signal avec bruit : 

$y_t=x_t+v_t$

on s'intéresse à l'estimation d'un signal $\(x_t\)$ de la forme : 

$\hat x_t=\sum_{r=-\infty}^\infty a_ry_{t-r}$

On cherche les coefficients du filtre $\(a_t\)$ qui minimisent le moyenne de l'erreur quadratique : 

$MSE=E\bigg[\bigg(x_t-\sum_{r=-\infty}^\infty a_ry_{t-r}\bigg)^2\bigg]$

Cela revient à résoudre 

$\sum_{r=-\infty}^\infty a_r\gamma_{yy}(s-r)=\gamma_{xy}(s)$

qui équivaut à :

$A(\omega)f_{yy}(\omega)=f_{xy}(\omega)$

$\(A(\omega)\)$ et le filtre optimal $\(a_t\)$ sont des transfomrations paires de Fourrier pour $\(B(\omega)\)$ et $\(\beta_t\)$. 

$f_{xy}(\omega)=B^*(\omega)f_{xx}(\omega)$

$f_{yy}(\omega)=|B(\omega)|^2f_{xx}(\omega)+f_{vv}(\omega)$

Impliquent que le filtre optimal serait la transformation de Fourrier de : 

$A(\omega)=\frac{B^*(\omega)}{\Big(|B(\omega)|^2+\frac{f_{vv}(\omega)}{f_{xx}(\omega)}\Big)}$

où le second terme du dénominateur est l'inverse du ratio signal-bruit : 

$SNR(\omega)=\frac{f_{xx}(\omega)}{f_{vv}(\omega)}$

Les filtres optimaux peuvent donc être calculés si les spectres du signal et du bruit sont connus ou bien si on suppose connu le ratio signal-bruit. Si le ratio signal-bruit est connu, le filtre optimal peut être calculé comme la transformation inverse de la fonction $\(A(\omega)\)$. On aurait alors la fonction d'estimation du filtre : 

$a_t^M=M^{-1}\sum_{k=0}^{M-1}A(\omega_k)e^{2\pi i\omega_kt}$

Avec $\(M<<n\)$ et en pratique pour $\(|t|>M/2\)$, $\(a_t=0\)$. 

Il est possible que la forme de la fréquence réponse affiche des transitions abruptes entre des régions où le ratio signal-bruit est fort et où le signal est faible. La fréquence réponse pourrait alors présenter des ondulations qui introduieraient des fréquences à différentes amplitudes. Une solution pour palier à ce problème est d'introduir le *tapering*. On utilise le filtre *tapered* $\(\tilde a=h_ta_t\)$ où $\(h_t\)$ est le *taper* cosinus. La fonction réponse au carré est alors $\(|\tilde A(\omega)|^2\)$ où : 

$\tilde A(\omega)=\sum_{t=-\infty}^\infty a_th_te^{-2\pi i\omega t}$

On garde un $\(L=3\)$ comme pour les analyses précédentes. On ajuste le choix de $\(M\)$ en fonction de l'allure de la fonction de fréquence réponse atteinte (Figure 7). Le $\(M\)$ chosi est 1000 (jours). 


<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" alt="Coefficients du filtre avec cosinus taper pour l'abondance de Synechococcus à Banyuls (haut). Fonctions de fréquence réponse désirée (en bleu, pointillé) et atteinte (en rouge) pour l'abondance de Synechococcus à Banyuls avec fréquence en cycle annuel (bas). $\(L=3\)$ et $\(M=1000\)$" width="672" />
<p class="caption"><span id="fig:unnamed-chunk-9"></span>Figure 7: Coefficients du filtre avec cosinus taper pour l'abondance de Synechococcus à Banyuls (haut). Fonctions de fréquence réponse désirée (en bleu, pointillé) et atteinte (en rouge) pour l'abondance de Synechococcus à Banyuls avec fréquence en cycle annuel (bas). $\(L=3\)$ et $\(M=1000\)$</p>
</div>

On voit que le filtre accorde encore une part d'importance au cycle annuel. Néanmoins il va complètement filtrer toutes les fréquences supérieures à 1,5. Le signal filtré obtenu est illsutré en Figure 8.


<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" alt="Série originale de l'abondance de Synechococcus à Banyuls (haut) comparée à la série filtrée (bas)." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-10"></span>Figure 8: Série originale de l'abondance de Synechococcus à Banyuls (haut) comparée à la série filtrée (bas).</p>
</div>

Le problème qui se présente de manière évidente et que le choix d'un $\(M=1000\)$ nous prive de 500 valeurs en début et 500 valeurs en fin de série. Néanmoins, ce filtre nous permet d'avoir une idée de la tendance générale en identifiant clairement les années et les cycles annuels. 

# Séries multiples et spectres croisés 

Cohérence...
