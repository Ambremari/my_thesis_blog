---
title: "Banyuls"
author: Mathilde COUTEYEN CARPAYE
date: 2022-02-20
categories: ["Visualisation"]
tags: ["Visualisation", "Banyuls"]
math: true
description : Visualisation des données de température de la station de Banyuls.
---



# Profils verticaux de Température

On dispose des profils de Température de 0 à 50m de profondeur à Banyuls de 2019 à 2021. Les observations sont réalisées environ tous les deux mois. La sonde CTD permet un relevé plus ou moins tous les 0.10-0.20m. 

Pour obtenir des observations aux mêmes profondeurs, on réalise une interpolation sur les données existantes pour obtenir des valeurs tous les 0.5m de 0.5 à 50m de profondeur. 



On propose de réaliser une ACP fonctionnelle pour les profils de température pour chaque mois. On choisi de faire un liassage avec une spline cubique avec des noeuds localisés aux observations. On fixe notre paramètre de lissage $\(\lambda=10e-3\)$.   



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" alt="a. Cinq premières composantes principales de l'ACP fonctionnelle réalisée sur les profils de température en janvier de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-2"></span>Figure 1: a. Cinq premières composantes principales de l'ACP fonctionnelle réalisée sur les profils de température en janvier de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en janvier, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-3"></span>Figure 2: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en janvier, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-4-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en janvier. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-4"></span>Figure 3: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en janvier. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en février de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-5"></span>Figure 4: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en février de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en février, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-6"></span>Figure 5: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en février, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en février. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-7"></span>Figure 6: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en février. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en mars de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-8"></span>Figure 7: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en mars de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en mars, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-9"></span>Figure 8: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en mars, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en mars. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-10"></span>Figure 9: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en mars. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en avril de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-11"></span>Figure 10: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en avril de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en avril, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-12"></span>Figure 11: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en avril, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-13-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en avril. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-13"></span>Figure 12: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en avril. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-14-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en mai de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-14"></span>Figure 13: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en mai de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-15-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en mai, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-15"></span>Figure 14: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en mai, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-16-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en mai. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-16"></span>Figure 15: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en mai. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-17-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en juin de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-17"></span>Figure 16: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en juin de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-18-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en juin, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-18"></span>Figure 17: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en juin, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-19-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en juin. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-19"></span>Figure 18: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en juin. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-20-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en juillet de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-20"></span>Figure 19: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en juillet de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-21-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en juillet, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-21"></span>Figure 20: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en juillet, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-22-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en juillet. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-22"></span>Figure 21: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en juillet. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-23-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en août de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-23"></span>Figure 22: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en août de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-24-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en août, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-24"></span>Figure 23: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en août, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-25-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en août. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-25"></span>Figure 24: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en août. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-26-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en septembre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-26"></span>Figure 25: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en septembre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-27-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en septembre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-27"></span>Figure 26: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en septembre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-28-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en septembre. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-28"></span>Figure 27: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en septembre. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-29-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en octobre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-29"></span>Figure 28: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en octobre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-30-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en octobre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-30"></span>Figure 29: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en octobre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-31-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en octobre. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-31"></span>Figure 30: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en octobre. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-32-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en novembre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-32"></span>Figure 31: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en novembre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-33-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en novembre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-33"></span>Figure 32: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en novembre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-34-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en novembre. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-34"></span>Figure 33: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en novembre. Les années sont représentées par un gradient de couleur.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-35-1.png" alt="a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en décembre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales." width="576" />
<p class="caption"><span id="fig:unnamed-chunk-35"></span>Figure 34: a. Cinq premières composantes principales de l'ACp fonctionnelle réalisée sur les profils de température en décembre de 0.5 à 50m de profondeur à banyuls. b. Eboulis des valeurs propres associées aux composantes principales.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-36-1.png" alt="Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en décembre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-36"></span>Figure 35: Deux premières composantes principales de l'ACP fonctionelle réalisée sur les profils verticaux de température en décembre, exprimées comme perturbation de la moyenne. La ligne noire représente la moyenne. La ligne grise représente l'effet d'une perturbation positive. La ligne orange représente l'effet d'une perturbation négative.</p>
</div>

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-37-1.png" alt="Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en décembre. Les années sont représentées par un gradient de couleur." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-37"></span>Figure 36: Représentation bi-dimensionnelle des scores de l'ACP réalisée sur les profils de température à banyuls en décembre. Les années sont représentées par un gradient de couleur.</p>
</div>


<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-38-1.png" alt="Comparaison des scores des 2 premières composantes principales pour chaque mois." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-38"></span>Figure 37: Comparaison des scores des 2 premières composantes principales pour chaque mois.</p>
</div>




<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-40-1.png" alt="Comparaison des profils verticaux moyens de températures en noir avec intervalles de confiance à 95% associés en gris pour chaque mois." width="672" />
<p class="caption"><span id="fig:unnamed-chunk-40"></span>Figure 38: Comparaison des profils verticaux moyens de températures en noir avec intervalles de confiance à 95% associés en gris pour chaque mois.</p>
</div>



<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-42-1.png" alt="Boxplot des valeurs de températures en fonction de la profondeur par intervalle de 3m pour chaque mois. " width="672" />
<p class="caption"><span id="fig:unnamed-chunk-42"></span>Figure 39: Boxplot des valeurs de températures en fonction de la profondeur par intervalle de 3m pour chaque mois. </p>
</div>
