## STEP 1

Après observations des images sources, le riz basmati est fin et allongé le riz camarguais de taille moyenne et ovale et enfin le riz japonais est plutôt rond.

<img src="img\basmati.png" style="zoom:20%;" /><img src="img\japonais.png" style="zoom:24%;" /><img src="img\camargue.png" style="zoom:20%;" />

## STEP 2

Nombre de grains de riz compté en analysant les fichiers pgm fournis.

| type de riz | nombres (4_8 / 8_4) | nombre après éliminations aux bordures (4_8 / 8_4) |
| ----------- | ------------------- | -------------------------------------------------- |
| basmati     | 141 / 116           | 124/100                                            |
| camarguais  | 132 / 111           | 112/93                                             |
| japonais    | 147 / 135           | 137/125                                            |

Exemple d'image obtenue pour le riz basmati :

![boundaryCurve_basmati_4_8_boundariesElimination](img\boundaryCurve_basmati_4_8_boundariesElimination.png)

## STEP 3&4

Voilà le résultat obtenu sur le riz japonais. 
On retrouve sur cette image à la fois l'inter-pixel boundaries et la GreedySegmentation.

![CurvesAndDss_japonais_4_8](D:\IMAC\IMAC3\S5\IMAC_DiscreteGeometry\src\TD02\img\CurvesAndDss_japonais_4_8.png)

---

![CurvesAndDss_japonais_4_8Bis](img\CurvesAndDss_japonais_4_8Bis.png)



## STEP 4, 5 & 6

### Aire

Pour calculer l'aire j'ai seulement utilisé dans ce TD la méthode 2 proposée en utilisant l'aire du polygone représentant l'objet digital (la méthode 1 ayant déjà été utilisée lors du premier TD).

Cette méthode converge bien pour une résolution multigride en ce qui concerne des objets convexes comme c'est le cas ici pour des trains de riz.

Ci dessous la répartition par histogramme des mesures effectuées sur les trois set de données.

![area](img\area.png)

On observe une répartition assez étalée pour les trois types de grains de riz. Cette mesure de l'aire ne permet donc pas de classifier efficacement les différents types de riz.

On peut malgré cela observer des distributions légèrement différentes en regardant la valeur médiane pour chaque type de riz mais ce n'est pas suffisant à mon sens.

| basmati | camargue | japonais |
| ------- | -------- | -------- |
| 8388    | 7094     | 6468,5   |

---

### Périmètre

Concernant le périmètre, j'ai utilisé dans ce TD la méthode 2 proposée en utilisant un polygone représentant l'objet digital (la méthode 1 ayant déjà été utilisée également lors du premier TD).

Cette méthode converge également bien par augmentation de la résolution multigride sur des objets convexes.

Ci dessous la répartition par histogramme des mesures effectuées sur les trois set de données.

![perimeter](img\perimeter.png)

On peut remarquer ici que les grains de riz japonais se distinguent par un périmètre plus élevé ce qui permet de les caractériser et les différencier ici de manière efficace. 

---

### Circularité

Pour caractériser la circularité des grains de riz j'ai utilisé la formule suivant : 

```cpp
const double circularity = (4 * M_PI * area) / (perimeter * perimeter);
```

Cette formule à la particularité d'être égale à $1$ pour un cercle  et s'exprime sous forme d'un quotient de l'air et du périmètre de la forme considérée.

Ci dessous la répartition par histogramme des mesures effectuées sur les trois set de données.

![circularity](img\circularity.png)

la circularité est étalée pour chaque type de grains de riz. Cela peut venir de ma manière d'analyser les grains de riz pas assez poussée et qui introduit des aberrations de mesure (lorsque que deux grains sont très proches ou alors coupé en morceaux).

Malgré cela et en observant les valeurs médianes pour les trois types de grains, on remarque que les grains de riz japonais se détachent et la circularité permet également ici de les distinguer.

| basmati | camarguais | japonais |
| ------- | ---------- | -------- |
| 1, 624  | 1,543      | 2,363    |

Malheureusement, les grains de riz basmati et camarguais se ressemblent beaucoup et les différentes valeurs d'aire, de périmètre et de circularité ne semble pas être suffisante pour les caractériser efficacement.

## annexes

Vous pourrez retrouver les différentes valeurs utilisées pour réaliser les diagrammes  dans les fichiers au format csv joints avec ce rapport.

Code disponible sur github via le liens suivant : https://github.com/dsmtE/IMAC_DiscreteGeometry