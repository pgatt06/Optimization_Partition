# Résolution du problème BCPk : Méthodes et Modélisations

## Contexte général

Le **BCPk (Balanced Connected Partition into k Parts)** est un problème d’optimisation combinatoire consistant à partitionner un graphe pondéré $G = (V, E, W)$ en $k$ sous-ensembles $V_1, \dots, V_k$ :
- tels que chaque sous-graphe induit $G[V_i]$ soit **connexe**,  
- et que le poids total $W(V_i) = \sum_{v \in V_i} W[v]$ soit **équilibré** entre les classes.

L’objectif est de **maximiser le poids minimal** parmi les $k$ classes :
$$
\max \ \min_{i=1,\dots,k} W(V_i)
$$
C’est un problème NP-difficile. Trois formulations mathématiques ont donc été proposées et résolues avec **Gurobi** via **JuMP**, utilisant des techniques de flots, de coupes et de séparation incrémentale de contraintes.

---

## ⚙️ Méthode 1 : Formulation à Flots (Q1)

### Principe
Cette méthode modélise le problème sous la forme d’un **réseau de flots** :
- On introduit $k$ **sources fictives**, chacune alimentant une composante connexe.
- Chaque sommet réel doit **recevoir un flux égal à son poids $W[v]$**.
- Les arcs du graphe transportent le flux sous contrainte de capacité, activés par des variables binaires.

### Modèle
- **Variables** :
  - $f_{uv} \ge 0$ : flux sur l’arc $(u,v)$,
  - $y_{uv} \in \{0,1\}$ : activation de l’arc.
- **Objectif** :  
  Maximiser le flux total sortant de la première source (le poids minimal).
- **Contraintes** :
  - Conservation du flux aux sommets réels : entrée – sortie = poids du sommet.
  - Un seul arc sortant par source et un seul « père » par sommet.
  - Activation des arcs : $f_{uv} \le w(G) \, y_{uv}$.
  - Ordre non décroissant des flots entre sources.

### Intérêt
Cette formulation garantit **la connectivité** de chaque composante via la circulation du flux.  
Elle est intuitive mais peut devenir **lourde** pour des graphes denses à cause du grand nombre d’arcs et de variables.

---

## ✂️ Méthode 2 : Formulation par Coupes et Séparation (Q2)

### Principe
La deuxième approche repose sur un **modèle compact** avec des variables binaires $x[v,i]$ indiquant l’appartenance du sommet $v$ à la classe $i$, complété par une **procédure de séparation** pour assurer la connectivité.

### Modèle de base
- **Variables** : $x[v,i] \in \{0,1\}$.
- **Objectif** :  
  Maximiser le poids total de la première classe $W(V_1)$, sous contrainte d’équilibre croissant :
  $$
  \sum_v W[v] x[v,i] \le \sum_v W[v] x[v,i+1}.
  $$
- **Contrainte d’unicité** : chaque sommet appartient à **au plus une classe**.

### Séparation de connectivité
- À chaque itération, on résout le modèle actuel.
- Si deux sommets non reliés $(u,v)$ appartiennent à la même classe avec $x[u,i] + x[v,i] > 1$,  
  un **problème de flot maximum / min-cut** est résolu pour détecter une violation de connectivité.
- Une **contrainte de coupe** est ajoutée pour interdire cette configuration :
  $$
  x[u,i] + x[v,i] - \sum_{z \in S} x[z,i] \le 1
  $$
  où $S$ est l’ensemble de sommets séparant $u$ et $v$.

### Intérêt
Cette approche est plus **scalable** que la formulation à flots, car elle ajoute dynamiquement les contraintes nécessaires à la connectivité.  
Cependant, la convergence peut être lente, surtout sans heuristique de filtrage.

---

## ⚡ Méthode 3 : Formulation par Coupes BOOSTÉE (Q3)

### Principe
Cette méthode améliore la précédente sur deux points :
1. **Boost de la séparation** :  
   Les tests de connectivité ne sont effectués que lorsque les deux sommets ont une probabilité élevée d’appartenir à la même classe  
   ($x[u,i] > 0.5$ et $x[v,i] > 0.5$), ce qui réduit considérablement le nombre de coupes inutiles.
2. **Ajout d’inégalités croisées ("cross inequalities")** :  
   Pour tout **4-cycle sans cordes** $(a,b,c,d)$, on impose :

   $$
   x[a,i_1] + x[c,i_1] + x[b,i_2] + x[d,i_2] \le 3, \quad \forall i_1 \ne i_2
   $$
   Ces contraintes évitent certaines configurations symétriques et non connexes.

### Intérêt
Cette version est **plus robuste et plus rapide** :
- Moins de coupes inutiles grâce au **BOOST**.
- Meilleure relaxation linéaire grâce aux **cross inequalities**.
- Généralement, elle produit des partitions connexes optimales plus efficacement.

---

## 🔍 Comparaison synthétique

| Méthode | Type de formulation | Connexité assurée | Vitesse | Robustesse |
|----------|--------------------|-------------------|----------|-------------|
| **Flow (Q1)** | Réseau de flots + binaire | Implicite | Lente sur grands graphes | Fiable |
| **Cuts (Q2)** | Variables d’affectation + coupes | Séparation dynamique | Moyenne | Bonne |
| **Cuts Boostée (Q3)** | + Boost + inégalités croisées | Séparation ciblée | Rapide | Excellente |

---

## 📈 Vérification et Visualisation

Le script prévoit :
- une fonction `check_connectivity(E, classes)` pour vérifier que chaque sous-graphe induit est bien **connexe** ;
- et une fonction `plot_partition(E, classes)` pour visualiser la **partition colorée** du graphe obtenu.

---
