# Glucocorticoids negatively relate to body mass on the short-term in a free-ranging ungulate

Lucas D. Lalande<sup>1</sup>, Emmanuelle Gilot-Fromont<sup>1,2</sup>, Jeffrey Carbillet<sup>1,3</sup>, François Débias<sup>1</sup>, Jeanne Duhayer<sup>1,4</sup>, Jean-Michel Gaillard<sup>1</sup>, Jean-François Lemaître<sup>1</sup>, Rupert Palme<sup>5</sup>, Sylvia Pardonnet<sup>1,6</sup>, Maryline Pellerin<sup>7</sup>, Benjamin Rey<sup>1</sup>, Pauline Vuarin<sup>1</sup>

<sup>1</sup> Université de Lyon, Université Lyon 1, CNRS, Laboratoire de Biométrie et Biologie Evolutive UMR 5558, F-69622 Villeurbanne, France\
<sup>2</sup> Université de Lyon, VetAgro Sup, Marcy l'Etoile, France\
<sup>3</sup> Current institution: Institute of Ecology and Earth Sciences, University of Tartu, Tartu, 51014, Estonia\
<sup>4</sup> Current institution: Office Français de la Biodiversité, Direction de la Recherche et de l'Appui Scientifique, Service Départemental de la Haute-Garonne, Villeneuve-de-Rivière, 31800, France\
<sup>5</sup> Unit of Physiology, Pathophysiology and Experimental Endocrinology, Department of Biomedical Sciences, University of Veterinary Medicine, Veterinärplatz 1, 1210 Vienna, Austria\
<sup>6</sup> Current institution: Université de Strasbourg, CNRS, IPHC UMR 7178, Strasbourg, France\
<sup>7</sup> Office Français de la Biodiversité, Direction de la Recherche et de l'Appui Scientifique, Service Conservation et Gestion Durable des Espèces Exploités, Châteauvillain, 52210, France\

Glucocorticoids negatively relate to body mass on the short-term in a free-ranging ungulate\
*Oikos*, **2023**\
DOI: [10.1111/oik.09769](https://doi.org/10.1111/oik.09769)

Corresponding author:\
[lalande.luke\@gmail.com](mailto:lalande.luke@gmail.com)

## Project description

Using longitudinal data on roe deer *Capreolus capreolus* from two populations facing markedly different environmental contexts, we tested whether baseline GC levels negatively correlate with body mass – a trait positively associated with demographic individual performance – on the short- to long-term. In support, higher baseline GC concentrations were asso- ciated to lighter body mass, both measured during the same capture event, in adults of both populations. Overall, we showed that despite the marked environmental and demographic differences between populations and despite the between-sex differences in life history (i.e. reproductive tactics), the relationship between body mass and GCs is consistent across environmental contexts, but might differ according to the life his- tory stage of an individual.

## Content

The deposit contains:

- A description of the dataset used and a link to download it,
- R scripts

### Dataset

The dataset can be found on Dryad, DOI: [10.5061/dryad.pnvx0k6sk](https://doi.org/10.5061/dryad.pnvx0k6sk) 

The anonymised dataset has to be used to replicate the article mentioned above.

It includes:

-   idkit: the ID of the kit used for sampling

-   numind: Roe deer ID

-   pop: Population (3F: Trois-Fontaines, C: Chizé)

-   cohorte: cohort

-   qualite_cohorte: cohort environmental quality defined as the mean fawn mass standardised for the date of capture

-   sexe: sex

-   annee: year of capture

-   datejulienne: Julian date of capture

-   qualite_an: year of capture environmental quality defined as the mean fawn mass standardised for the date of capture

-   ageannee: age (years)

-   masse: mass (kg)

-   id_bioch: ID for biochemical analyses

-   dureemin: time between capture and sampling (min)

-   neutroP: Proportion of neutrophil

-   lymphoP: Proportion of lymphocyte

-   FGM_Priority: Whether faecal sampling is prioritize towards stress analyses or other analyses (nutrition, parasitism, ...)

-   FGMngg: FGM levels (ng/g)

-   NA values indicates that no measurements could be performed on the field or that the age or year of birth (cohort) of an individual was unknown because it has not been captured as fawn or juvenile

### R Scripts

Below you will find a description of all scripts available and used for this paper.

- `00_repeatability.R`: Script to explore data, calculate repeatability and produce figures (Supp.) for data distribution
- `01_figure1.R`: Produces Figure 1 of the article (according to script `01_short-term_stress_bm`)
- `01_short-term_stress_bm`: Script to explore relationship between FGM and body mass in juveniles and adults
- `02_figure2.R`: Produces Figure 2 of the article (according to script `02_medium-term_stress_bm.R`)
- `02_medium-term_stress_bm.R`: Script to explore relationship between FGM and the change of body mass between two consecutive years for juveniles, growing and prime-age individuals
- `03_early-late_stress_bm.R`: Script to explore the relationship between juvenile FGM and adult body mass
- `03_figure3.R`: Produces Figure 3 of the article (according to script `03_early-late_stress_bm.R`)

## Configuration 

Necessitate `R (>= 4.2.2)` and `RStudio`

## License

This repository contains the R code used to reproduce analyses from the article:

> Lalande et al. **(2023)**. Glucocorticoids negatively relate to body mass on the short-term in a free-ranging ungulate. *Oikos.*

-   Code and scripts are released under the **MIT** license.
-   The data are available on Dryad (DOI: [10.5061/dryad.pnvx0k6sk](https://doi.org/10.5061/dryad.pnvx0k6sk)) and are governed by their own license.


## Citation

Cite the following if you refer to the article:

> Lucas D. Lalande, Emmanuelle Gilot-Fromont, Jeffrey Carbillet, François Debias, Jeanne Duhayer, Jean-Michel Gaillard, Jean-François Lemaître, et al. **(2023)**. Glucocorticoids Negatively Relate to Body Mass on the Short-Term in a Free-Ranging Ungulate. *Oikos* 2023(10):e09769. <https://doi.org/10.1111/oik.09769>.

Cite the following if using the dataset:

> Lucas D. Lalande, Emmanuelle Gilot-Fromont, Jeffrey Carbillet, François Debias, Jeanne Duhayer, Jean-Michel Gaillard, Jean-François Lemaître, et al. **(2023)**. Data from: Glucocorticoids Negatively Relate to Body Mass on the Short-Term in a Free-Ranging Ungulate. *Dryad*. <https://doi.org/10.5061/DRYAD.PNVX0K6SK>.

Cite the following if using the codes and scripts made available here:

> ICI INSERER UN LIEN ZENODO GITHUB (release) POUR LA CITATION

## Contact: 

For all queries:

Lucas D. Lalande - [lalande.luke\@gmail.com](mailto:lalande.luke@gmail.com)
