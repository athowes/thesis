
# Bayesian spatio-temporal methods for small-area estimation of HIV indicators

> \[!WARNING\]  
> Thesis under construction! Check back in a few months. Currently
> 136/~200 pages, and many/~tending to zero mistakes.

PhD in [Modern Statistics and Statistical Machine
Learning](https://statml.io/) at [Imperial College
London](https://www.imperial.ac.uk/).

Supervised by: [Seth Flaxman](https://sethrf.com/) and [Jeff
Eaton](https://www.imperial.ac.uk/people/jeffrey.eaton).

Progress towards ending AIDS as a public health threat by 2030 is not
being made fast enough. Effective public health response requires
accurate, timely, high-resolution estimates of epidemic and demographic
indicators. Limitations of available data make obtaining these estimates
difficult. I developed and applied Bayesian spatio-temporal methods to
meet this challenge. First, I examined models for area-level spatial
structure. <!-- , finding that... --> Second, I estimated district-level
HIV risk group proportions, enabling behavioural prioritisation of
prevention services, as put forward in the Global AIDS Strategy.
Finally, I developed a novel deterministic Bayesian inference method,
combining adaptive Gauss-Hermite quadrature with principal component
analysis, motivated by the Naomi district-level model of HIV indicators.
Together, the contributions in this thesis help to guide precision HIV
policy in sub-Saharan Africa, as well as advancing Bayesian methods for
spatio-temporal data.

![](figures/hiv-aids/naomi-continent.png)

## Chapters

|     | Title                                                                                   | GitHub repository                                             | Journal                                                                                                           |
|-----|-----------------------------------------------------------------------------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| 1   | [Introduction](https://athowes.github.io/thesis/introduction.html)                      |                                                               |                                                                                                                   |
| 2   | [The HIV/AIDS epidemic](https://athowes.github.io/thesis/hiv-aids.html)                 |                                                               |                                                                                                                   |
| 3   | [Bayesian spatio-temporal statistics](https://athowes.github.io/thesis/bayes-st.html)   |                                                               |                                                                                                                   |
| 4   | [Spatial structure](https://athowes.github.io/thesis/beyond-borders.html)               | [`beyond-borders`](https://github.com/athowes/beyond-borders) | In preparation!                                                                                                   |
| 5   | [A model for risk group proportions](https://athowes.github.io/thesis/multi-agyw.html)  | [`multi-agyw`](https://github.com/athowes/multi-agyw)         | [PLOS Global Public Health](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0001731) |
| 6   | [Fast approximate Bayesian inference](https://athowes.github.io/thesis/naomi-aghq.html) | [`naomi-aghq`](https://github.com/athowes/naomi-aghq)         | In preparation!                                                                                                   |
| 7   | [Future work and conclusions](https://athowes.github.io/thesis/conclusions.html)        |                                                               |                                                                                                                   |

## Citation

If you would like to cite this work, please use:

    @phdthesis{howes23,
      author = {Howes, Adam},
      school = {Imperial College London},
      title = {Bayesian spatio-temporal methods for small-area estimation of HIV indicators},
      year = {2023}
    }

## Frequently asked questions

> How can I read the thesis?

Thanks for being interested! You can read either the
[HTML](https://athowes.github.io/thesis/) or
[PDF](https://athowes.github.io/thesis/main.pdf) version. I know, an
overwhelming choice. It’s still very much a work in progress at the
moment though, so I’d recommend checking back in a few months unless you
love non-sequiturs.

> How did you format this thesis?

I used the R package
[`thesisdown`](https://github.com/ismayc/thesisdown), inspired by
[`bookdown`](https://github.com/rstudio/bookdown). So far it has been
working relatively seamlessly, so I’d recommend it!

> Are there any resources you’d recommend for an introduction to this
> area of research?

I’d recommend something like [Spatial and Spatio-temporal Bayesian
models with `R-INLA`](https://sites.google.com/a/r-inla.org/stbook/) by
Marta Blangiardo and Michela Cameletti. I have a repository with further
miscellaneous recommended
[resources](https://github.com/athowes/resources), if you are
interested.

> I’m a statistician: which parts of the thesis might interest me?

If you use spatial random effects to model areal spatial structure,
[Chapter 4](https://athowes.github.io/thesis/beyond-borders.html). If
you’re interested in modelling multinomial data using the
multinomial-Poisson transformation and structured random effects,
[Chapter 5](https://athowes.github.io/thesis/multi-agyw.html). If you
have a complicated model which is not compatible with `R-INLA`, but
would still like to use INLA-like methods, [Chapter
6](https://athowes.github.io/thesis/naomi-aghq.html). If you have such a
model, get in touch! All of the methods are compatible with any model
written in the (very general) [Template Model
Builder](https://kaskr.github.io/adcomp/Introduction.html) R package
(`TMB`).

> I’m a HIV epidemiologist: which parts of the thesis might interest me?

Primarily [Chapter 5](https://athowes.github.io/thesis/multi-agyw.html)
will be of interest to you. The [Global AIDS
Strategy](https://www.unaids.org/en/Global-AIDS-Strategy-2021-2026) sets
out goals for prioritisation of prevention programming for adolescent
girls and young women according to risk behaviour and epidemic setting.
To enable implementation of the strategy, I estimated risk group
specific population sizes, prevalence and incidences at a
district-level. I also used these estimates to evaluate the extent to
which risk varies by age, behaviour, and geographic area. You might also
be interested in the analysis in [Chapter
6](https://athowes.github.io/thesis/naomi-aghq.html) applying the Naomi
small-area estimation model to data from Malawi.
