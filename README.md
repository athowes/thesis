
# Bayesian spatio-temporal methods for small-area estimation of HIV indicators

> \[!NOTE\]  
> If you’re interested in my advice on doing a PhD, see this [blog
> post](https://athowes.github.io/posts/2024-04-01-phd-advice/)!

> \[!TIP\]  
> I have made the completed
> [corrections](https://athowes.github.io/thesis/correct.pdf) for this
> thesis available. There are various other
> [enhancements](#enhancements) I would ideally like to make, but as
> they say “a good thesis is a done thesis”.

PhD in [Modern Statistics and Statistical Machine
Learning](https://statml.io/) at [Imperial College
London](https://www.imperial.ac.uk/).

Available as: [HTML](https://athowes.github.io/thesis/) or
[PDF](https://athowes.github.io/thesis/main.pdf).

Supervised by: [Seth Flaxman](https://sethrf.com/) and [Jeff
Eaton](https://www.imperial.ac.uk/people/jeffrey.eaton).

Progress towards ending AIDS as a public health threat by 2030 is not
being made fast enough. Effective public health response requires
accurate, timely, high-resolution estimates of epidemic and demographic
indicators. Limitations of available data and statistical methodology
make obtaining these estimates difficult. I developed and applied
Bayesian spatio-temporal methods to meet this challenge. First, I used
scoring rules to compare models for area-level spatial structure with
both simulated and real data. Second, I estimated district-level HIV
risk group proportions, enabling behavioural prioritisation of
prevention services, as put forward in the UNAIDS Global AIDS Strategy.
Third, I developed a novel deterministic Bayesian inference method,
combining adaptive Gauss-Hermite quadrature with principal component
analysis, motivated by the Naomi district-level model of HIV indicators.
In developing this method, I implemented integrated nested Laplace
approximations using automatic differentiation, enabling use of this
algorithm for a wider class of models. Together, the contributions in
this thesis help to guide precision HIV policy in sub-Saharan Africa, as
well as advancing Bayesian methods for spatio-temporal data.

![](figures/hiv-aids/naomi-continent.png)

## Chapters

|     | Title                                                                                                                        | GitHub repository                                             | Journal                                                                                                           |
|-----|------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| 1   | [Introduction](https://athowes.github.io/thesis/introduction.html)                                                           |                                                               |                                                                                                                   |
| 2   | [The HIV/AIDS epidemic](https://athowes.github.io/thesis/hiv-aids.html)                                                      |                                                               |                                                                                                                   |
| 3   | [Bayesian spatio-temporal statistics](https://athowes.github.io/thesis/bayes-st.html)                                        |                                                               |                                                                                                                   |
| 4   | [Models for areal spatial structure](https://athowes.github.io/thesis/beyond-borders.html)                                   | [`beyond-borders`](https://github.com/athowes/beyond-borders) | In preparation!                                                                                                   |
| 5   | [A model for risk group proportions](https://athowes.github.io/thesis/multi-agyw.html)                                       | [`multi-agyw`](https://github.com/athowes/multi-agyw)         | [PLOS Global Public Health](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0001731) |
| 6   | [Fast approximate Bayesian inference](https://athowes.github.io/thesis/naomi-aghq.html)                                      | [`naomi-aghq`](https://github.com/athowes/naomi-aghq)         | In preparation!                                                                                                   |
| 7   | [Conclusions](https://athowes.github.io/thesis/conclusions.html)                                                             |                                                               |                                                                                                                   |
| A   | [Appendix to models for areal spatial structure](https://athowes.github.io/thesis/models-for-areal-spatial-structure.html)   |                                                               |                                                                                                                   |
| B   | [Appendix to a model for risk group proportions](https://athowes.github.io/thesis/a-model-for-risk-group-proportions.html)   |                                                               |                                                                                                                   |
| C   | [Appendix to fast approximate Bayesian inference](https://athowes.github.io/thesis/fast-approximate-bayesian-inference.html) |                                                               |                                                                                                                   |

## Citation

If you would like to cite this work, please use:

    @phdthesis{howes23,
      author = {Howes, Adam},
      school = {Imperial College London},
      title = {Bayesian spatio-temporal methods for small-area estimation of HIV indicators},
      year = {2023}
    }

## Slides

Slides for my thesis defense are available
[here](https://athowes.github.io/thesis/slides.pdf). They may be useful
to provide a brief overview of the research. For more detailed slides,
see this
[presentation](https://athowes.github.io/thesis/slides-long.pdf).

<p align="center">
<a href="https://athowes.github.io/thesis/slides.pdf">
<img src="docs/slide-1.png" alt="The title slide!" height="300"/> </a>
</p>

## Frequently asked questions

> How can I read the thesis?

Thanks for being interested! You can read either the
[HTML](https://athowes.github.io/thesis/) or
[PDF](https://athowes.github.io/thesis/main.pdf) version. I know, an
overwhelming choice. Depending on my appetite, it may be improved
post-defense in March 2024.

> How did you format this thesis?

I used the R package
[`thesisdown`](https://github.com/ismayc/thesisdown), inspired by
[`bookdown`](https://github.com/rstudio/bookdown). So far it has been
working relatively seamlessly, so I’d recommend it! I’ve found `gt`
together with `knitr::is_html_output` and `knitr::is_latex_output` works
well to present tables nicely in multiple formats.

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
[Chapter 4](https://athowes.github.io/thesis/beyond-borders.html). This
chapter uses a variety of model comparison techniques (scoring rules,
coverage assessments, information criteria) that may also be of
interest. If you’re interested in modelling multinomial data using the
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
which risk varies by age, behaviour, and geographic area. This work
forms that basis for the “sub-national HIV estimates in priority
populations” UNAIDS [tool](https://hivtools.unaids.org/shipp/). You may
also be interested in the analysis in [Chapter
6](https://athowes.github.io/thesis/naomi-aghq.html) applying the Naomi
small-area estimation model to data from Malawi.

## Enhancements

Here are a collection of enhancements I think would improve this thesis:

### [Chapter 2](https://athowes.github.io/thesis/hiv-aids.html)

- The sections about challenges and statistical approaches used to
  overcome those challenges could be 1) better connected to the work
  done in the thesis, and 2) better integrated with existing literature.
  Doing so is relatively challenging 1) because this chapter precedes
  proper introduction of the methods used in this thesis, and is instead
  meant to provide a high-level overview, and 2) because the statistical
  approaches described e.g. “borrowing information” are relatively
  general and would be difficult, though not impossible, to credit to
  any particular works.

### [Chapter 3](https://athowes.github.io/thesis/bayes-st.html)

- The writing throughout Chapter 3 is weak in places. Particularly the
  sections on 1) deterministic Bayesian inference methods (difficult to
  know how much to say given focus of later chapters on this
  material) 2) properties of spatio-temporal data, 3) aspects of the
  survey section.
- In Chapter 3 I follow other authors in using the notation
  $u_k(w_{ki})$ to refer to random effects. It would be good to connect
  up this functional notation with specifying random effects as $u_i$
  rather than some function of some covariates.
- In Chapter 3, it would be nice to include a figure illustrating the
  DHS sampling procedure.

### [Chapter 4](https://athowes.github.io/thesis/beyond-borders.html) and [Appendix A](https://athowes.github.io/thesis/models-for-areal-spatial-structure.html)

- The simulation study was run using 250 replicates. As you can see from
  the plots showing the mean and standard errors, this sample size was
  insufficient to distinguish between models in some cases. All the more
  so zooming into single areas. It would be relatively simple to
  increase the sample size here, but this wasn’t done in the interests
  of time.
- For the simulation study on the four vignette geometries, the
  lengthscale priors are mis-specified with respect to the true
  lengthscale. This seems like an odd choice. Likely these experiments
  should be rerun simulating data from a more suitable lengthscale than
  the value 2.5 used currently.
- In Chapter 4, it would be useful to frame the Besag model (and BYM2,
  if possible) in terms of an equivalent kernel. I believe that the
  technical vignette [Paciorek
  (2008)](https://www.stat.berkeley.edu/users/paciorek/research/techVignettes/techVignette2.pdf)
  does this.
- In Chapter 4, calculating the DIC and WAIC values for each of the
  fitted models would be informative as to the possible benefits of the
  other model comparison techniques used in the chapter. This would
  require writing a function to take a model fitted using `TMB` or
  `aghq` and output the model comparison criteria. Likely the best
  approach would be to use samples, as this is the most transferable
  way.

### [Chapter 5](https://athowes.github.io/thesis/multi-agyw.html) and [Appendix B](https://athowes.github.io/thesis/a-model-for-risk-group-proportions.html)

- Too little emphasis is placed on the HIV prevalence and HIV incidence
  results, as compared to the HIV risk group results. For example,
  continental choropleths could be produced for the these
  epidemiological quantities as well.
- Chapter 5 could benefit from more discussion of the statistical
  results and conclusions from the work. Some of this work is already
  done in my [retrospective blog
  post](https://athowes.github.io/posts/2023-04-21-risk-group-retrospective/)
  about the work.

### [Chapter 6](https://athowes.github.io/thesis/naomi-aghq.html) and [Appendix C](https://athowes.github.io/thesis/fast-approximate-bayesian-inference.html)

- A more thorough description of the approximations to the Laplace
  approximation used by [Rue, Martino, and Chopin
  (2009)](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2008.00700.x)
  and [Wood
  (2020)](https://academic.oup.com/biomet/article-abstract/107/1/223/5572662).
  It would be instructive to implement these approximations for a simple
  example. Additionally, a more complete description of the “augmenting
  the latent field” issue (perhaps a simple example would be instructive
  here too).
- For all figures showing the use of a quadrature rule, it could be
  informative to compute and display the resulting integral estimate.
  When compared to a known truth, this would make demonstrate the value
  of e.g. adaption.
- Inclusion of some broader discussion of the value of automatic
  differentiation for INLA-like inference strategies. See the
  conversation I began
  [here](https://groups.google.com/g/r-inla-discussion-group/c/avPWD5ED0NM/m/b94kaaUVAgAJ)
  on the `R-INLA` Google group.
- Further detail about how automatic differentiation works might be
  helpful. This could include a simple example. See this nice blog post
  [“symbolic differentiation in a few lines of
  code”](https://reside-ic.github.io/blog/symbolic-differentiation-in-a-few-lines-of-code/)
  by Rich FitzJohn.
- Although the epilepsy example shows that the INLA results from `TMB`
  are highly comparable to `R-INLA`, they are not exactly the same. As
  such it would be valuable to provide an explanation for the possible
  causes. There are things that `R-INLA` does that I have not talked
  about. The best source of information about this is [Osgood-Zimmerman
  and Wakefield
  (2022)](https://onlinelibrary.wiley.com/doi/abs/10.1111/insr.12534).
- Running NUTS via `tmbstan` for Appendix A I found that some of the
  chains hung for a very long time. This is suggestive of the posterior
  geometries being tricky. In part this is confusing because I have
  previously run some of these models, implemented in `rstan` directly,
  without trouble. (See the `tmbstan` implementations
  [here](https://github.com/athowes/arealutils/blob/main/R/tmbstan.R),
  and the `rstan` implementations
  [here](https://github.com/athowes/arealutils/tree/main/dev/tempstan),
  all as part of the `arealutils` R package). It could be beneficial
  to 1) try to understand why it is that sometimes the chains hang 2)
  repeat the comparison using `rstan` versions of the models. The
  challenge in doing 2) is that then the guarantee that the models are
  the same is lost. That said, I have previously overcome this by taking
  parameters, evaluating their log-posterior under the `TMB` and `rstan`
  C++ templates, and making sure that they are identical up to a
  additive constant (this is on the log-scale, remember). For example, I
  did this in the case of the epilepsy GLM. Thinking about it more, this
  would be a valuable exercise in any case, to ensure that the
  implementations in the package are consistent. (I don’t think there is
  a way to evaluate the objective function corresponding to an `R-INLA`
  model, but if there were that would be great to do too.)
- For the figure in Chapter 6 showing the CCD grid, I think that these
  points
  ([produced](https://github.com/athowes/thesis/blob/master/figures/naomi-aghq/inla-grid-demo.R)
  using `rsm::ccd`) should have associated weights (and therefore be
  different `(aes(size = ...))`) but I am unsure as to how to generate
  the weights. There is a section in [Rue, Martino, and Chopin
  (2009)](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2008.00700.x)
  which could be useful.
- I ran Laplace marginals with emprirical Bayes for the Naomi ELGM but
  did not have enough time to present the results and integrate them
  into the discussion. This would be of interest to do.
- For the section about Naomi as an ELGM, it would be valuable to
  note 1) which features are possible in `R-INLA` (and exactly why or
  why not) 2) whether the feature is an important non-linearity of the
  model or more of a trivial notation issue.
- At one point I was interested in the kernel stein discrepancy as a way
  to measure distances between (samples from) distributions. I would
  have been interested to read more about these measures and their
  relation to MMD.

### General

- Writing of the results and conclusions sections for Chapters 4 and 6
  was relatively rushed. As such, it’s likely that a more thorough job
  could be done interpreting the results and connecting them to key
  takeaways.
- Rendering the PDF version, figures and tables tend to move around a
  lot, especially in the appendix. It would be good to have been control
  over this, but I know that this can be challenging in LaTeX.
- Note to check that references are displayed in a systematic way. Are
  there any better settings than the one that I have currently?
