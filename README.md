# Master's Thesis
This project contains the preliminary version of the *Thesis* presented in fulfillment of the requirements for the degree of **_Master of Science in Statistics and Data Science for Social, Behavioral and Educational Sciences_** at [KU Leuven](https://onderwijsaanbod.kuleuven.be/2020/opleidingen/e/CQ_50550147.htm#activetab=diploma_omschrijving), during the academic year 2020-2021.


## Status
**_In progress_**

The preliminary document can be donwloaded from [**here**](https://raw.githubusercontent.com/jriveraespejo/thesis/master/thesis_JoseRivera.pdf)


## General information

### Generalized Linear Latent and Mixed Models: _method, estimation procedures, advantages, and applications to educational policy_

Local independence is one of the key assumptions of Item Response Theory (IRT) models, and it is comprised of two parts: (i) local item independence, and (ii) local individual independence [3, 21]. In the former case, the assumption entails that the individual’s response to an item does not affect the probabilty of endorsing another item, after conditioning on the individual’s ability. While in the case of the latter, the assumption considers that a individual’s response to an item is independent of another person’s response to that same item [44].

The literature has shown that IRT models are not robust to the violation of local independence. The transgression of the assumption affects model parameter estimates, inflates measurement reliabilities and test information, and underestimates standard errors (see Yen [56], Chen and Thissen [9], and Jiao et al. [24]).

However, item response data arising from educational assessments often display several type of dependencies, e.g. testlets, where items are constructed around a common stimulus [53]; the measurement of multiple latent traits within individuals [44]; cluster effects, where correlation among individuals results from the sampling and measurement mechanism used to gather the data [43]; among others. A good motivating example, that will permeate this research, is the reading comprehension sub-test, from the Peruvian public teaching career national assessment. The test is designed to measure three hierarchically nested sub-dimensions of reading comprehension: literal, inferential and reflective abilities. Furthermore, the items are bundled together in testlets related to a common text or passage. Finally, multiple cluster effects are present, e.g. at the region, and district level, just to mention a few. 

Recent studies have proposed IRT Dual Dependence Models (DDM) to deal with the teslet and individual clustering dependencies observed in the data [15, 14, 13, 24, 11, 12, 44, 7]. The majority of these representations have been developed under the bayesian framework, and they are similar in parametrization to multilevel models. On the other hand, an almost independent line of research, the Generalized Linear Latent and Mixed Models (GLLAMM) [38, 40, 48, 41], have extended the capabilities of hierarchical models on the estimation of multiple latent traits at different hierarchical levels. These
developments have been motivated mostly under the frequentist framework, and they are similar in parametrization to a hierarchical Structural Equation Model (SEM).

While the initial sense is that both developments are independent of each other, following their literature, one can easily notice that they share more than a resemblance. Both follow a multilevel/hierarchical multidimensional approach to account for the clustering of persons within the samples and the item bundles (DDM), or the latent structures within the individuals (GLLAMM). However, it is important to point out that in some cases the model parametrization between the two developments differs in a way, that some of them appear to be useful only under their specific contexts. Fortunately, their integration under the bayesian framework is not only trivial, but it can be motivated under either
type of models.

The benefits of the integration revolves around two facts: (i) the educational data often presents all of the aforementioned dependencies and more, as in the motivating example; and (ii) as it was point out in the second paragraph, in order to reach appropriate conclusions from the parameter estimates, IRT models need to account for all of these dependencies. The latter is particularly important as, more often than not, a researcher is interested in produce inferences at the structural level of the model, i.e. how a different set of manifest variables explain the variability in the latent variables, or how the latent variables explain other manifest or latent variables, at different levels. As an example, one
might be interested in finding evidence if the latent “abilities” of the teachers are explained  y their initial educational conditions, i.e. if they were educated in an institute, university or both. The main purpose of this would be to identify the type of teacher that might benefit more from the in-service training, offered by the national educational authorities, making the intervention cost-effective.

From the previous description, one can infer that the proposed IRT representation would be complex and highly dimensional. Moreover, as educational assessments are usually scored in a binary way (the individual either endorse or not the item), and because not all individuals are assessed by all the items, the model will be trained on sparse data. From the modeling perspective, neither of previous points presents a challenge for the bayesian framework. However, it has long been recognized that complex parametrizations, that allow this powerful modeling schemes, introduce pathologies that make Markov Chain
Monte Carlo methods (MCMC) face performance challenges [16, 17, 33, 34, 4]. This is highly relevant, because in order to make inferences about the posterior distribution of the parameters, the chains need to achieve three requirements, highly related to the performance of the method: stationarity, convergence and good mixing [30].

Throughout the bayesian IRT literature one often finds that four solutions are offered to ensure the fulfillment of the previous requirements, and they can be classified in two broad groups: (i) solutions that involves changing the settings of the MCMC methods, and (ii) solutions that involves readjusting the bayesian model.

In the first category we find two proposals: (a) increasing the number of iterations per chains, with large burn-in and thinning processes, and (b) designing model-specific MCMC algorithms. The easiest to implement and more prevalent in the literature is the former, e.g. Fujimoto [14] used chains with 60, 000 iterations, where 15, 000 were discarded and the remaining were thinned in jumps of 3; while Fujimoto [13] used 225, 000 iterations, with burn-in of 30, 000 and thinning with jumps of 15. Among the drawbacks
of this solution are the large computational times; the user involvement on deciding the specific setting for the process, which could be different for different parameters in the same model; and finally, the lack confidence that larger chain iterations actually produce a proper posterior investigation. On the other hand, several authors have developed high-tech MCMC algorithms that aim to optimize their performance within a particular class of models [34]. In these cases, the developers re-evaluate not only the use of the programming language, with the purpose of speeding and improving performance (e.g. Fujimoto [14]); but also the inclusion of ad-hoc model assumptions, like the ones used in staple software
developments like Mplus [31] or Stata [39]. It is clear from the previous that this solution is not accessible to all researchers, either because of the lack of programming skills, or the restrictive cost of access involved in acquiring the software. But more importantly, this solutions are not always applicable to a wider framework of similar models [34].

In the second category we also find two proposed solutions: (a) re-write the model in an alternative parametrization, and (b) encode prior information through the prior distributions. On both solutions the purpose is to ensure the identification of the parameters within the model, which helps to stabilize the MCMC procedure [18]. An example of the former is Fujimoto [14], who decomposed the items’ discriminatory parameters into overall and specific item discriminations. For the latter, Fujimoto [15] used informative
priors also for the items’ discrimination parameters.

More often than not, researchers use two or more of the aforementioned solutions to achieve a good performance in the chains. However, as point out by Betancourt and Girolami [4], even the most simple hierarchical models present formidable pathologies, that no simple correction can be performed to visit the posterior distribution properly. This is true no matter the rotation/rescaling of the parameter, or the amount of data. In this context, several authors [16, 17, 33, 34, 4] showed that prior information can be included in the model, not only through the prior distributions but also by encoding it in the model itself, changing the posterior sampling geometries by removing the dependence of the parameters on other sampled parameters, therefore favoring the performance of MCMC chains.

Given all of the above, the present research will focus its attention on showing how easy is to account for all of the dependencies that educational data often display, under the GLLAMM framework. Furthermore, given that only the literature related to gaussian hierarchical models have shown the benefits of changing the posterior sampling geometries, through the use of the non-centered parametrization [16, 17, 33, 34, 4], it seems sensible to provide a similar assessment for nonlinear hierarchical models, and in particular the ones with latent stochastic processes like IRT [34]. Finally, the research will apply the newfound knowledge to a data coming from a large Teacher’s standardized educational assessments from Peru.

### Technology
The computational implementation of the method will be developed in Stan [65] and R [46, 64].


### Supervisors
* [Geert Molenberghs](https://www.kuleuven.be/wieiswie/nl/person/00056633)
* [Wim Van den Noortgate](https://www.kuleuven.be/wieiswie/nl/person/00006844)


## Contact
Created by [Jose Manuel Rivera Espejo](http://linkedin.com/in/jriveraespejo)
