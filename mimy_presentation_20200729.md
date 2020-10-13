How to measure selection on pathogen evasion
========================================================
author: Richèl J.C. Bilderbeek
date:
autosize: true

[https://github.com/richelbilderbeek/mimy_presentation_20200729](https://github.com/richelbilderbeek/mimy_presentation_20200729) ![](CC-BY-NC-SA.png)

Libraries
========================================================

The Bianchi, Bilderbeek, Bogaart Question:


```r
library(bbbq)
library(pureseqtmr)
library(testthat)
```


```r
cat(1 + 1)
```

```
2
```


***

![](frans_and_geert.jpg)

Arms race: detection
========================================================

![](er_was_eens_het_leven_resized.jpg)

Research question
========================================================

Do pathogens evolve to avoid host detection?

![](evasion.jpg)

Natural selection in pathogens
========================================================

There are loci that are strongly selected upon [1]

***

![](F4_400.png)

No natural selection on host evasion
========================================================

We cannot detect natural selection on host evasion
from full genomes [2].

***

![](han2018_1b.png)

Sure!
========================================================

 * Vital loci are rare
 * Vital loci differ
 * Selection on detecting rare$^2$ loci is unlikely

Therefore, detection of vital loci is unlikely

***

![](waiting_to_small.jpg)

Transmembrane helices
========================================================

 * TMHs are general structures that all hosts and pathogens have
 * Hence, selection on TMHs is more likely

Can we detect this?

***

![](tmh_50.png)

Detect
========================================================

Transition|$n_{obs}$|$p$
----------|---------|----
D -> D    |?        |?
D -> U    |?        |?
U -> D    |?        |?
U -> U    |?        |?

***

 * $n_{obs}$: number observed
 * $p$: chance we observe this

If $p < 0.05$, we claim there is selection working upon it

Measure oserved number
========================================================

 * Take a strong dataset
 * Count the state transitions

Transition|$n_{obs}$|$p$
----------|---------|---
D -> D    |350      |?
D -> U    |200      |?
U -> D    |100      |?
U -> U    |350      |?

***

![](covid_75.png)

Measure the chance we observe this
========================================================

Given a focal TMH AA sequence


```r
f <- paste0(
  "CMIGFVIYLLFGFI",
  "LSLMCVFVLFVILILI"
)
expect_true(is_tmh(f))
```

and haplotype ...


```r
h <- "HLA-B*27:05"
```

***

... we can predict if it is detected:
































```
Error in bbbq::check_mhc_haplotype_name(mhc_haplotype = mhc_haplotype,  : 
  argument "ic50_prediction_tool" is missing, with no default
```
