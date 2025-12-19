# CATFISH : GWAS pathway analysis pipeline

This Markdown is structured into:

- [**INTRODUCTION**](#introduction) — what the pipeline is and why multiple tests are needed  
- [**METHODS**](#methods) — paper‑ready subsections with explicit equations and assumptions  
- [**USAGE**](#usage) — installation + reproducible example commands and R snippets
- [**RESULTS**](#results) — an example analysis for soil Nitrogen GWAS

---

## INTRODUCTION

### What CATFISH is

**CATFISH** (Combining <ins>**C**</ins>auchy (ACAT) <ins>**A**</ins>daptive TFisher (soft), <ins>**F**</ins>isher, m<ins>**I**</ins>n-P, and <ins>**S**</ins>touffer for <ins>**H**</ins>olistic pathway analysis) is a multi-test pathway framework built on LD-aware MAGMA gene-level GWAS statistics adjusted for gene length and SNP density that combines ACAT, soft TFisher, Fisher, Stouffer and minP. It then uses an omnibus test based on permutation-calibrated minP or ACAT to collapse these multiple pathway tests into a single, correlation-robust enrichment p-value that is sensitive to both sparse and polygenic pathway patterns. In short, CATFISH casts a wide net across complementary tests and reels in a single pathway p-value.

CATFISH uses:

1. **MAGMA** for LD-aware **SNP → gene** inference (gene-level p-values).
2. **Multiple gene → pathway combination tests** (ACAT, Fisher, soft TFisher, Stouffer, minP).
3. A **correlation-robust LD-aware omnibus test** (permutation-calibrated minP and ACAT-O) that aggregates these tests into a single pathway-level p-value.

### Why multiple tests are needed

Pathways can be significant for statistically different reasons:

- one driver gene dominates,
- many genes show coordinated moderate enrichment,
- diffuse polygenic shift, and
- hybrids of these patterns.

No single gene-set statistic is uniformly most powerful across these various possibilities. Instead of relying on a single test, CATFISH runs several complementary pathway tests and combines them into one omnibus p-value. We summarize these patterns using a set of **pathway signal archetypes** (sparse driver, coordinated moderate enrichment, diffuse polygenic shift, hybrid driver–support, single-gene proxy), which describe different ways a pathway can appear significant.

# Pathway signal archetypes

We divide the pathway signals into a set of archetypes that describe different ways a pathway can be important in a biological sense explained in detail below. We then use a combination of statistical tests chosen to be representative to each of these behaviors and a provide a biological example for each archetype.

## Archetype I — Sparse Driver Architecture (SDA)

![Archetype I — Sparse Driver Architecture](Figures/Fig_archetypesI.png)
*Archetype I: Sparse Driver Architecture (SDA).*

**Signature:** A small number of genes are extremely significant; most genes look null.

**Gene-level p-value pattern:**

- There exists a small $K \ll G$ such that the top $K$ gene p-values are extremely small, $p_{(1)}, \dots, p_{(K)} \ll \alpha$ (e.g. $p_{(1)} \sim 10^{-7}$ or smaller),
- The remaining genes are approximately null, $p_{(K+1)}, \dots, p_{(G)} \sim \mathrm{Uniform}(0,1)$.
- This produces a sharp “elbow” in the ranked p-values (a few tiny hits followed by a long flat tail).

**Interpretation:**  
An SDA pathway is significant because **a small set of driver genes dominates the signal**, rather than broad involvement of most pathway members. This can occur when the trait-relevant biology passes through a **bottleneck** (committed step, rate-limiting enzyme, key regulator, or essential transporter) so that genetic variation concentrates its effect at a few control points. In contrast, pathway annotations typically include many additional enzymes, modifiers, and general “support” genes that may be necessary for pathway operation but do not carry strong association for the trait. Under SDA, association is therefore concentrated in the top $K$ genes, yielding very small ordered p-values $p\_(1), …, p_(K)$ followed by a long tail $p_(K+1), …, p_(G)$ that is close to uniform.

**Biological example:**  
**Aspartokinase in the aspartate-derived amino acid pathway**

The aspartate-derived amino-acid biosynthesis pathway converts **aspartate** into essential amino acids, such as **lysine**, **threonine**, **methionine**, and **isoleucine**. In plants and bacteria, the first step is catalyzed by **aspartokinase (AK)**, which phosphorylates aspartate to **aspartyl-phosphate** that feeds multiple branched network to produces several end-products. Due to its position at the starting point of the pathway, AK acts as a **flux-controlling bottleneck**. Variations in AK activity alter carbon and nitrogen flow through the whole network. Downstream enzymes, including tailoring steps, dehydrogenases, and transaminases, typically exhibit dispersed or buffered roles. Thus, the gene-level pattern aligns with SDA. A single AK gene or a few genes downstream exhibit very low $p$ values (e.g. $10^{-8}\text{-}10^{-10}$), while the majority are dispersed across the interval $(0,1)$.

**Best detectors in CATFISH:**
- **ACAT** — sensitive to a few extremely small p-values (RECOMMENDED).  
- **minP / Tippett** — targets the minimum p-value; optimal when one gene dominates.


---

## Archetype II — Coordinated Moderate Enrichment (CME)

![Archetype II — Coordinated Moderate Enrichment (CME)](Figures/Fig_archetypesII.png)
*Archetype II — Coordinated Moderate Enrichment (CME)*


**Signature:** many genes show moderate association; no single gene is extremely low.

**Gene-level p-value pattern:**
- A non-trivial fraction of genes have moderately small p-values, e.g.
  $$p_i \in [10^{-3},\,0.05]\quad\text{for many } i.$$
- The top signal is not orders-of-magnitude beyond the rest (no single-gene spike), e.g.
  $$p_{(1)} \not\ll p_{(k)}\quad\text{for small }k\ (\text{e.g., }k=5,10,20).$$
- The ranked p-values show a **broad shoulder** (many good genes) but not a sharp elbow.

**Interpretation:**  
CME signifies **collective functional engagement**. The pathway behaves like a coordinated module wherein numerous components exert small-to-moderate effects. This is expected when phenotypes emerge from distributed regulation, redundancy, and buffering/feedback as dramatic single-gene effects are rare. Statistically, enrichment arises from numerous mildly informative genes rather than single driver.

**Biological example:**  
**Cytokine / immune signaling cascades**

In numerous immunological pathways, the output is regulated not by a singular "master gene," but through **distributed modulation across various tiers** of a signaling circuit. A clear example is **TNFα / IL-1β → NF-κB**, where the engagement of upstream receptors ultimately activates **IKK complexes**, which phosphorylate **IκB inhibitors**, facilitating the nuclear translocation of NF-κB family members. Notably, the **temporal dynamics** of NF-κB activation (rapid/transient versus slower/sustained) are significantly influenced by the **stimulus class** and receptor context  (Zhao et al., 2018).The significance of these dynamics lies in the variability of transcriptional outputs influenced by stimulus, NF-κB family composition, and cell type. Core feedback and marker targets encompass genes such as **NFKBIA (IκBα)** and **TNFAIP3 (A20)**, underscoring the notion that pathway behavior is governed by numerous regulatory nodes rather than a singular switch (Zhao et al., 2018).

This “many-knobs” architecture is also apparent one layer upstream in **TNFRSF signaling**. TNFRSF receptors bind trimeric TNFSF ligands, however, increasing evidence suggests that **a single trimeric ligand–receptor complex fails** to elicit complete signaling output. Rather, for certain TNFRSF pathways (notably the **classical NF-κB pathway**), successful activation may necessitate **secondary interactions or clustering of multiple trimeric receptor complexes** (Medler et al., 2019). Mechanistically, this indicates that pathway output relies on the coordinated effects of receptor assembly/avidity, adaptor recruitment, kinase activation thresholds, and the strength of negative feedback—precisely the type of system where typical genetic variation is anticipated to produce **numerous modest perturbations** rather than a singularly significant driver.

In CATFISH terminology, this yields a CME pattern: within a cytokine/immune circuit, gene-level p-values exhibit a surplus of moderate signals (e.g., numerous genes around 10⁻³–10⁻²) without a singular extreme outlier:
$$ (p_{(1)}, p_{(2)}, \ldots) \text{ comprises numerous values approximately in the range of } 10^{-3} \text{--}10^{-2} \text{ without an extreme such as } 10^{-12}.Consequently, CME pathways are optimally represented by **evidence-accumulating tests** (e.g., Fisher, Stouffer/mean-Z, and mild-truncation/soft-TFisher), whose efficacy is enhanced when *numerous* route members exhibit moderate associations, rather than depending on a singular peak.

**Best detectors in CATFISH:**
- **Fisher’s method** (aggregates evidence across many moderately small p-values) (RECOMMENDED).
- **Stouffer / mean-Z** (gains power when many genes shift together).
- Optionally **wFisher / weighted Z** if you later add biologically informed weights.



---

## Archetype III — Diffuse Polygenic Shift (DPS)

![Archetype III — Diffuse Polygenic Shift (DPS)](Figures/Fig_archetypesIII.png)
*Archetype III — Diffuse Polygenic Shift (DPS)*


**Signature:** the pathway’s genes are, on average, slightly more associated than the genome-wide background, but almost none cross a conventional significance threshold.

**Gene-level p-value / Z pattern:**
- Most gene p-values satisfy:
  $$p_i > 0.05\quad\text{for most } i.$$
Yet the pathway shows a small but consistent shift in adjusted gene-level Z-scores:

$$
\overline{Z}_{S,\mathrm{adj}} \;=\; \frac{1}{G}\sum_{i\in S} Z_{i,\mathrm{adj}} \;\neq\; 0,
$$

often with a coherent sign (bias in one direction).
- Equivalent distributional statement:
  $$\{p_i: i\in S\}\ \text{is subtly enriched toward smaller values relative to Uniform}(0,1),$$
  but without extreme outliers.
- No sharp elbow, instead, the ranked p-values show a gentle, global downward bend relative to null.

**Interpretation:**  
DPS indicates a **global pathway bias aligned with polygenicity**: numerous genes individually exert minimal impacts in a uniform direction, resulting in pathway enrichment due to the collective subtle shift of the entire module rather than the influence of a singular "star" gene. This is the regime characterized by numerous little issues and the absence of significant challenges. Spike-hunting tests, such as minP/Tippett or highly aggressive truncation, are generally underpowered in this context due to the absence of a singular extreme p-value to leverage. Conversely, mean-/distribution-sensitive tests (Stouffer/mean-Z, Fisher, and competitive regression) are specifically formulated to identify this subtle, pervasive divergence from the null hypothesis.

**Biological example:**  

**Biological example – human height as a diffuse polygenic shift**

The height of adult humans exemplifies a highly polygenic characteristic. Initial GIANT meta-analyses involving over 180,000 individuals identified 180 loci and "hundreds of variants" associated with height. Nevertheless, they accounted for merely 10% of the phenotypic variance, despite height exhibiting approximately 80% heritability. Subsequent meta-analyses involving about 700,000 Europeans identified thousands of linked SNPs and hundreds of locations, validating the perspective that height is regulated by several common variants of minimal effect rather than a limited number of high-impact genes. A recent investigation of a "saturated map" identified over 12,000 genome-wide significant SNPs across more than 7,000 genomic segments, encompassing around 21% of the genome, thereby reaffirming Fisher’s original polygenic hypothesis for height proposed in 1918.

These variants are classified into several growth-related programs, such as chondrocyte proliferation and hypertrophy in the growth plate, extracellular matrix and cartilage organization, growth hormone and IGF-1 signaling, and morphogen pathways including TGF-β and Hedgehog, as well as overarching developmental and endocrine regulators. Analyses of height GWAS loci through pathway and gene-set evaluations have demonstrated enrichment for signaling pathways such as TGF-β and Hedgehog, as well as for genes associated with skeletal growth, growth plate regulation, and pertinent Mendelian growth disorders. Lango Allen et al. (2010) discovered that genes adjacent to height-associated variants congregate in biologically coherent pathways, including TGF-β signaling, Hedgehog signaling, and histone and growth/development gene sets. Furthermore, several SNPs near these pathway genes "narrowly miss" genome-wide significance, suggesting numerous additional sub-threshold contributors within the same modules. Guo et al. (2018) demonstrate that genes adjacent to height GWAS loci are enriched in processes and tissues pertinent to growth, including growth plate cartilage.

When examined at the level of an individual pathway (e.g., TGF-β signaling, Hedgehog signaling, or growth-plate extracellular matrix genes), this structure inherently produces a **Diffuse Polygenic Shift (DPS)** pattern. Throughout the gene set, numerous genes possess one or more prevalent variations that have minor impacts on height. Certain genes attain definitive genome-wide relevance, while many others exhibit relatively mild or nominal associations. The outcome indicates that, in contrast to random gene sets, the distribution of gene-level test statistics within these pathways exhibits a shift towards more robust evidence overall characterized by an increased number of genes with small or moderate p-values and a decreased number of genes appearing entirely null—despite the fact that most individual genes would not, on their own, substantiate a strong association claim. This scenario exemplifies the optimal application of CATFISH’s DPS-oriented detectors (Stouffer/mean-Z on adjusted gene-level Z-scores, and optionally MAGMA-style competitive regression tests). They assess whether the **average** association signal across a biologically coherent pathway is subtly yet consistently heightened in comparison to the genome-wide background.

**Best detectors in CATFISH:**
- **Stouffer / mean-Z on** `Z_adj` (unweighted; permutation-calibrated) (RECOMMENDED).
- Optionally **competitive regression-style gene-set models** (e.g., **MAGMA competitive**) if included as a component test (NOT INCLUDED IN CATFISH)


---

## Archetype IV — Hybrid Driver–Support (HDS)

![Archetype IV — Hybrid Driver–Support (HDS)](Figures/Fig_archetypesIV.png)
*Archetype IV — Hybrid Driver–Support (HDS)*


**Signature:** a few very strong genes plus a some moderately associated genes.

**Gene-level p-value pattern:**
- One or a few top genes are extremely significant, e.g.
  $$p_{(1)},\,p_{(2)} \ll 10^{-4}\quad(\text{often much smaller}).$$
- Beyond the top hits, several additional genes show moderate evidence:
  $$p_{(k)} \in [10^{-3},\,0.05]\quad\text{for multiple }k \text{ (support genes).}$$
- The remaining genes are near-null:
  $$p_{(j)} \sim \mathrm{Uniform}(0,1)\quad\text{for most other }j.$$
- A small “spike” at the top (drivers) plus a clear “shoulder” of moderately small p-values (support), then a flat tail.

**Interpretation:**  
HDS exhibits a **hierarchical pathway structure**. A small number of “driver” genes exert the most significant effects, whilst a group of genes provides modest yet persistent associations. This is prevalent in pathways where flow or signal is regulated by a limited number of control points, whereas effective route output also relies on the synchronized activity of various downstream components. This architecture is statistically positioned between Sparse Driver (single spike) and Coordinated Moderate Enrichment (wide shoulder). A distinct driver signal exists, although the significance of the pathway is augmented by supplementary moderate signals.

**Biological example** 
**LDL cholesterol as a hybrid driver–support pathway**

The regulation of LDL-cholesterol (LDL-C) exemplifies a scenario in which a limited number of genes serve as primary "drivers" within a more extensive polygenic framework. Extensive GWAS and sequencing investigations consistently identify *LDLR*, *APOB*, and *PCSK9* as fundamental genes associated with Mendelian hypercholesterolemia. Infrequent coding or splice-altering variants in these genes can alter LDL-C by approximately half a standard deviation or more and are responsible for classical familial hypercholesterolemia, significantly elevating the risk of coronary artery disease. Recent whole-genome sequencing investigations involving over 60,000 individuals reveal that even infrequent *non-coding* variations next to *LDLR* and *PCSK9* can exert effects comparable to clinically recognized exonic FH variants, hence supporting their significance as primary regulators of LDL-C homeostasis.

A supporting network of lipoprotein and cholesterol metabolism genes surrounds these drivers. GWAS of traditional lipids and nuclear magnetic resonance (NMR)-based lipoprotein characteristics have revealed numerous new loci affecting LDL particle size, concentration, and composition, including apolipoprotein clusters (*APOE/APOC*, *APOA1/A5*), hepatic lipase (*LIPC*), and transporters such as *ABCG5/ABCG8*. Individually, prevalent variants at these loci typically elucidate only a minor proportion of LDL-C variance, however, collectively, they constitute a significant segment of the genome-wide polygenic signal. Recent extensive multi-ancestry meta-analyses now identify hundreds of lipid-associated loci distributed throughout this extensive lipoprotein network.

Translating this biology into gene-level association statistics for an LDL-related trait reveals that the route is neither exclusively "sparse driver" nor entirely "coordinated moderate". Typically, one observes several robust gene-level signals at *LDLR*, *APOB*, *PCSK9*, and a few linked loci, supported by a wider array of modestly correlated genes implicated in lipoprotein assembly, remodeling, and cholesterol transport. The **Hybrid Driver–Support (HDS)** pattern is characterized by a pathway whose enhancement is indicative of both a limited number of predominant, high-impact genes and a legitimate, albeit subtler, influence from the broader metabolic framework. In CATFISH terminology, this refers to the regime when **soft TFisher** (which prioritizes the lower tail while still considering moderate p-values), along with **Fisher** and the **omnibus combination**, aligns effectively with the underlying biology.

**Best detectors in CATFISH:**
- **Soft TFisher** (tail-focused; gains power from a few strong hits *plus* additional modest hits)
- **Fisher** (accumulates evidence across the moderate support set)
- **Omnibus combination** (e.g., ACAT across ACAT/Fisher/soft-TFisher/Stouffer; or permutation-calibrated minP across methods)

---

## Archetype V — Single-Gene Proxy Pathway (SGP)

![Archetype V — Single-Gene Proxy Pathway (SGP)](Figures/Fig_archetypesV.png)
*Archetype V — Single-Gene Proxy Pathway (SGP)*


**Signature:** the pathway looks significant only because it contains one very strong gene; the remaining members look essentially null.

**Gene-level pattern:**
- The top gene has an extremely small p-value (e.g. $10^{-7}\text{--}10^{-12}$).
- The rest of the genes have p-values that look noisy / $\mathrm{Uniform}(0,1)$, with no clear excess of small p’s.
- If you ignore the top gene, there is no obvious enrichment left in the pathway.

**Interpretation:**  
The pathway effectively serves as a proxy for a singular gene-level relationship. This frequently occurs when:
- an annotation term is very small (one gene plus a couple of weakly related neighbors), or
- numerous pathway definitions redundantly incorporate the same driver gene, resulting in several "distinct" pathways activating while all reference the same underlying gene.

Biologically, the driver gene remains significant (and may also align with the SDA), although the pathway-level assertion provides no additional information beyond indicating that "this gene is strongly associated." In CATFISH, we consequently designate SGP patterns as a cautionary note: view these as gene-centric findings with associated pathway designations, rather than as proof that the entire pathway is collectively activated.

**Biological example:**  
**PAH in phenylalanine metabolism**

In humans, phenylalanine metabolism is primarily regulated by a singular bottleneck enzyme, phenylalanine hydroxylase (PAH), a hepatic monooxygenase that catalyzes the conversion of phenylalanine to tyrosine in a process dependent on tetrahydrobiopterin (BH₄) (Scriver, 2007; Elhawary et al., 2022). Classical phenylketonuria (PKU) and associated forms of hyperphenylalaninemia (HPA) occur when phenylalanine hydroxylase (PAH) activity is significantly diminished, resulting in markedly elevated blood phenylalanine levels, relatively decreased tyrosine levels, and the accumulation of neurotoxic metabolites. This culminates in the distinctive untreated PKU phenotype characterized by extremely high phenylalanine, low tyrosine, and progressive neurological impairment (Elhawary et al., 2022). Extensive clinical and molecular studies indicate that the predominant cause of HPA/PKU cases is pathogenic variants in PAH, whereas a minority results from deficiencies in BH₄ synthesis, recycling, or the PAH co-chaperone DNAJC12 (Blau et al., 2014; Himmelreich et al., 2021; Elhawary et al., 2022). In summary, within the overarching "phenylalanine metabolism" pathway, PAH represents the pivotal flux-controlling step. Both common and rare variations at PAH significantly influence systemic phenylalanine levels, while the majority of other pathway components (transporters, minor side-enzymes, cofactor recycling genes) have considerably weaker or less frequent effects at the population level.

**Best detectors in CATFISH:**
- **minP / Tippett** 
- **ACAT** 

---

## Archetype VI — Competitive Enrichment Above Background (CEAB)

![Archetype VI — Competitive Enrichment Above Background (CEAB)](Figures/Fig_archetypesVI.png)
*Archetype VI — Competitive Enrichment Above Background (CEAB)*


**Signature:** the pathway is not merely “associated”; it is enriched above the genome-wide polygenic background. It passes a MAGMA competitive test ($\beta_s > 0$), indicating that genes inside the set exhibit, on average, greater associated than those outside the set.

**Gene-level pattern:**

Let $Z_g$ be the gene-level $Z$ used by MAGMA, derived from gene p-values via a probit transform (higher $Z$ = stronger association).  
In a CEAB pathway $s$:

- The distribution of $\{Z_g : g \in s\}$ is elevated in comparison to $\{Z_g : g \notin s\}$. Equivalently, $\mathrm{mean}(Z_{\mathrm{in\_set}}) > \mathrm{mean}(Z_{\mathrm{outside\_set}})$, rather than merely $\mathrm{mean}(Z_{\mathrm{in\_set}}) > 0$.
- The signal is generally not represented by a single-gene proxy; rather, one may see multiple modestly strong genes or a small top-tail alongside an elevated mean. The crucial aspect is that the *average* in-set relationship surpasses the background level
- The pathway retains its significance following MAGMA's default covariate adjustment for confounding gene characteristics (e.g., gene size and gene density, including log-transforms), indicating that it is not "large genes/dense genes" responsible for the observed effect.

**Interpretation:**  
CEAB is the most challenging route narrative to fabricate, as it signifies authentic enrichment rather than: (i) generic polygenicity (where multiple sets appear related under independent testing), (ii) artifacts from annotation overlap, or (iii) confounding effects due to gene size or density.  

Competitive tests possess a more generalized null hypothesis because they explicitly adjust for baseline associations inherent in polygenic characteristics. In MAGMA’s Crohn’s disease case study, numerous gene sets seemed related through self-contained testing; however, only one maintained significance under the competitive hypothesis—demonstrating that CEAB identifies sets with associations that exceed expectations based on polygenicity.

**Biological example**
**From MAGMA’s Crohn’s Disease analysis**  

In the WTCCC Crohn’s disease dataset, MAGMA’s self-contained analysis identified $39$ related gene sets. However the competitive analysis recognized only one of those $39$ as enriched above background: “Regulation of AMPK activity via LKB1 (REACTOME).”

This exemplifies a standard CEAB pattern: the pathway meets the "above-background enrichment" criterion. Two supplementary sets (e.g., *Cell adhesion molecules*; *ECM-receptor interaction*) attained competitive significance solely when gene size/density correction was disabled, which MAGMA reads as (at least partially) confounding-induced inflating rather than genuine enrichment.

**Best detectors in CATFISH:**

- MAGMA competitive p-value (from `*.gsa.out`) as an external “above-background” anchor:  
  if $p_{\mathrm{comp}}$ is small and $\beta_s > 0$, interpret as CEAB-supported enrichment.

---

### Bias warning
Pathway-based analysis for assessing over-representation or enrichment is influenced by various biases, including gene size, pathway size, SNP coverage density, and linkage disequilibrium (LD) patterns, all of which must be addressed explicitly (White et al., 2020; PMC6391732). In CATFISH, we tackle these issues in three phases: (i) The SNP to gene analysis is conducted using MAGMA’s LD-aware multi-SNP model, ensuring that gene-level Z/P values account for local LD structure and SNP density; (ii) we subsequently regress MAGMA gene Z-scores against log(gene length) and log(number of SNPs), utilizing the residual-based $P_adj$ for all subsequent gene to pathway analyses, thereby eliminating any remaining dependence on gene size and SNP density; and (iii) at the pathway level, we eschew simplistic over-representation tests, opting instead to calibrate our omnibus statistics through gene-label LD-aware permutations that maintain each gene’s adjusted p-value and the observed distribution of pathway sizes, yielding enrichment p-values that are resilient to these established biases.

Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/

---

## METHODS

### Notation
## Gene-level covariate adjustment (gene length & SNP density)

Let a pathway (gene set) be denoted by $S$, containing $G = |S|$ genes indexed by $g = 1,\dots,G$.

For each gene $g \in S$, let:
- $Z_g$ denote the raw gene-level $Z$-statistic, and
- $p_g$ denote the corresponding two-sided $p$-value.

We collect these into vectors:

$$
\mathbf{Z}_S = (Z_1, Z_2, \dots, Z_G),
\qquad
\mathbf{p}_S = (p_1, p_2, \dots, p_G).
$$

To remove residual dependence on gene length and SNP density, we fit a linear model:

$$
Z_g = \beta_0 + \beta_1 \log(L_g) + \beta_2 \log(\mathrm{SNPs}_g) + \varepsilon_g.
$$

We define the adjusted $Z$-scores as the regression residuals:

$$
Z_{g,\mathrm{adj}} = \hat{\varepsilon}_g.
$$

The corresponding adjusted two-sided $p$-values are:

$$
p_{g,\mathrm{adj}} = 2\,\Phi\!\left(-\left|Z_{g,\mathrm{adj}}\right|\right),
$$

where $\Phi(\cdot)$ is the standard normal CDF.

We denote the adjusted vectors by:

$$
\mathbf{Z}_{S,\mathrm{adj}} = (Z_{1,\mathrm{adj}}, \dots, Z_{G,\mathrm{adj}}),
\qquad
\mathbf{p}_{S,\mathrm{adj}} = (p_{1,\mathrm{adj}}, \dots, p_{G,\mathrm{adj}}).
$$


---

## 1) Gene-level association statistics (SNP → gene)

For each gene $g$, MAGMA generates a gene‑level association p‑value $p_g$ by aggregating SNP‑level signals within or within a window of the gene while accounting for local LD using a reference panel. We employed MAGMA’s `multi=snp-wise` model, which integrates a SNP-wise mean test (effective for numerous small effects) and a SNP-wise top test (effective for a single strong SNP) into a unified LD-aware omnibus statistic per gene, hence ensuring the robustness of gene p-values against varying within-gene causal structures.

The workflow emcompasses:

- SNPs are mapped to genes (gene boundaries with optional windows).
- A multi‑marker gene model accounts for LD among SNPs in the genic region.
- MAGMA generates gene statistics (e.g., $Z_g$ and $p_g$).

In CATFISH, the p-values (or Z statistics) of the MAGMA gene are utilized as inputs for all pathway analyses.

---

## 2) Gene-level adjustment for gene size and SNP density

Despite LD-aware gene testing, gene-level signals may still demonstrate residual dependency on gene size and SNP density. CATFISH executes a post-hoc correction at the gene level.

Let:

- $Z_g$ be the MAGMA gene Z‑statistic,
- $L_g$ be gene length (bp),
- $S_g$ be number of SNPs mapped to the gene (e.g., $$NSNPS$$).

We fit a regression line as:

$$
Z_g = \beta_0 + \beta_1 \log(L_g) + \beta_2 \log(S_g) + \varepsilon_g.
$$

We define adjusted residual Z as:

$$
Z^{\mathrm{adj}}_g = Z_g - \widehat{Z}_g,
\quad \widehat{Z}_g = \widehat{\beta}_0 + \widehat{\beta}_1 \log(L_g) + \widehat{\beta}_2 \log(S_g).
$$

We then convert to a two‑sided adjusted p-value:

$$
p^{\mathrm{adj}}_g = 2\Phi\left(-|Z^{\mathrm{adj}}_g|\right),
$$

where $\Phi(\cdot)$ is the standard normal CDF.

---

## 3) Pathway-level test statistics (gene → pathway)

CATFISH computes multiple pathway statistics from $\{p_g\}_{g \in S}$ $\{Z_g\}_{Z \in S}$.

### 3.1 ACAT (Aggregated Cauchy Association Test)

Define the Cauchy‑transformed score for each gene:

$$
t_g = \tan\left(\pi\left(\tfrac{1}{2} - p_g\right)\right).
$$

Define non‑negative weights $w_g \ge 0$ with $\sum_{g \in S} w_g = 1$ (default $w_g = 1/G$).

The ACAT statistic is:

$$
T_{\mathrm{ACAT}}(S) = \sum_{g \in S} w_g\, t_g
= \sum_{g \in S} w_g \tan\left(\pi\left(\tfrac{1}{2} - p_g\right)\right).
$$

The combined p-value is:

$$
p_{\mathrm{ACAT}}(S) = \tfrac{1}{2} - \frac{1}{\pi}\arctan\left(T_{\mathrm{ACAT}}(S)\right).
$$

**Key property (interpretation):** ACAT is asymptotically dominated by the smallest p-values, and is therefore sensitive to **sparse driver** architectures.

---

### 3.2 Fisher’s method

Fisher’s statistic is:

$$
T_{\mathrm{Fisher}}(S) = -2\sum_{g \in S} \log(p_g).
$$

Under independence,

$$
T_{\mathrm{Fisher}}(S) \sim \chi^2_{2G},
\quad
p_{\mathrm{Fisher}}(S) = 1 - F_{\chi^2_{2G}}\left(T_{\mathrm{Fisher}}(S)\right),
$$

where $F_{\chi^2_{2G}}(\cdot)$ is the $\chi^2$ CDF with $2G$ degrees of freedom.

**Key property (interpretation):** Fisher is sensitive to **coordinated moderate enrichment** (many moderately small p-values).

---

### 3.3 Adaptive Soft TFisher (data-adaptive tail focusing via oTFisher)

A practical limitation of fixed- $\tau$ soft TFisher is that the appropriate tail-focus is contingent upon the *unknown* pathway topology (dense versus sparse, weak versus strong). TFisher proposes a **omnibus, data-adaptive** selection of truncation and weighting parameters, termed **oTFisher**, which autonomously identifies the most advantageous configuration for the observed p-value distribution (ref).

#### Soft TFisher family (recap)

For a pathway $S$ with adjusted gene-level p-values $\{p_g\}_{g\in S}$ and a soft-threshold parameter $\tau\in(0,1]$, the **soft TFisher** statistic is

$$W^{\mathrm{soft}}(S;\tau)=\sum_{g\in S}\left[-2\log(p_g)+2\log(\tau)\right]_{+},\qquad (x)_+=\max(x,0)$$

This is the $\tau_1=\tau_2=\tau$ special case of the TFisher family and implements a *continuous* down-weighting near the cutoff (soft vs hard truncation). 

#### Data-adaptive $\tau$ (oTFisher)

Let $\mathcal{T}=\{\tau_1,\dots,\tau_m\}$ be a small grid of candidate thresholds (e.g. a few small/medium/large values; TFisher shows that a sparse grid is usually sufficient). 

For each $\tau_j\in\mathcal{T}$, compute:

1) the statistic $W^{\mathrm{soft}}(S;\tau_j)$, and  
2) its **null p-value**
$$p_{\tau_j}(S)=\Pr\!\left(W^{\mathrm{soft}}(\cdot;\tau_j)\ge W^{\mathrm{soft}}(S;\tau_j)\mid H_0\right),$$
using the TFisher null calculation for $W_n(\tau_1,\tau_2)$ with $\tau_1=\tau_2=\tau_j$. 

Then define the **adaptive soft TFisher omnibus** as the minimum across the grid:

$$p_{\mathrm{aTF}}(S)=\min_{\tau\in\mathcal{T}} p_{\tau}(S)$$

This is precisely the oTFisher concept (“select the most substantial TFisher p-value among candidate $(\tau_1,\tau_2)$ configurations”), tailored to the **soft-thresholding** line $\tau_1=\tau_2$. Due to the dependence of the $p_{\tau_j}(S)$ (originating from the same ordered p-values), oTFisher offers an analytic calibration for the minimum across the grid by employing a multivariate normal approximation of the vector comprising component TFisher statistics (calculated through multivariate normal probabilities). In the soft scenario when $\tau_1=\tau_2=\tau$, the mean and covariance are simplified, allowing for efficient computation of the multivariate normal (MVN) probability, such as by Genz-style MVN cumulative distribution function (CDF). 

In CATFISH terminology, this provides a *adaptive tail sensor* without the necessity of hard-coding a single $\tau$ and without the need for extensive permutation for the *within-method* calibration (we still overlay our MVN/perm framework for the final omnibus across methods).

A minimal grid that usually works well in practice is something like:

$$\mathcal{T}=(\{0.01\,0.05\,0.5\,1\})$$

which oTFisher explicitly uses in its soft-thresholding omnibus examples, and which covers “rare-ish hits”, “moderate tail”, and “nearly Fisher”.


---

### 3.4 Stouffer's method (mean-Z; diffuse polygenic shift)

Stouffer's technique consolidates gene-level **Z** statistics instead of p-values and exhibits increased sensitivity to DPS. In CATFISH, the gene-level Z input is sourced directly from MAGMA’s gene output (`ZSTAT`). Significantly, MAGMA’s Z-scale is optimally understood as a **association-strength score** (a probit transformation of the gene p-value), indicating that greater positive values signify stronger evidence of association, rather than indicating the direction of effect as trait-increasing or trait-decreasing. Consequently, the natural pathway-level Stouffer test in this context is **one-sided (greater)**, assessing the enrichment of positive association strength inside the pathway.

**Default (unweighted) Stouffer:**

For a pathway $$S$$ with $$G=|S|$$ genes,

$$Z_{\mathrm{stouffer}}(S) = \frac{1}{\sqrt{G}} \sum_{g \in S} Z_g$$

$$p_{\mathrm{stouffer}}(S) = 1 - \Phi\!\left(Z_{\mathrm{stouffer}}(S)\right),$$

where $$\Phi(\cdot)$$ is the standard normal CDF.

**Optional (weighted) Stouffer:**

CATFISH optionally supports nonnegative per-gene weights $$w_g$$ (e.g., based on SNP count or other gene-level quantities). In that case,

$$Z_{\mathrm{stouffer}}^{(w)}(S)=\frac{\sum_{g\in S} w_g\, Z_g}{\sqrt{\sum_{g\in S} w_g^2}},\qquad p_{\mathrm{stouffer}}^{(w)}(S)=1-\Phi\!\left(Z_{\mathrm{stouffer}}^{(w)}(S)\right)$$

A **two-sided** option can be reported for completeness when a genuinely *signed* gene-level Z is available:

$$p_{\mathrm{stouffer}}^{(2\text{-sided})}(S)=2\,\Phi\!\left(-\left|Z_{\mathrm{stouffer}}(S)\right|\right),$$

but this is not the default for MAGMA-style association-strength Z scores.

**Permutation / MVN calibration (dependence-aware):** Because genes within pathways can be correlated (LD and genomic proximity), CATFISH can calibrate Stouffer p-values by recomputing $$Z_{\mathrm{stouffer}}$$ (or $$Z_{\mathrm{stouffer}}^{(w)}$$) under a null generated by either (i) global resampling of gene Z’s, or (ii) MVN simulation using a pathway-specific gene–gene correlation matrix. For one-sided enrichment, the empirical p-value is computed as the fraction of null draws with $$Z_{\mathrm{stouffer}}^{\mathrm{null}} \ge Z_{\mathrm{stouffer}}^{\mathrm{obs}}$$.





---

### 3.5 minP / Tippett (single-gene proxy diagnostic)

Define the pathway's minimum gene p-value:

$$
p_{\min} = \min_{g \in S} p_g
$$

Under independence, Tippett's combined p-value is:

$$
p_{\mathrm{tippett}} = 1 - (1 - p_{\min})^G
$$

However, due to the correlation of gene p-values (LD, shared biology), CATFISH generally considers **minP as a diagnostic** and/or adjusts it through permutation if needed. 

Let:
- $$p_{\min}^{\mathrm{obs}}$$ be the observed minimum p-value
- $$p_{\min}^{(b)}$$ be the minimum p-value in permutation $$b$$
- $$B$$ be the total number of permutations

Then the permutation p-value is:

$$
p_{\min}^{\mathrm{perm}} = \frac{1 + \text{count}(p_{\min}^{(b)} \leq p_{\min}^{\mathrm{obs}})}{B + 1}
$$

where $$\text{count}(\cdot)$$ counts the number of permutations satisfying the condition.

**Notes:**
- Highly sensitive to **single extreme genes**, so it flags **Sparse Driver (SDA)** and especially **Single-Gene Proxy (SGP)** pathways.
- Use "remove top gene and recompute" as the practical SGP check.
---

## 4) Correlation among pathway tests

All pathway statistics above are functions of the same gene-level p-values $$\{p_g\}$$ and are therefore intrinsically correlated. For a given pathway $$S$$, define $$T_{\mathrm{ACAT}},\$$ $$T_{\mathrm{Fisher}},\$$ $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau),\$$ $$T_{\mathrm{Stouffer}},\$$ $$T_{\mathrm{minP}},\$$
where

- $$T_{\mathrm{ACAT}}$$ is the ACAT statistic on $$S$$ (Section 3.1),
- $$T_{\mathrm{Fisher}}$$ is the Fisher / sum–log–p statistic (Section 3.2),
- $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau)$$ is the soft–TFisher statistic with truncation parameter $$\tau$$ (Section 3.3),
- $$T_{\mathrm{Stouffer}}$$ is the (possibly weighted) Stouffer / mean-$$Z$$ statistic, and
- $$T_{\min} = \min_{g \in S} p_g$$ is the Tippett / minP statistic.

Under the null, this vector is **not** jointly independent, because each component is a deterministic function of the same set of p-values $$\{p_g\}$$ (or equivalently the same ordered p-values $$\{p_{(k)}\}$$).

Intuitively:

- **Fisher** and **soft TFisher** are both based on $$-\log(p_g)$$, with TFisher acting like a truncated / reweighted Fisher focused on $$p_g \le \tau$$. Whenever the sum of $$-\log(p_g)$$ is large (many small p-values), both statistics tend to be extreme in the same direction.

- **Stouffer** first maps p-values to $$Z$$-scores via $$Z_g = \Phi^{-1}(1 - p_g/2)$$ and then averages or sums them. This is a monotone transform of the same $$\{p_g\}$$, so pathways with many small p-values will simultaneously inflate Fisher, TFisher, and Stouffer.

- **ACAT** and **soft TFisher** both emphasize the **lower tail**: ACAT uses $$\tan\{\pi(0.5 - p_g)\}$$ (very heavy weight on the smallest p’s), while TFisher heavily weights small $$p_g$$ through $$-\log(p_g)$$ below $$\tau$$. A pathway with one or a few very small p-values will push both $$T_{\mathrm{ACAT}}$$ and $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau)$$ in the same direction.

- **minP** is literally the smallest p-value $$p_{(1)}$$. It is almost perfectly aligned with the extreme-tail behavior of ACAT and the most aggressive truncation versions of TFisher: whenever $$p_{(1)}$$ is tiny, ACAT, TFisher, and often Fisher / Stouffer will also look unusually extreme.

On top of this, linkage disequilibrium and shared biology induce correlation among the gene-level inputs themselves (the $$\{p_g\}$$ and corresponding $$\{Z_g\}$$ are not independent across genes), which further strengthens dependence among all of the pathway-level statistics.

**Implication.** Because $$T_{\mathrm{ACAT}},\$$ $$T_{\mathrm{Fisher}},\$$ $$T_{\mathrm{TF}}^{\mathrm{soft}}(\tau),\$$ $$T_{\mathrm{Stouffer}},\$$ $$T_{\mathrm{minP}},\$$
is strongly dependent, naïve combination rules that assume independence between component tests (e.g. simple Bonferroni or analytic minP under independence) would generally be mis-calibrated and often anti-conservative. This is why CATFISH uses **permutation-calibrated minP across methods** (and treats ACAT-O as a model-based sensitivity analysis), so that the joint null distribution of these correlated statistics is learned empirically rather than assumed.


---

## 5) Omnibus pathway p-value across methods

Each pathway in CATFISH is evaluated by a panel of complementary gene → pathway tests (ACAT, Fisher, soft TFisher, Stouffer, and minP), each tuned to a different pattern of gene-level signal. Because these five statistics are all functions of the same adjusted gene p-values, they are strongly correlated and can disagree in which pathways they call “most” enriched. Rather than choosing a single favorite test a priori, we combine their evidence into a **single omnibus pathway p-value** using two complementary strategies: (i) an analytic ACAT-O combination across method-level p-values, and (ii) a permutation-calibrated minP across methods that provides a dependence-robust, “best-of-tests” summary. The following subsections describe these two omnibus layers.


### 5.1 Omnibus ACAT (ACAT-O across methods)

Let $$p_1, \dots, p_5$$ denote the five component p-values for pathway $$S$$, and let weights $$v_j \ge 0$$ satisfy $$\sum_{j=1}^5 v_j = 1$$ (default $$v_j = 1/5$$).

Define the ACAT-O Cauchy statistic

$$T_{\mathrm{omni,ACAT}}(S) = \sum_{j=1}^5 v_j \,\tan\{\pi(0.5 - p_j)\},$$
which under the global null is approximately standard Cauchy. The analytic omnibus p-value is
$$p_{\mathrm{omni,ACAT}}(S)= 0.5 - \frac{1}{\pi}\arctan\{T_{\mathrm{omni,ACAT}}(S)\}.$$

This ACAT-O layer is most sensitive when **at least one** component test (e.g. ACAT for sparse drivers, Fisher/Stouffer for coordinated enrichment, soft TFisher for hybrid patterns, or minP for hard single-gene hits) is strongly significant, even if the others are only modest or null.

---

### 5.2 Omnibus minP across methods with permutation 

To obtain a complementary, more conservative summary that explicitly accounts for correlation between component tests, we also compute a minimum-p omnibus statistic across methods: 

$$T_{\mathrm{omni,min}}(S) = \min_{j \in \{1,\dots,5\}} p_j= \min\({ p_{\mathrm{ACAT}}(S), p_{\mathrm{Fisher}}(S),p_{\mathrm{TF}}(S;\tau), p_{\mathrm{Stouffer}}(S),p_{\mathrm{minP}}(S) \})$$ 

Because the $$p_j$$ are correlated (they are all computed from the same set of gene-level p-values), we calibrate $$T_{\mathrm{omni,min}}(S)$$ by gene-label permutation. 

For each of $$B$$ permutations: 

- randomize gene labels across pathways (preserving the empirical distribution of gene-level p-values and pathway sizes),
- recompute all five component tests for each pathway,
- record, for permutation $$b = 1,\dots,B$$, $$T_{\mathrm{omni,min}}^{(b)}(S) = \min_j p_j^{(b)}(S)$$

The permutation-based omnibus p-value is then 

$$\hat p_{\mathrm{omni,min}}(S)= \frac{1 + \text{count}\{\,b : T_{\mathrm{omni,min}}^{(b)}(S)\le T_{\mathrm{omni,min}}(S)\,\}}{B + 1}$$ 

We use $$\hat p_{\mathrm{omni,min}}(S)$$ as the **primary** omnibus p-value in the main analyses because it 

1. protects nominal type-I error under *arbitrary* dependence between component tests, and
2. explicitly targets the "best" component test for each pathway while accounting for the fact that this best test is effectively chosen post hoc via the min across methods.

The analytic ACAT-O p-value $$p_{\mathrm{omni,ACAT}}(S)$$ is reported alongside as a **higher-power, model-based sensitivity analysis**, highlighting pathways that are consistently strong across methods or dominated by a single very informative test.

---

## 6) Multiple testing correction

Across all pathways, omnibus p-values $\{p_{\mathrm{omni}}(S)\}$ are adjusted using Benjamini–Hochberg FDR:

$$
q_{\mathrm{BH}}(S) = \mathrm{BH}\left(p_{\mathrm{omni}}(S)\right).
$$

Because each pathway yields a **single** omnibus p-value, no additional penalty is required for the number of component tests.

---

## 7) Omnibus across methods with resampling calibration (global + LD-aware MVN)

For each pathway $$S$$, we compute multiple *component* pathway p-values from gene-level evidence and then combine them into a single **omnibus** pathway p-value. Because component tests reuse the same gene-level inputs and are therefore correlated, the omnibus p-value is optionally calibrated under the null using either (i) **global gene-set resampling** or (ii) an **LD-aware MVN** simulation driven by MAGMA gene–gene correlations.

### 7.1 LD-aware SNP $$\rightarrow$$ gene inputs via MAGMA

GWAS summary statistics were aggregated to gene-level association statistics using MAGMA’s LD-aware SNP-to-gene model (multi-model SNP-wise gene analysis) with a reference LD panel. SNPs were mapped to genes using a symmetric $$\pm 25$$ kb window. This yields per-gene association measures such as $$P_g$$ and $$Z_g$$ that are calibrated while accounting for within-gene LD.

To reduce confounding by gene size and SNP density, gene-level evidence can be adjusted by regressing on $$\log(\text{gene length})$$ and $$\log(\(\text{SNPs})$$, and using residual-based size-adjusted gene p-values as pathway inputs (when provided as an alternative p-value column). In the implementation, p-based component tests use the preferred p-value column (e.g., an adjusted column if present, otherwise raw MAGMA gene p-values).

### 7.2 Component pathway tests

For each pathway $$S$$ with member genes $$g \in S$$, we compute the following component pathway p-values:

1. **ACAT (gene-level)**  
2. **Fisher (gene-level)**  
3. **Adaptive soft TFisher (gene-level)**  
4. **Gene-level minP (Sidák)**  
5. **Stouffer Z (optional; from gene Z-scores)**  
6. **MAGMA competitive gene-set p-value (optional; reported separately)**  

All p-values are numerically stabilized by clipping to $$[p_{\min},\,1-p_{\min}]$$ (e.g., $$p_{\min}=10^{-15}$$)
to avoid numerical issues in log/tangent transforms.

### 7.3 Omnibus across methods (ACAT-O or minP-O)

Let the component p-value set be

$$\mathcal{P}(S)=\{p_{\mathrm{ACAT}}(S),\,p_{\mathrm{Fisher}}(S),\,p_{\mathrm{TFisher}}(S),\,p_{\mathrm{minP}}(S),\,p_{\mathrm{Stouffer}}(S)\},$$

where the Stouffer term is included only when gene Z-scores are available. (MAGMA competitive is reported separately and is excluded from the resampling-calibrated omnibus by default.)

We define two omnibus options:

**(i) ACAT-omnibus (`omnibus="ACAT"`)**

$$
T_{\mathrm{omni,ACAT}}(S)
=\frac{1}{K}\sum_{j=1}^K \tan\!\Big(\pi\big(\tfrac{1}{2}-p_j(S)\big)\Big),
\qquad
p_{\mathrm{omni,ACAT}}(S)
=\tfrac{1}{2}-\frac{1}{\pi}\arctan\!\Big(T_{\mathrm{omni,ACAT}}(S)\Big).
$$

where $$\{p_j(S)\}_{j=1}^K$$ are the available component p-values.

**(ii) minP-omnibus (`omnibus="minP"`)**

$$T_{\mathrm{omni,min}}(S) \;=\; \min_j p_j(S),\qquad p_{\mathrm{omni,min}}(S) \;=\; 1-\big(1-T_{\mathrm{omni,min}}(S)\big)^{K}$$

Because component tests are correlated, the analytic omnibus p-value is treated as a convenient summary and
(optionally) calibrated via resampling (Sections 7.4–7.5).

---

## 7.4 Global gene-set resampling calibration (`perm_mode="global"`)

Global resampling calibrates the omnibus under a conservative null that preserves pathway size and the empirical
marginal distributions of gene-level statistics while breaking pathway membership.

### 7.4.1 Gene pool construction

Let $$\mathcal{G}$$ be the set of genes with non-missing pathway inputs:
- a p-value pool $$\{P_g : g\in \mathcal{G}\}$$ from the preferred p-value column, and
- if Stouffer is enabled, a matched Z pool $$\{Z_g : g\in \mathcal{G}\}$$ from `z_col`,
filtered using the same gene-level mask so $$(P_g, Z_g)$$ remain paired by gene identity.

### 7.4.2 Resampling scheme (per pathway)

For each pathway $$S$$ with $$|S|=d$$ and for $$b=1,\dots,B$$:

1. **Sample indices (gene-label resampling).**  
   Draw a length-$$d$$ index vector $$I^{(b)} = (i_1,\dots,i_d)$$ uniformly from $$\{1,\dots,|\mathcal{G}|\}$$.  
   Sampling is **without replacement** when $$d \le |\mathcal{G}|$$; if $$d > |\mathcal{G}|$$ it falls back to
   sampling **with replacement**.

2. **Construct resampled gene evidence (paired sampling).**  
   Build

   $$P^{(b)} = (P_{i_1},\dots,P_{i_d}),$$

   and if Stouffer is enabled,

   $$Z^{(b)} = (Z_{i_1},\dots,Z_{i_d}),$$

   using the *same* sampled indices $$I^{(b)}$$.  
   **Paired $$(P_g, Z_g)$$ resampling:** this enforces that Stouffer and the p-based components are driven by
   the *same* resampled genes in each replicate.

3. **Recompute component tests on the resampled set.**  
   Using $$P^{(b)}$$ compute $$p_{\mathrm{ACAT}}^{(b)}(S)$$, $$p_{\mathrm{Fisher}}^{(b)}(S)$$,
   $$p_{\mathrm{TFisher}}^{(b)}(S)$$ (including the same $$\tau$$-grid search), and $$p_{\mathrm{minP}}^{(b)}(S)$$.  
   Using $$Z^{(b)}$$, compute $$p_{\mathrm{Stouffer}}^{(b)}(S)$$ under the same alternative as the observed test.
   (For simplicity and stability, the global Stouffer null is typically treated as **unweighted** even if
   weights are used in the observed Stouffer statistic.)

4. **Compute the resampled omnibus.**  
   Combine the resampled component p-values using the same omnibus rule (ACAT-O or minP-O) to obtain
   $$p_{\mathrm{omni}}^{(b)}(S)$$.

### 7.4.3 Empirical calibration

Let $$p_{\mathrm{omni}}(S)$$ be the observed (analytic) omnibus p-value for pathway $$S$$. The global-resampling
calibrated omnibus p-value is computed with the +1 correction:

$$\hat p_{\mathrm{omni,global}}(S)=\frac{1+\left|\{\,b : p_{\mathrm{omni}}^{(b)}(S)\le p_{\mathrm{omni}}(S)\,\}\right|}{B+1}$$

This resampling directly captures dependence among component tests because all component statistics are
recomputed on the *same* resampled gene sets.

---

## 7.5 LD-aware MVN calibration using MAGMA gene–gene correlations (`perm_mode="mvn"`)

MVN calibration models correlation among gene-level statistics within each pathway using a Gaussian copula
parameterized by MAGMA-derived gene–gene correlations.

### 7.5.1 Pathway-specific correlation matrix $$R_S$$

For each pathway $$S=\{g_1,\dots,g_d\}$$, we construct a $$d\times d$$ matrix $$R_S$$ from a sparse MAGMA
correlation pairs file (three columns: `gene1`, `gene2`, `r`):
- set $$(R_S)_{ii}=1$$,
- fill $$(R_S)_{ij}=r_{ij}$$ when the pair $$(g_i,g_j)$$ is present,
- default missing pairs to 0,
- clip correlations to a safe bound (e.g., $$|r|\le 0.999$$).

When requested (`make_PD=TRUE`), $$R_S$$ is projected to the nearest positive definite correlation matrix
(e.g., via a nearPD procedure) before simulation.

### 7.5.2 MVN simulation and copula p-values (per pathway)

For $$b=1,\dots,B$$, simulate

$$Z^{(b)} \sim \mathcal{N}(0, R_S)$$

**Single latent draw drives both Z- and p-based components.**  
We derive p-values for p-based component methods from the *same* simulated $$Z^{(b)}$$ using a Gaussian-copula
mapping:

- **Uniform marginals (default; `mvn_marginal="uniform"`):**

  $$U^{(b)}=\Phi(Z^{(b)}),\qquad P^{(b)}=2\min(U^{(b)},1-U^{(b)}),$$

  yielding marginal Uniform$$(0,1)$$ p-values while preserving dependence via $$R_S$$.

- **Empirical marginals (optional; `mvn_marginal="empirical"`):**
  map $$U^{(b)}=\Phi(Z^{(b)})$$ through an empirical null quantile function derived from a genome-wide p-value
  pool, allowing conservative/deflated empirical marginals to be reflected while preserving dependence.

### 7.5.3 Recomputing component tests and omnibus under MVN

For each replicate $$b$$:
- compute p-based components (ACAT/Fisher/adaptive TFisher/minP) from the simulated p-vector $$P^{(b)}$$,
  using the same numerical clipping and the same `tau_grid` search for TFisher;
- compute Stouffer from the simulated $$Z^{(b)}$$ directly, using the observed pathway weights $$w_g$$ if weights
  are enabled (otherwise equal weights) and the same alternative hypothesis;
- combine the resulting component p-values using the selected omnibus rule (ACAT-O or minP-O) to obtain
  $$p_{\mathrm{omni}}^{(b)}(S)$$.

### 7.5.4 Empirical calibration

The MVN-calibrated omnibus p-value is computed as

$$\hat p_{\mathrm{omni,mvn}}(S)=\frac{1+\left|\{\,b : p_{\mathrm{omni}}^{(b)}(S)\le p_{\mathrm{omni}}(S)\,\}\right|}{B+1}$$

---

## 7.6 Final omnibus p-value, MAGMA competitive, and multiple testing correction

Depending on the calibration mode, we report:
- $$p_{\mathrm{omni,analytic}}(S)$$: analytic omnibus (ACAT-O or minP-O),
- $$\hat p_{\mathrm{omni,global}}(S)$$: global-resampling calibrated omnibus (if run),
- $$\hat p_{\mathrm{omni,mvn}}(S)$$: MVN-calibrated omnibus (if run),
- `magma_pvalue`: MAGMA competitive gene-set p-value (reported separately).

The **final** omnibus p-value prioritizes the most informative calibration available:

$$
p_{\mathrm{omni,final}}(S)=
\begin{cases}
\hat p_{\mathrm{omni,mvn}}(S), & \text{if MVN calibration was performed;}\\
\hat p_{\mathrm{omni,global}}(S), & \text{else if global resampling was performed;}\\
p_{\mathrm{omni,analytic}}(S), & \text{otherwise.}
\end{cases}
$$

Finally, we apply Benjamini–Hochberg FDR correction across all tested pathways to obtain
$$q_{\mathrm{omni,final}}(S)$$ (reported as `omni_p_final_BH`), and analogously provide BH-adjusted versions for
analytic/global/MVN omnibus p-values when present.


### 7.7 Treatment of MAGMA competitive in the omnibus

We additionally compute (and report) the MAGMA competitive gene-set p-value for each pathway as a separate
summary statistic (`magma_pvalue`). However, by default this competitive p-value is **not included in the
resampling-calibrated omnibus** (`include_magma_in_perm=FALSE`) because there is no cheap, principled null
generator for the MAGMA competitive statistic within the gene-set resampling layer.

Consequently, when resampling calibration is performed, the omnibus $$p_{\mathrm{omni}}(S)$$ (analytic and
calibrated) is computed over the **five gene-derived component tests only**
$$\{p_{\mathrm{ACAT}}, p_{\mathrm{Fisher}}, p_{\mathrm{TFisher}}, p_{\mathrm{minP}}, p_{\mathrm{Stouffer}}\}$$,
and the MAGMA competitive result is interpreted alongside the omnibus rather than being embedded within it.


---

## USAGE

# CATFISH (R package wrapper)

**CATFISH** is the R interface implementation used to run CATFISH‑style workflows on top of MAGMA, and to compute ACAT/Fisher/TFisher + omnibus pathway statistics.

---

## Installation

### 1) Install MAGMA (external dependency)

Download MAGMA from the official site and make the `magma` executable available on your `PATH`:

- https://ctg.cncr.nl/software/magma

### 2) Install CATFISH in R

```r
# install.packages("devtools")  # if needed
devtools::install_github("nirwan1265/MAGCAT")
library(MAGCAT)
```

### 3) Optional: set MAGMA path

```r
MAGCAT::magma_set_path("/full/path/to/magma")
```

---

## Conceptual workflow (end-to-end)

1. **SNP → gene (MAGMA)**
   - Prepare SNP locations (`*.snp.loc`) and gene locations (`*.genes.loc`).
   - Run MAGMA annotation and gene analysis to get gene Z and p.

2. **Gene-level adjustment (optional)**
   - Regress $Z_g$ on `log(gene_length)` and `log(NSNPS)`; derive $p^{adj}_g$.

3. **Gene → pathway tests**
   - Compute pathway p-values from adjusted gene p-values using:
     - ACAT,
     - Fisher,
     - soft TFisher (tail-focused),
     - Stouffer's test,
     - minP.

4. **Omnibus**
   - Combine pathway p-values using minP or ACAT to produce $p_{\mathrm{omni}}$.

5. **Multiple testing**
   - BH FDR (and optional Storey q-values).
  
6. Permutation
   - Use either random sampling or MVN (RECOMMENDED).

---

## MAGMA commands (typical)

```bash
# 1) Annotate SNPs to genes
magma \
  --annotate \
  --snp-loc  <snp.loc> \
  --gene-loc <genes.loc> \
  --out      <prefix>

# 2) Gene analysis (LD-aware)
magma \
  --bfile      <LD_reference_panel_prefix> \
  --pval       <gwas.pval.txt> N=<N> \
  --gene-annot <prefix>.genes.annot \
  --gene-model multi=snp-wise \
  --out        <prefix>
```

---

## Quick R example (CATFISH)

> **Note:** This is the exact end-to-end pipeline used (MAGMA → gene adjustment → pathway tests → omnibus). Paths, filenames, and column mappings should be edited to match your local files.

```r
############################################################
## CATFISH PIPELINE: MAGMA → SIZE-ADJUSTED GENE P → PATHWAYS
## (with brief explanation of the key parameters)
############################################################

## CATFISH gives you:
##  - R wrappers around the MAGMA binary (magma_annotate, magma_gene)
##  - Parallel chromosome-wise runs (chroms, n_threads)
##  - Plant-ready GFF3 → gene loc helpers
##  - Automatic gene-size/#SNP adjustment
##  - A suite of pathway tests + an omnibus layer
##
## Below is the end-to-end flow, with parameter explanations in comments.
############################################################


############################################################
## 1. Build MAGMA gene location file from GFF3
############################################################

gff_path     <- "/Users/.../Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
gene_loc_out <- "/Users/.../maize.genes.loc"

gff3_to_geneloc(
  gff        = gff_path,
  out        = "inst/extdata/maize.genes.loc",
  chr_prefix = "chr"  # strips "chr" prefix so MAGMA’s CHR matches PLINK CHR
)

## - You can call species = "maize" in CATFISH , it automatically
##   uses inst/extdata/maize.genes.loc 

## - You can also use any gff3 file for your organism of choice
## - Parses the GFF3, extracts gene features,
## and writes a MAGMA-ready-to-use gene location file (GENE, CHR, START, STOP, STRAND…).


############################################################
## 2. SNP → gene annotation with MAGMA (magma_annotate wrapper)
############################################################

stats_file <- "/Users/.../raw_GWAS_MLM_3PC_N.txt"

magma_annotate(
  stats_file     = stats_file1,
  rename_columns = c(
    CHR    = "chr",    # your column "chr" → MAGMA expects "CHR"
    SNP    = "rs",     # your column "rs"  → "SNP"
    POS    = "ps",     # your column "ps"  → "POS"
    PVALUE = "p_wald"  # your p-value column → "PVALUE"
  ),
  species    = "maize",        # uses built-in maize.genes.loc from step 1
  out_prefix = "N_maize_MLM",  # prefix for MAGMA output files
  out_dir    = "annot",        # directory to write MAGMA .genes.annot
  window     = c(25, 25)       # +/- 25 kb around each gene for SNP mapping
)

## - magma_annotate() builds the MAGMA command line and calls the MAGMA binary for you.
## - You just supply stats_file and remember to rename your columns or do it here


############################################################
## 3. Gene-level MAGMA (multi = snp-wise) with R wrapper
############################################################

# Plink based bed/bim.fam files
bfile      <- "/Users/.../all_maize2"   # PLINK basename: .bed/.bim/.fam

## 3A. Chromosome-wise run using NMISS (per-SNP sample size)
## NMISS is the number of missing genotypes
magma_gene(
  bfile      = bfile,
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = stats_file,
  n_total    = 3539,             # optional global N if NOBS/NMISS absent
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NMISS  = "n_miss"            # MAGMA interprets this as per-SNP sample size
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_genes_by_chr",
  gene_model = c("multi=snp-wise"),  # multi-parameter SNP-wise model. This works best for CATFISH. Combines top and mean SNPs. 
  chroms     = 1:10,                 # run MAGMA separately for chr 1–10
  n_threads  = 10                    # run up to 10 MAGMA jobs in parallel
)

## 3B. Chromosome-wise run using NOBS (per-SNP N)
## NOBS is the total number of genotypes used
magma_gene(
  bfile      = bfile,
  gene_annot = "annot/N_maize_MLM.genes.annot",
  stats_file = stats_file2,
  rename_columns = c(
    CHR    = "Chr",
    SNP    = "SNP",
    POS    = "Pos",
    PVALUE = "P.value",
    NOBS   = "nobs"             # per-SNP N; preferred over n_total
  ),
  out_prefix = "N_maize_MLM",
  out_dir    = "magma_multi_snp_wise_genes_by_chr_N_maize",
  gene_model = c("multi=snp-wise"),
  chroms     = 1:10,
  n_threads  = 10
)

## - magma_gene() is a R wrapper around the MAGMA binary.
##   * Handles all flags, temp files, and error checking.
##   * Additional can parallelize the analysis
##   * chroms + n_threads to run multiple chromosomes in parallel.
##   * Only requires a minimal rename_columns spec instead of reformatting


############################################################
## 4. Combine per-chromosome MAGMA gene outputs
############################################################

# Load the file
files <- sprintf(
  "/Users/.../magma_multi_snp_wise_genes_by_chr_N_maize/N_maize_MLM_chr%d.multi_snp_wise.genes.out",
  1:10
)

gene_list <- lapply(files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

genes_all_raw <- do.call(rbind, gene_list)

# Ensure the gene p column is named "P"
colnames(genes_all_raw)[9] <- "P"

# For genes appearing on multiple chromosomes or windows, keep smallest P
o         <- order(genes_all_raw$GENE, genes_all_raw$P)
genes_all <- genes_all_raw[o, ]
genes_all <- genes_all[!duplicated(genes_all$GENE), ]

# Optionally sort by CHR and START
if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

write.table(
  genes_all,
  file      = "/Users/.../magma_N_maize.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

## - This combines 10 per-chromosome MAGMA outputs into a single file for ease of use


############################################################
## 5. Gene length extraction + SNP density bias adjustment
############################################################

## 5A. Extract gene lengths from GFF3
gff3_path <- system.file(
  "extdata", "GFF3",
  "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3",
  package = "MAGCAT"
)

maize_gene_len <- get_gene_lengths(
  gff3_file  = gff3_path,
  output     = TRUE,
  output_dir = "inst/extdata",
  file_name  = "Zea_mays_gene_lengths.tsv"
)

## Output includes:
##   gene_id, chr, start, end, length


## 5B. Adjust gene-level Z/P for gene length & #SNPs

genes_all <- read.table(
  "/Users/.../magma_N_maize.txt", # from previous step
  header = TRUE, stringsAsFactors = FALSE
)

# Adjsut the pvalue based on gene length and number of snps
adj_out <- magcat_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = maize_gene_len,
  gene_col     = "GENE",
  p_col        = "P",
  z_col        = "ZSTAT",    # raw MAGMA Z
  len_gene_col = "gene_id",
  len_col      = "length"
  # nsnp_col   = "NSNPS"     # #SNPs per gene
)

genes_adj <- adj_out$genes    # includes Z_raw, Z_adj, P_adj, log_gene_length, log_nsnp
lm_fit    <- adj_out$fit      # lm(Z_raw ~ log_gene_length + log_nsnp)

write.csv(genes_adj, "genes_adj.csv", row.names = FALSE)

## - MAGMA gene p-values are known to correlate weakly with gene size and SNP density.
## - magcat_adjust_gene_p():
##   * Fits a linear model: Z_raw ~ log(gene length) + log(#SNPs).
##   * Uses residuals (Z_adj) as “size/SNP-adjusted” Z-scores.
##   * Converts Z_adj back to P_adj = 2*pnorm(-|Z_adj|).
## - You then use P_adj in gene→pathway tests, reducing bias toward big genes.


############################################################
## 6. Load pathway definitions from PMN/CornCyc or use saved files
############################################################

maize_pw <- magcat_load_pathways(
  species  = "maize",
  gene_col = "Gene-name"  # column in the PMN gene-set file that matches MAGMA gene IDs
)

############################################################
## 7. Omnibus combining methods (ACAT or minP)
############################################################

## All tests use gene_results + species/pathways to:
##  - find genes per pathway,
##  - take their p-values (raw P or adjusted P_adj),
##  - compute a pathway-level p per method.

omni_minp <- omni_pathways(
  gene_results      = genes_adj,
  species           = "maize",
  gene_col          = "GENE",
  p_col             = "P_adj",
  effect_col        = "Z_adj",
  is_onetail        = FALSE,
  ptrunc            = 0.05,
  min_p             = 1e-15,
  do_fix            = TRUE,
  omnibus           = "ACAT",      # minP or "ACAT"
  B_perm            = 10000L,
  seed              = 123,
  perm_mode    = "mvn",       # mvn or 
  magma_genes_out = "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_N_maize.txt",
  remove_singletons = TRUE,
  output            = TRUE,
  out_dir           = "magcat_omni_full"
)

## In a single call omni_pathways gives you:
##   - acat_p       : ACAT pvalue per pathway
##   - fisher_p    : Fisher pvalue per pathway
##   - tpm_p        : truncated Fisher (soft) pvalue per pathway
##   - stouffer_p   : Stouffer pvalue per pathway
##   - minp_gene_p  : minP pvalue per pathway
##   - omni_p       : combination of these methods (ACAT or minP) pvalue per pathway
##   - omni_perm_p  : permutation-calibrated omnibus pvalues
##   - BH FDR for omni_p and each component
##   - omni_perm_BH  : permutation-calibrated BH pvalues using perm pvalues

############################################################
## 8. Pathway-level tests (gene → pathway)
############################################################

## You can run each pathway individually as well

### 8A. ACAT per pathway

pw_res_acat_adj <- magcat_acat_pathways(
  gene_results = genes_adj,     # adjusted Pvalue
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  B            = 10000L,            # 10,000 is good enough
  seed         = NULL,
  output       = TRUE,
  out_dir      = "acat_results"
)

### 8B. Fisher pathways

wf_res_raw <- magcat_wfisher_pathways(
  gene_results = genes_adj,   # raw P
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  effect_col   = "ZSTAT",
  #weight_col   = NULL, # If you have any other weights. 
  is_onetail   = FALSE
)

### 8C. Truncated Fisher (soft) / TFisher(soft)

soft_tf_res_adj <- magcat_soft_tfisher_pathways(
  gene_results     = genes_adj,
  species          = "maize",
  gene_col         = "GENE",
  p_col            = "P_adj",
  tau1             = 0.05,     # soft truncation threshold
  B_perm           = 10000L,  # permutations for empirical pvalue
  seed             = 123,
  analytic_logical = TRUE,
  output           = TRUE,
  out_dir          = "magcat_tfisher_soft"
)

### 8D. Stouffer (sum of Z across genes)

stouf_res <- magcat_stouffer_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  weight_col   = NULL,    # equal weights for all genes
  B_perm       = 10000L,    # permutations for empirical pvalue
  seed         = 123,
  output       = TRUE,
  out_dir      = "magcat_stouffer"
)

### 8E. Gene-level minP per pathway

minp_res <- magcat_minp_pathways(
  gene_results = genes_adj,
  species      = "maize",
  gene_col     = "GENE",
  p_col        = "P_adj",
  B_perm       = 10000L,      # permutations for empirical pvalue
  min_p        = 1e-15,
  do_fix       = TRUE,
  output       = TRUE,
  out_dir      = "magcat_minp_maize"
)

```

---


## References

- White MJ et al. *Strategies for Pathway Analysis using GWAS and WGS Data*. Current Protocols in Human Genetics (2020)
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6391732/
- de Leeuw CA et al. *MAGMA: Generalized Gene‑Set Analysis of GWAS Data*. PLoS Comput Biol (2015).
  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219
- Liu Y et al. *ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies*.  The American Journal of Human Genetics (2019).
  https://pubmed.ncbi.nlm.nih.gov/30849328/
- Zhang H et al. *TFisher: A powerful truncation and weighting procedure for combining p-values*. Annals of Applied Statistics (2020)
https://projecteuclid.org/journals/annals-of-applied-statistics/volume-14/issue-1/TFisher--A-powerful-truncation-and-weighting-procedure-for-combining/10.1214/19-AOAS1302.full
- Yoon S et al. *Powerful p-value combination methods to detect incomplete association*. Nature (2021)
  https://www.nature.com/articles/s41598-021-86465-y
- Tippett, L. H. C. *The Methods of Statistics*. London:Williams & Norgate (1931)
- Westfall, P. H., & Young, S. S. *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment*. New York: Wiley(1993)
  


