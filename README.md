Spillover risk is multidimensional and complicated. Within a pool of competent
hosts, it is driven by a multiplicity of factors [@Plowright2017PatZoo].
Although proxies for the local richness of hosts is commonly analysed [see
*e.g.* @Anthony2017GloPat for coronaviruses], there is an argument to be made
that species who are not (known as) competent hosts of a specific virus genus
may not factor into this [@Plowright2015EcoDyn], calling for species-level
information. This is especially true as competence data increases predictive
power when the taxonomic scope of hosts of a viral family increases
[@Becker2020InfInt; @Mull2022VirIso]. Similarly, host species who share viruses
at different rates should be weighted accordingly [@Albery2020PreGlo]. In
mammals, key functional traits (for which phylogeny is a reasonable proxy) are
determinants of the spillover potential [@Olival2017HosVir]; these include,
notably, body mass, and affinity for urban environments [@Albery2022UrbMam].
Finally, especially when the pool of potential hosts spans the entire globe,
there may be local host pools that are highly unique; not having been observed
in other locations, these can act on the overall risk either by providing novel
contact opportunities, reflecting unique host-environment combinations
[@Engering2013PatHos], or facilitating rapid evolutionary changes in specialism
of their pathogens [@Agosta2010HowSpe]. In the specific case of generalist
pathogens [which is the case many viruses in the betacoronavirus genus, see
*e.g.* @MacLean2021NatSel], there is conceptual and empirical support to the
idea that these community-level mechanisms are even more important in driving
the overall risk [@Power2004PatSpi].

In this paper, we examine the biogeographic structure of bats-betacoronaviruses
associations, based on a curated dataset of all confirmed bat hosts of
betacoronaviruses. By drawing on concepts from the Geographic Mosaic Theory of
Coevolution [GMTC; @Thompson2005GeoMos], we turn these associations into a
spatially explicit additive mapping of zoonotic risk components, revealing the
extreme heterogeneity of risk at the global scale. Explicitely framing the
notion of spillover risk based on propositions from the GMTC (which is to say,
based on a framework linking interactions between species to change within
species) is a novel idea, that should be relatively general. Indeed, it only
assumes the action of well described evolutionary mechanisms. The benefit of
this approach is to provide the potential for a more dynamic and nuanced
understanding of risk: not only on ecological timescales, but also by providing
clues about which areas can change over micro-evolutionary timescales. This
provides a way to look at spatial structure by accounting for more notions than
species richness/similarity, but also a way to identify spatial areas of higher
risk.

We identify the Amazon and South-Eastern Asia as hotspots where the phylogenetic
distinctiveness of betacoronaviruses is the highest [@Anthony2017GloPat];
surprisingly, current data suggest that viral sharing between hosts is high in
the Amazon and low in South-Eastern Asia, which has the potential to result in
different evolutionary dynamics between these two regions, hinting at different
futures for their viral communities. This work is important both as a
description of the bats-betacoronavirus complex, but also because more broadly,
bats are known reservoirs for a variety of emerging viruses and other pathogens
[@Calisher2006BatImp;@Melaun2014BatPot], making balancing the needs for bat
conservation and disease prevention most likely very difficult and a source of
human-wildlife conflicts, especially in more densely populated areas
[@Rego2015AssHum; @Stone2015ManCon].

# Methods

## Known betacoronavirus hosts

We downloaded the data on bats hosts of betacoronaviruses assembled by
@Becker2022OptPre from `https://www.viralemergence.org/betacov` on Apr. 2022,
and filtered it to "known" hosts (established before the emergence of
SARS-CoV-2) and "novel" hosts (confirmed through sampling since the emergence of
SARS-CoV-2). The original database was assembled by a combination of data mining
and literature surveys, including automated alerts on the "bats" and
"coronavirus" keywords to identify novel empirical evidence of
bats-betacoronaviruses associations.

## Bats occurrences

We downloaded the rangemap of every extant bat species that was either
classified as an empirically documented host of beta-coronaviruses from the
previous step, according to recent IUCN data [@IUCN2021IucRed]. The range maps
were subsequently rasterized using the `rasterize` function from `GDAL`
[@RouaultEven2022GdaOgr] at a resolution of approximately **TK TP**. For every
pixel in the resulting raster where at least one bat host of betacoronavirus was
present, we extract the species pool (list of all bat species), which was used
to calculate the following risk assessment components: phylogenetic diversity,
bat compositional uniqueness, and predicted viral sharing risk.

## Bats phylogeography

For every pixel, we measured Faith’s Phylogenetic Diversity [@Faith1992ConEva]
based on a recent synthetic tree with robust time calibration, covering about
6000 mammalian species [@Upham2019InfMam]. Faith’s PD measures the sum of unique
branches from an arbitrary root to a set of tips, and comparatively larger
values indicate a more phylogenetic diverse species pool. We measured
phylogenetic diversity starting from the root of the entire tree (and not from
Chiroptera); this bears no consequences on the resulting values, since all
branches leading up to Chiroptera are only counted one per species pool, and (as
we explain when describing the assembly of the composite risk map), all
individual risk components are ranged in [0,1]. This measure incorporates a
richness component, which we chose not to correct for; the interpretation of the
phylogenetic diversity is therefore a weighted species richness, that accounts
for phylogenetic over/under-dispersal in some places.

## Bats compositional uniqueness

For every species pool, we measured its Local Contribution to Beta-Diversity
[@Legendre2013BetDiv]; LCBD works from a species-data matrix (traditionally
noted as $\mathbf{Y}$), where species are rows and sites are columns, and a
value of 1 indicates occurrence. We extracted the Y matrix assuming that every
pixel represents a unique location, and following best practices
[@Legendre2019SpaTem] transformed it using Hellinger’s distance to account for
unequal bat richness at different pixels. The correction of raw community data
is particularly important for two reasons: first, it prevents the artifact of
richer sites having higher importance; second, it removes the effect of overall
species richness, which is already incorporated in the phylogenetic diversity
component. High values of LCBD indicate that the pixel has a community that is
on average more dissimilar in species composition than what is expected knowing
the entire matrix, i.e. a more unique community. Recent results by
@Dansereau2022EvaEco shows that LCBD measures are robust with regards to spatial
scale, and are therefore applicable at the global scale.

## Viral sharing between hosts

For all bat hosts of betacoronaviruses, we extracted their predicted viral
sharing network [@Albery2020PreGlo]. This network stores pairwise values of
viral community similarity. To project viral sharing values into a single value
for every pixel, we averaged the pairwise scores. High values of the average
sharing propensity means that this specific extant bat assemblage is likely to
be proficient at exchanging viruses.

## Composite risk map

To visualize the aggregated risk at the global scale, we combine the three
individual risk components (phylogenetic diversity, compositional uniqueness,
and viral sharing) using an additive color model [@Seekell2018GeoLak]. In this
approach, every risk component gets assigned a component in the RGB color model
(phylogenetic diversity is green, compositional uniqueness is red, and viral
sharing is blue). In order to achieve a valid RGB measure, all components are
re-scaled to the [0,1] interval, so that a pixel with no sharing, no
phylogenetic diversity, and no compositional uniqueness is black, and a pixel
with maximal values for each is white. This additive model conveys both the
intensity of the overall risk, but also the nature of the risk as colors diverge
towards combinations of values for three risk components. Out of the possible
combinations, the most risky in terms or rapid diversification and spillover
potential is high phylogenetic diversity and low viral sharing
[@Cavender-Bares2009MerCom; @Gomulkiewicz2000HotSpo], in that this allows
multiple independent host-virus coevolutionary dynamics to take place in the
same location. In the colorimetric space, this correspond to yellow -- because
the HSV space is more amenable to calculations for feature extraction [see
*e.g.* @Keke2010StuSki], we measured the risk level by calculating the angular
distance of the hue of each pixel to a reference value of 60, and weighted this
risk level by the value component. Specifically, given a pixel with colorimetric
coordinates $(h,s,v)$, its ranged weighted risk value is

$$
v\times\left[1-\frac{\left|\text{atan}\left(\text{cos}(\text{rad}(h)), \text{sin}(\text{rad}(h))\right) - X\right|}{2\pi}\right]\,,
$$

where X is $\text{atan}\left(\text{cos}(\text{rad}(60)),
\text{sin}(\text{rad}(60))\right)$, a constant approximately equal to $0.5235$.

## Viral phylogeography and evolutionary diversification

We used the following query to pull all betacoronavirus sequence data from the
GenBank Nucleotide database except SARS-CoV-2; ("Betacoronavirus"[Organism] OR
betacoronavirus[All Fields]) NOT ("Severe acute respiratory syndrome coronavirus
2"[Organism] OR sars-cov-2[All Fields]). We added a single representative
sequence for SARS-CoV-2 and manually curated to remove sequences without the
RNA-dependent RNA polymerase (RdRp) sequence or that contained words indicating
recombinant or laboratory strains including “patent”, “mutant”, “GFP”, and
“recombinant”. We filtered over-represented taxa including betacoronavirus 1,
hCoV-OC43, Middle East respiratory syndrome coronavirus, Murine hepatitis virus,
and hCoV-HKU1. Curated betacoronavirus RdRp sequences were then aligned using
MAFFT v 1.4.0 (**Katoh and Standley 2013**, parameters in text?) and a maximum likelihood
tree reconstructed in IQ-TREE v 1.6.12 (Nguyen et al. 2015) with ModelFinder
(Kalyaanamoorthy et al. 2017) ultrafast bootstrap approximation (Hoang et al.
2018) and the following parameters (**STEPH WILL ADD**, parameters in text?).

We first tested the hypothesis that hotspots of viral diversification would
track hotspots of bat diversification. To do so, we plotted the number of known
bat hosts (specifically only those included in the phylogeny, so there was a 1:1
correspondence between data sources) against the “mean evolutionary
distinctiveness” of the associated viruses. To calculate this, we derived the
fair proportions evolutionary distinctiveness [@Isaac2007MamEdg] for each of the
viruses in the tree, then averaged these at the bat species level, projected
these values onto their geographic distributions, and averaged across every bat
found in a given pixel. As such, this can be thought of as a map of the mean
evolutionary distinctiveness of the known viral community believed to be
associated with a particular subset of bats present.

## Co-distribution of hosts and viral hotspots

Subsequently, we tested the hypothesis that the biogeography of bat
betacoronaviruses should track the biogeography of their hosts. To test this
idea, we loosely adapted a method from [@Kreft2007GloPat; @Kreft2010FraDel], who
proposed a phylogenetic method for the delineation of animal biogeographic
regions. In their original method, a distance matrix - where each row or column
represents a geographic raster’s grid cell, and the dissimilarity values are the
“beta diversity similarity” of their community assemble - undergoes non-metric
multidimensional scaling (NMDS); the first two axes of the NMDS are projected
geographically using a four-color bivariate map. Here, we build on this idea
with an entirely novel methodology. First, we measure the phylogenetic distance
between the different viruses in the betacoronavirus tree by using the
cophenetic function in `ape` [@Paradis2019ApeEnv]; subsequently, we take a
principal components analysis of that distance matrix (readily interchangeable
for NMDS in this case) to project the viral tree into an n-dimensional space. We
then take the first two principal components and, as with the evolutionary
distinctiveness analysis, aggregated these to a mean host value and projected
them using a four-color bivariate map.

# Results and discussion

## Host richness does not predict virus distinctiveness

Bats are found worldwide and are both one of the most diverse groups among
mammals (**Moratelli & Calisher, 2015**), and one of the main animal reservoir
for different strains of betacoronaviruses [@Drexler2014EcoEvo]. This has
attracted attention to areas where high diversity of bats, and therefore
presumably high diversity of betacoronaviruses, can be an important issue for
human health [@Calisher2006BatImp; @Moratelli2015BatZoo]. By overlaying the IUCN
rangempas for confirmed bat hosts of betacoronaviruses [@fig:richness; top], we
see that the the main hotspots of host richness are primarily South-Eastern
Asia, parts of Southern Europe, and to a lesser extent parts of Africa in the
-25-0 range of latitudes. The description of host richness is an important first
step towards understanding risk, as previous research
[@Anthony2017GloPat;@Mollentze2020VirZoo] states that locally diverse bat
communities could maintain more viruses and hence, a higher probability of
having a pathogen that could represent a risk for human health.

![Top panel: relative diversity of known bat hosts of betacoronaviruses. This
map shows that the region with the largest number of possible hosts is
South-Eastern Asia. Bottom panel: congruence between the evolutionary
distinctiveness of the hosts (grey to blue) and the viruses (grey to red). By
contrast to the richness map, this reveals that South America has the most
evolutionary distinct hosts *and* viruses, whereas South-Eastern Asia has mostly
distinct viruses. This is congruent with know results about New World bats being
evolutionary distinct, and suggests that they similarly have distinct
viruses.](figures/combined_richness.png){#fig:richness}

Nevertheless, locally diverse and virus-rich bat communities could represent an
increased risk of spillover under climate change through the creation of novel
interactions [@Carlson2022CliCha], and therefore the diversity of
betacoronavirus strains should similarly be ccounted for. In @fig:richness
(bottom), we contrast the evolutionary distinctiveness of bats and viruses --
this reveals a slightly different portrait than bat richness alone. Chiropterans
can be classified, from a macro-evolutionary standpoint, as microchiroptera and
macrochiroptera, where macrochiroptera have an older history from an
evolutionary perspective compared to macrochiroptera [@Springer2013PhyBat;
@Teeling2005MolPhy]. Specifically, we would expect that the so-called "New
World" group of bats, being more evolutionary distinct, would also have
evolutionary distinct viruses. Indeed @fig:richness (bottom) reveals it to be
the case, and this region harbors a distinct bat-betacoronavirus complex. By
contrast, South-Eastern Asia has a lot of non-evolutionary distinct bats, but
evolutionary-distinct viruses.

It is noteworthy that outside of South America, viral evolutionary
distinctiveness does not accurately tracks host diversity, with some areas
having over-distinct viruses (southern China but, oddly, not the rest of
southeast Asia). There are a number of likely explanations. First, given the
richness of bats in southeast Asia, many betacoronaviruses likely remain to be
discovered in this region. Indeed, global predictions by @Becker2022OptPre
highlight that southeast Asia is a likely hostpot of unconfirmed hosts of
betacoronaviruses, which would likely result in additiona viral discoveries.
This idea is unsurprising given the growing realization, around the emergence of
SARS-CoV-2, that a unique lineage of similar viruses are widespread in bats but
still mostly undescribed. The most distinct bats/betacoronavirus complex is
found in South America, a region with a comparatively lower number of hosts;
this matches with the isolation through viariance of the host group, and may
highlight a different co-evolutionary dynamic. Alternatively, this
distinctiveness hostpot may be a product of under-sampling: South-America is one
of the places where the fewest betacoronaviruses have been discovered
[@Anthony2017GloPat;@Allen2017GloHot; @Olival2017HosVir], resulting in sparser
phylogenetic tree, thereby artificially inflating distinctiveness. Adding more
viruses would bring the distinctiveness of known sequences down.

## The phylogeographic regions of hosts and their viruses overlap

Despite the difference in evolutionary distinctiveness globally, there are
reasons to expect that the phylogeography of bats and betacoronaviruses should
show some degree of congruence. High density of hosts sharing the same virus
(albeit possibly different strains) can drive or result from evolution of the
bat antiviral immune system, resulting in spatially distinct immunological
responses, as evidenced in several bat species [@Banerjee2020NovIns]. Immune
characteristics that allow bats to be better adapted to infection by emerging
viruses [@Gorbunova2020WorGoe; @Irving2021LesHos], in addition to being
hardcoded in their genome [@Jebb2020SixRef], may be related to a wide variety of
diets [@Banerjee2020NovIns; @MorenoSantillan2021LarGen; @Muylaert2021PreFut],
themselves likely to be driven by spatial effects, especially at the local scale
-- bats, indeed, occupy a variety of environments, and therefore display a
variety of adaptations to these environments.

![Phylogeographic regions of bats (top) and viruses (bottom) based on the joint
analysis of their occurrence and phylogenetic relatedness. The different colors
show tendencies to separate alongside the first two components of a PCoA. Note
that the PCoA for the bats and viruses are independent, and so cannot be
compared directly -- that being said, the regions can be compared across
maps.](figures/combined_biogeo.png){#fig:biogeo}

In @fig:biogeo, we show a projection of the phylogeographic signal of bats (top)
and viruses (bottom) in space; the distinct groupings (represented by different
colors symbolizing positions in the subspace formed by the first two axes of the
PCoA) are essentially equivalent between the two groups, and can be coarsely
delineated as southeast Asia, Eurasia above a northing of 25, and Africa and
south America. These results suggest that, although the evolutionary
distinctiveness of the bat/betacoronavirus complex varies spatially, the system
shows an important degree of spatial consistency, with a reduced number of
bioregions. Available information describing the spillover of zoonotic
betacoronaviruses of bat origin where data was available before and up through
the COVID-19 pandemic puts spillover events of SARS-CoV-2 in Wuhan, China;
SARS-CoV in Guangdong, China based on the presence of closest known viruses
circulating in nature, and a nearby location where serological (antibody)
evidence has indicated human exposure to SARS-like viruses [@Wang2018SerEvi];
MERS-CoV in Saudi Arabia based on index cases available from a
recently-published compendium of cases [@Ramshaw2019DatGeo]. For the latest
event, most if not all index cases are presumed to be camel-to-human
transmission, and the precise origin point (if it exists) of MERS-CoV in bats is
uncertain. Recent recombinant canine coronavirus spillover events in Haiti
[@Lednicky2021IsoNov] and Europe [@Vlasova2022AniAlp] are not relevant here, as
bats’ involvement in these cycles of transmission have been supposed to be
non-existent. These index cases fall within different phylogeographic bioregions
(@fig:biogeo), which further highlight the issue that different host-virus
sub-systems may lead to widespread emergence.

## Coevolution-informed spillover risk is different in space

As host richness, joint disctinveness, or phylogeographic structure suggest that
the bat/betacoronavirus complex is globally fragmented enough to give rise to
both different levels of risk (as evidenced by the spatial location of spillover
events) and different types of co-evolutionary dynamics, we turn to the
Geographic Mosaic Theory of Coevolution **REF** to provide a measure of risk
accounting for multiple processes. In @fig:trivariate, we overlapped three
components of spillover risk: viral sharing, *i.e.* the chance that two bats
will share viruses overall; Local Contribution to Beta Diversity, *i.e.* the
fact that a bat community is compositionally unique compared to the average
compositional similarity across the entire system; finally, host phylogenetic
diversity, *i.e.* how dispersed the bats in a location are within the tree of
life. This approach leads to the definition of broad biogeographic regions of
risk, where the same color represents the same type of risk. By way of constrat
to figures @fig:richness and @fig:biogeo, these regions do not necessarilly
overlap with previous spatial partitions of the bat/betacoronavirus complex.

![Trivariate additive mapping of the components of risk in the red/green/blue,
where high virus sharing is encoded in the blue channel, host phylogenetic
diversity in the green channel, and compositional uniqueness in the red channel.
A pixel that would maximize all measures would be a pure white (specifically
RGB(1.0, 1.0. 1.0)), and a pixel with the lowest possible values would be pure
black (specifically RGB(0.0, 0.0, 0.0)). Therefore, lighter values (the sum of
the three channels gets closer to 3) indicate higher risk, and the color
indicates the proportional distribution of the factors making up the total
risk.](figures/risk_trivariate.png){#fig:trivariate}

From the perspective of spillover risk, the most important combination of
factors is a high phylogenetic diversity of hosts with low viral sharing; this,
essentially, means that very different betacoronavirus could co-exist within the
same place. This is particularly the case given that betacoronaviruses often
evolve and even achieve host shifts through recombination, which requires the
co-occurrence of sufficiently distinct viruses to be a major driver of
emergence. In @fig:trivariate, this corresponds to yellow to pale green areas,
which are essentially limited to South-Eastern Asia, and to some part of
Sub-Saharan Africa. Adopting a geographic mosaic theory perspective on risk,
other regions of the world are of lesser concern [@fig:risk]. Our risk
decomposition does not account for viral diversity or distinctiveness. The
simple rationale behind it is that the acquisition of viral data is rarely
disconnected from the acquisition of host data; furthermore, there are more
sources of information on hosts than on viruses, allowing to develop a
host-centric perspective on risk. Areas with high bat diversity and high
turnover *may* facilitate the evolutionary radiation of viruses, matching
previous findings that the diversification of bat coronaviruses is driven
largely by host shifts (inter-genus or higher levels of cross-species
transmission) and, to a lesser degree, cospeciation and sharing [intra-genus
cross-species transmission; @Anthony2017GloPat]. This diversification is not an
actual risk factor for spillover itself, but acts downstream of a spillover
event by increasing the random chance of the emergence of a virus with the raw
genomic components required for the potential to infect humans.

![Extraction of a measure of risk based on the colorimetric space from
@fig:trivariate. The risk is a composite measure of the color value and angular
distance to the yellow hue, as defined in the methods, ranged in the unit space.
Based on this analyses, regions at high risk of spillover are southeast Asia and
Madagascar.](figures/risk_map.png){#fig:risk}

From another perspective, areas of high host uniqueness and virus sharing
(red-to-pink) could provide hotspots of betacoronavirus risk through mixing of
unique viruses (via codivergence) and in turn recombination. Under our
framework, such a hotspot was identified in Madagascar, where most bat species
are endemic following evolutionary divergence from sister species in both
African and Asian continents [*e.g.* @Shi2014DeeDiv]. Recent surveillance
[@Kettenburg2022FulGen] has identified a novel betacoronavirus (in the subgenus
*Nobecovirus*) in Madagascar-endemic pteropid bat species (*Pteropus rufus*,
*Rousettus madagascariensis*), emphasizing strong proof of principle in model
predictions.

## Human occupancy drives different levels of effective risk globally

Based on the previous result, we extracted the risk component from the composite
map (see Methods), to provide a single measure of risk varying between 0 and 1.
This measure is presented in @fig:risk. However, this maps the potential risk,
which must be weighed by the potential for contacts with humans. As a proxy for
this measure, we used the proportion of build/urban land from the EarthEnv
dataset: this is a reasonable proxy for the density of humans per unit area,
which increases the probability of pathogen spread more widely
[@Hazarie2021IntPop]. Since human activity is required to amplify the frequency
of virus encounters and thus create areas of viral amplification, mapping the
potential risk against measures of land use is required to generate a more
actionable assessment of risk. This map is presented in @fig:compound. Most of
South America and Europe are at low risk, as although densely populated,
settlements tend to be in areas with lower potential risk. However, this mapping
reveals that South-East Asia, the Indian subcontinent, and parts of sub-Saharan
Africa, are at high risk due to the overlap between built areas and bat
communities representing more opportunities for cross-species transmission of
betacoronaviruses. In looking for the origins of SARS in China, @Xu2004EpiClu
present serological evidence that strongest human-animal contact results in
higher risk of virus exposure, regardless of the animal species, but that
different types of contact had different impacts. Ideally, finer-grained
information about human activity (rather than human presence through
anthropisation) could allow to partition this risk further.

![Overlap of the percent of each pixel occupied by urbanized structures,
representing the degree of settlement, on the risk map. Darker pixels correspond
to more risk, in that the GMTC-derived risk of @fig:risk is high *and* the pixel
is densely occupied by human populations. This approach increases the relative
risk of several regions in Africa, and highlights the risk in India, southeast
China, and the arabic peninsula where areas of high to moderate risk overlap
with areas of denser population.](figures/risk_compounded.png){#fig:compound}

# Conclusion

Our study focuses largely on the biogeography of hosts. Yet, we know that
viruses with high host plasticity, that is, the ability of a given virus to
adapt to various taxonomic orders and ecological groups
[@KreuderJohnson2015SpiPan], are more likely to amplify viral spillover,
followed by secondary human-to-human transmission, and geographical spread
[@Hazarie2021IntPop]. High viral host plasticity is an especially important
trait for RNA viruses like betacoronaviruses [@Haddad2021SarPos]. Indeed, our
analysis of viral sequences reveals that Latin America is a hotspot of viral
distinctiveness, suggesting that this part of the bats-betacoronavirus system
may be undergoing independent evolutionary dynamics (related species sharing
viruses that are different from the rest of the global pool). The other hotspot
of viral distinctiveness is S.E. Asia, in which richness is high but sharing is
low; this suggests a different type of evolutionary dynamics (unrelated viruses
coevolving with evolutionarily distinct hosts, generating high diversity
locally). Both of these areas should be priority areas for sampling, especially
since @Becker2022OptPre advance that they harbor undiscovered hosts of
beta-coronaviruses. This diversity of hosts, and the mechanisms by which the
exchange of viruses occurs between species, is largely affected by the local
environmental conditions and environmental change.

Bats are important reservoir hosts for different classes of microorganisms, many
of which a threat to human health [@Letko2020BatVir; @VanBrussel2022ZooDis].
Chiropterans emerged around 64 million years ago and are one of the most diverse
mammalian orders, with an estimated richness of more than 1400 species
[@Peixoto2018SynEco; @Simmons2020BatSpe]. They exhibit a broad variety of
habitat use, behaviour, and feeding strategies, putting them at key positions in
the delivery and provisioning of several ecosystem services, tied to important
ecosystem-derived benefits to human [@Kasso2013EcoEco]. For example, bats are an
essential component of many seed-dispersal networks [@Mello2011MisPar]. Over
two-thirds of bats are know to be either obligate or facultative insectivores,
therefore playing an important role in the regulation of insect pests that can
affect crops [@Voigt2016BatAnt; @Williams-Guillen2008BatLim], and vectors of
pathogens that put a risk on human health [@Gonsalves2013MosCon;
@Gonsalves2013MosInf]. Because bats are globally distributed and have a long
evolutionary history, phylogeographic and biogeographic approaches are required
to shed light on the contemporary distribution of coevolutionary processes
between bats and the pathogens they host. Not all areas in which bats, viruses,
and human are co-occuring are facing a risk of spillover towards human
populations, and the areas in which this risk exist may not be facing risks of
the same nature and magnitude.

There are several factors that drive changes in the diversity of bats
[@Alves2018GeoVar], but human activities' effects on the ecosystem (like
modifications of land use) could significantly decrease it. Therefore, it can be
suggested that changes in the diversity of betacovs in bats are linked to their
biogeographic variation, and human population density and other anthropogenic
factors are decisive moderators for its implications in public health. With the
increase of contact between humans and potential hosts, we also increase the
risk of emergence of novel diseases [@Johnson2020GloShi], as previous studies on
RNA viruses suggest the importance of host phylogeography at the time of virus
dispersal [@Gryseels2017WheVir]. One of these scenarios where interaction
between bats and humans can occur can be seed dispersal in tropical
agroecosystems. It opens the discussion of whether the fruits thrown by bats not
only disperse seeds but could also be a source of indirect interaction between
viruses of bat origin and humans [@Deshpande2022ForFru]. This represents a
challenge for conservation strategies and disease ecology since some areas can
haveboth potential for the acquisition of zoonotic viruses and bat-human
interactions; in particular, the challenge lies in the fact that actual exposure
must then be quantified accounting for several transmission scenarios, including
both direct and indirect bat - human interaction.

**Acknowledgements**: We acknowledge that this study was conducted on land
within the traditional unceded territory of the Saint Lawrence Iroquoian,
Anishinabewaki, Mohawk, Huron-Wendat, and Omàmiwininiwak nations. This work was
supported by funding to the Viral Emergence Research Initiative (VERENA)
consortium including NSF BII 2021909 and a grant from Institut de Valorisation
des Données (IVADO). This research was enabled in part by support provided by
Calcul Québec (www.calculquebec.ca) and Compute Canada (www.computecanada.ca).
NF is funded by the NSERC BIOS² CREATE program. TP and NF are funded by the
Courtois Foundation. RLM was supported by Bryce Carmine and Anne Carmine (née
Percival), through the Massey University Foundation.

# References