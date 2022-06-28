Disease emergence is complex, and is driven not only by animal-human contact,
but also by the underlying evolutionary dynamics in viral reservoirs
[@Plowright2017PatZoo]. Although host richness is often used as a superficial
proxy for spillover risk [@Anthony2017GloPat; @Ruiz-Aravena2022EcoEvo], these
approaches oversimplify the relevant interspecific heterogeneity in immunology,
behavior, and other traits, and therefore overlook unique host pools that allow
for the rapid evolution of highly divergent viruses [@Agosta2010HowSpe]. In the
case of generalist pathogens like betacoronaviruses, there is conceptual and
empirical support to the idea that these community-level mechanisms are even
more important [@Power2004PatSpi], particularly given that cross-species
transmission may, as a rule, structure viral evolution more than co-divergence
with hosts [@Geoghegan2017ComAna]. This creates a disconnect between coevolutionary theory
and most existing ecological frameworks for mapping spillover risk.

The Geographic Mosaic Theory of Coevolution (GMTC) attempts to explicitly
connect microevolutionary dynamics to the macroecology and biogeography of
symbiotic interactions [@Thompson2005GeoMos]. The GMTC posits that
coevolutionary processes among pairs [@Thompson1994CoePro] or complexes
[@Janzen1980WheIt] of species are structured in space by the rippling effects of
abiotic conditions onto evolutionary mechanism, giving rise to fragmented
systems with different structure and ecologically dynamics over large spatial
extents [@Price2002MacThe]. The GMTC predicts a spatial fragmentation of
coevolutionary dynamics under the joint action of three processes
[@Gomulkiewicz2007DosDon]: coevolutionary hot- and coldspots, which appear when
the intensity of *interaction* (in terms of reciprocal fitness consequences)
varies spatially; selection mosaics, wherein the intensity of *selection* varies
across space, driven by both the biotic complexity of the community (locally
diverse hosts and viruses are more biotically complex) and the local
favorability of the environment [@Thrall2007CoeSym]; and trait remixing, which
occurs when coevolutionary dynamics are driven by by the arrival (or departure)
of _functional traits_, through changes in community composition due to
invasions, meta-community dynamics, and disperal.

Here, we apply the GMTC to explore and explain the global biogeography of
betacoronaviruses, the group that includes SARS-CoV, MERS-CoV, and SARS-CoV-2.
In their bat reservoirs, coronaviruses evolve through a mix of host jumps,
recombination among disparate lineages, and, to a lesser degree, co-divergence
with their hosts [@Anthony2017GloPat]---a mix of mechanisms that creates a
complex and nonlinear relationship between host diversity and viral emergence.
Working from a recently published database of bat hosts of betacoronaviruses, we
test whether spatial structure in bat-betacoronavirus coevolution is
identifiable at a global scale. Aiming to explain these patterns, we develop a
generalized framework for applying the GMTC to host-virus interactions, with a
specific emphasis on the potential to create independent coevolutionary dynamics
(and therefore spatial fragmentation in risk) through heterogeneity. We develop
a trivariate risk assessment system that connects each GMTC mechanism to a
quantifiable aspect of host-virus interactions: (i) viral sharing rates in host
communities, representing the strength of potential interaction between viruses
and any one host (i.e., places where viruses undergo constant host switching may
be coevolutionary coldspots); (ii) the phylogenetic diversity of hosts, as a
proxy for variation in the immunological mechanisms that antagonize viruses
(i.e., the selection mosaic); and (iii) the local uniqueness of the bat
community, representing the potential for viruses to be exposed to novel host
traits (e.g., variation in receptor sequences). Together, we argue that these
can be used to identify and map the evolutionary drivers that---in conjunction
with transmission processes (e.g., viral prevalence in reservoirs and
animal-human contact rates)--- determine disease emergence risk. 

# Results and Discussion

## Bat and betacoronavirus biogeography are broadly consistent

Most previous work has assumed that the presence or richness of key groups of
bat hosts are predictive of coronavirus diversity [@Anthony2017GloPat;
@Ruiz-Aravena2022EcoEvo]. Projecting bat and betacoronavirus phylogeny over
space (@fig:biogeo), we find support for the idea that bat community assembly is
directly responsible for a global mosaic of viral evolution. The distinct
groupings (represented by different colors, symbolizing positions in a subspace
formed by the first two phylogenetic principal components) are essentially
equivalent between the two groups, and can be coarsely delineated as (1) south
and southeast Asia, (2) east Asia (including northern China), west Asia, and the
Mediterranean coast; (3) Eurasia above a northing of 40; and (4) Africa and
south America. In some cases, this diverges from expectations about coronavirus
biogeography: for example, previous work has rarely flagged India as a region of
interest, but for both bats and betacoronaviruses, the subcontinent falls into
the same regions as the southeast Asian peninsula (and indeed, the region is
home to known bat hosts of nobecoviruses, sarbecoviruses, and merbecoviruses
[@Ruiz-Aravena2022EcoEvo]).

![Phylogeographic regions of bats (top) and viruses (bottom) based on the joint
analysis of their occurrence and phylogenetic relatedness. The different colors
show tendencies to separate alongside the first two components of a PCoA. Note
that the PCoA for the bats and viruses are independent, and so cannot be
compared directly -- that being said, the regions can be compared across
maps.](figures/combined_biogeo.png){#fig:biogeo}

Overall, these results suggest that the boundaries of bat and betacoronavirus
biogeographic regions are largely consistent. This may be surprising, given that
cospeciation plays a minor role in coronavirus diversification
[@Anthony2017GloPat], a property that would theoretically allow for substantial
broad divergence in their biogeography. However, host jumps at the family level
or higher are relatively rare and significant events in coronavirus evolutionary
history [@Anthony2017GloPat; @Latinne2020OriCro]; as a result, the mosaic of
betacoronavirus phylogeography is assembled from a set of overlapping smaller
coevolutionary systems, superimposed in space and filtered by the importance of
different subgroups in local host communities. For example, the most speciose
and cosmopolitan family of bats, the vesper bats (Vespertilionidae), are
considered the primary reservoir of merbecoviruses [@Latinne2020OriCro;
@Ruiz-Aravena2022EcoEvo]; but in the Americas, where merbecoviruses are the only
lineage present, they have only been found in other bat taxa. At the coarsest
scale, these heterogeneities are lost, and betacoronavirus biogeography tracks
the deep rifts in bat evolutionary history---but within broad regions, the
component coevolutionary systems may have very different dynamics.

## Hotspots of bat and betacoronavirus biodiversity are distinct

Bats, the second most diverse groups of mammals, are found worldwide; gradients
in their species richness generally track broader patterns of mammal diversity,
with a striking Neotropical hotspot (especially in the Amazon basin) and a
secondary hotspot centered in the southeast Asian peninsula. These hotspots of
bat diversity are generally presumed to be hotspots of viral adaptive radiation,
and therefore areas of concern for human health.
[@Anthony2017GloPat;@Olival2017HosVir] However, the hotspots of known bat
betacoronavirus hosts show a distinct pattern, with primary hotspots (both in
terms of size and higher values) of host richness situated in southeast Asia,
parts of southern Europe, and to a lesser extent parts of Africa in the -25-0
range of latitudes (@fig:richness; top). Although hundreds of species likely
host undiscovered betacoronaviruses, machine learning predictions have suggested
that these undiscovered reservoirs should follow the same diversity gradient
[@Becker2022OptPre]. In principle, these hotspots of locally-diverse, virus-rich
bat communities should drive more adaptive diversification in their viruses. 

![Top panel: relative diversity of known bat hosts of betacoronaviruses. This
map shows that the region with the largest number of possible hosts is
South-Eastern Asia. Bottom panel: congruence between the evolutionary
distinctiveness of the hosts (grey to blue) and the viruses (grey to
red).](figures/combined_richness.png){#fig:richness}

However, we find that the global pattern of betacoronavirus phylogenetic
distinctiveness is quite distinct from both bat host richness and phylogenetic
distinctiveness (@fig:richness; bottom). In contrast to the sparsity of
Neotropical betacoronavirus hosts, South America has the most evolutionary
distinct hosts *and* viruses, followed by secondary hotspots in southeast Asia
and the Rift Valley region have mostly distinct viruses. Some degree of sampling
bias may contribute to these patterns: for example, South-America is one of the
places where the fewest bat betacoronavirus sequences have been generated
[@Anthony2017GloPat;@Allen2017GloHot; @Olival2017HosVir], resulting in a sparser
phylogenetic tree, and artificially inflating distinctiveness; conversely,
disproportionate research effort in eastern China [@Cohen2022SamStr] may have
led to a more complete inventory of the local diversity of coronaviruses, again
inflating these metrics relative to underlying patterns. Even accounting for
these potential biases, though, there is obvious heterogeneity in
betacoronavirus evolutionary distinctiveness that is distinct from overall bat
diversity.

Overall, these patterns recapitulate the evolutionary history of both the order
Chiroptera and the genus _Betacoronavirus_. Horseshoe bats (Rhinolophidae) are
both the reservoirs of the SARS-like viruses  (subgenus _Sarbecovirus_), the
group of pandemic threats that have been of the greatest interest to researchers
[@Latinne2020OriCro] (and so have been sampled most intensively
[@Cohen2022SamStr]). The hotspots of host richness and viral diversity in
southeast Asia---both of which are disproportionately high, considering the
global landscape of bat species richness---are almost entirely driven by viral
adaptive radiation through host switching within this clade[@Becker2022OptPre;
@Ruiz-Aravena2022EcoEvo]. In contrast, the Neotropical hotspot of viral
distinctiveness is driven by isolation by host vicariance. Out of the four main
groups of betacoronaviruses, only the subgenus _Merbecovirus_ (MERS-like
viruses) has been found in animals in the Americas---an introduction that is
generally presumed to be ancient [@Ruiz-Aravena2022EcoEvo]. While comparatively
understudied, New World merbecoviruses have been found in the ghost-faced bats
(Mormoopidae), New World leaf-nosed bats (Phyllostomidae), and free-tailed bats
(Molossidae) [@Olival2020PosZoo]. The former two groups are endemic to the
Neotropics, while the explosive adaptive radiations of the latter two (and
particularly the phyllostomids) are responsible for the hotspot of bat diversity
in the Amazon. Together, these clades of New World bats play host to a distinct
regime of betacoronavirus coevolution. 

## Coevolution-informed emergence risk is different in space

As host richness, joint distinctiveness, or phylogeographic structure suggest
that the bat-betacoronaviruses complex is globally fragmented enough to give
rise to both different levels of risk (as evidenced by the spatial location of
spillover events) and different types of co-evolutionary dynamics, we turn to
the Geographic Mosaic Theory of Coevolution to provide a measure of risk
accounting for multiple processes. In @fig:trivariate, we overlapped three
components of spillover risk: viral sharing, *i.e.* the chance that two bats
will share viruses overall; Local Contribution to Beta Diversity, *i.e.* the
fact that a bat community is compositionally unique compared to the average
compositional similarity across the entire system; finally, host phylogenetic
diversity, *i.e.* how dispersed the bats in a location are within the tree of
life. This approach leads to the definition of broad biogeographic regions of
risk, where the same color represents the same type of risk. By way of constrat
to figures @fig:richness and @fig:biogeo, these regions do not necessarilly
overlap with previous spatial partitions of the bat-betacoronaviruses complex.

![Trivariate additive mapping of the components of risk in the red/green/blue,
where high virus sharing is encoded in the blue channel, host phylogenetic
diversity in the green channel, and compositional uniqueness in the red channel.
A pixel that would maximize all measures (highest possible risk) would be a pure
white (specifically RGB(1.0, 1.0. 1.0)), and a pixel with the lowest possible
values would be pure black (specifically RGB(0.0, 0.0, 0.0)). Therefore, lighter
values (the sum of the three channels gets closer to 3) indicate higher risk,
and the color indicates the proportional distribution of the factors making up
the total risk.](figures/risk_trivariate.png){#fig:trivariate}

From the perspective of spillover risk, the most important combination of
factors is a high phylogenetic diversity of hosts with low viral sharing; this,
essentially, means that very different betacoronaviruses could co-exist within
the same place. This is particularly the case given that betacoronaviruses often
evolve and even achieve host shifts through recombination, which requires the
co-occurrence of sufficiently distinct viruses to be a major driver of
emergence. In @fig:trivariate, this corresponds to yellow to pale green areas,
which are essentially limited to South-Eastern Asia, and to some part of
Sub-Saharan Africa. Adopting a geographic mosaic theory perspective on risk,
other regions of the world are of lesser concern (@fig:risk). Our risk
decomposition does not account for viral diversity or distinctiveness. The
simple rationale behind it is that the acquisition of viral data is rarely
disconnected from the acquisition of host data. There are more sources of
information on hosts than on viruses, allowing to develop a host-centric
perspective on risk (although this estimate would more accurate with viral
traits related to *e.g.* ability to switch hosts or pathogenic potential). Areas
with high bat diversity and high turnover *may* facilitate the evolutionary
radiation of viruses, matching previous findings that the diversification of bat
coronaviruses is driven largely by host shifts (inter-genus or higher levels of
cross-species transmission) and, to a lesser degree, cospeciation and sharing,
representing intra-genus cross-species transmission [@Anthony2017GloPat]. This
diversification is not an actual risk factor for spillover itself, but acts
downstream of a spillover event by increasing the random chance of the emergence
of a virus with the raw genomic components required for the potential to infect
humans.

![Extraction of a measure of *Betacoronavirus* spillover risk from bat hosts
based on the colorimetric space from @fig:trivariate. The risk is a composite
measure of the color value and angular distance to the yellow hue, as defined in
the methods, ranged in the unit space. Based on this analyses, regions at high
risk of spillover are southeast Asia and
Madagascar.](figures/risk_map.png){#fig:risk}

From another perspective, areas of high host uniqueness and virus sharing
(red-to-pink) could provide hotspots of *Betacoronavirus* risk through mixing of
unique viruses (via codivergence) and in turn recombination. Under our
framework, such a hotspot was identified in Madagascar, where most bat species
are endemic following evolutionary divergence from sister species in both
African and Asian continents [@Shi2014DeeDiv]. Recent surveillance
[@Kettenburg2022FulGen] has identified a novel *Betacoronavirus* (in the
subgenus *Nobecovirus*) in Madagascar-endemic pteropid bat species (*Pteropus
rufus*, *Rousettus madagascariensis*), emphasizing strong proof of principle in
model predictions.

## Human landscapes filter the geography of emergence risk

The relationship between the underlying pathogen pool and emergence risk is
mediated by both human-wildlife interfaces (the probability of spillover) and
opportunities for onward transmission (the probability that spillovers become
epidemics)[@Plowright2017PatZoo]. As a proxy for both, we finally overlaid the
risk component from the composite map (see above) with the proportion of built
land, as a proxy for a mix of habitat disturbance, potential for bat synanthropy
or contact with bridge hosts like livestock [@Rulli2021LanUse; @Cui2019OriEvo],
and human population density and connectivity [@Plowright2017PatZoo;
@Muylaert2022PreFut; @Hassell2017UrbDis] (@fig:compound). Accounting for these
factors, most of South America and Europe are at comparatively lower risk,
as--although densely populated--settlements tend to be in areas with lower
potential risk. Conversely, regions like Malaysia and the northern coast of
Australia have a high evolutionary risk component, but should represent a
relatively lower effective risk due to low human density. However, southeast
Asia, the Indian subcontinent, and scattered hotspots in sub-Saharan Africa are
at high risk due to the overlap between human populations and natural
opportunities for cross-species transmission of betacoronaviruses. 

![Overlap of the percent of each pixel occupied by urbanized structures,
representing the degree of settlement, on the spillover risk map (where the risk
comes only from wildlife, and ignores multi-hosts chains of transmissions
including non-bats hosts). Darker pixels correspond to more risk, in that the
GMTC-derived risk of @fig:risk is high *and* the pixel is densely occupied by
human populations. This approach increases the relative risk of several regions
in Africa, and highlights the risk in India, southeast China, and the Arabian
peninsula where areas of high to moderate risk overlap with areas of denser
population.](figures/risk_compounded.png){#fig:compound}

Reassuringly, these predictions correspond to the geographic origins of the
three bat-origin coronaviruses that have recently emerged in human populations.
While available information puts the spillover of SARS-CoV-2 in a live animal
market in Wuhan, China, the ultimate origin of the virus is almost certainly in
a divergent lineage of sarbecoviruses from the Indochinese peninsula that was
poorly characterized prior to the pandemic [@Worobey2022HuaMar;
@Temmam2022BatCor; @Boni2020EvoOri]. Similarly, the SARS-CoV outbreak began in
Guangdong province in 2002, reaching humans through small carnivore bridge
hosts, but was eventually traced back to a set of likely progenitor viruses
found in cave-dwelling horseshoe bats in Yunnan province [@Hu2017DisRic];
nearby, antibody evidence has indicated human exposure to SARS-like viruses
[@Wang2018SerEvi].  MERS-CoV was originally detected in Saudi Arabia,
accompanied by a nearly identical virus sequenced from an Egyptian tomb bat
(_Taphozous perforatus_)[@Memish2013MidEas], but is widespread in camels in East
Africa and the Middle East, and may have reached its bridge host decades earlier
than originally supposed [@Muller2014MerCor]; as a result, the geography of the
original bat-to-camel transmission is still widely regarded as uncertain. All of
these are broadly consistent with the risk factors we identify. Notably, India
and west Africa are additional hotspots that have yet to experience the
emergence of a bat coronavirus into human populations, but may still be at
risk---particularly given known gaps in bat surveillance [@Cohen2022SamStr], and
a dense population in both regions with global connectivity. In any of these
regions, surveillance on viral reservoirs can be paired with targeted monitoring
of high-risk human populations (i.e., those with regular wildlife contact
[@Xu2004EpiClu]) for maximum impact.

# Conclusion

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
therefore actively contributing for agricultural pest control [@Voigt2016BatAnt;
@Williams-Guillen2008BatLim], and vectors of pathogens that put a risk on human
health [@Gonsalves2013MosCon; @Gonsalves2013MosInf]. Because bats are globally
distributed and have a long evolutionary history, phylogeographic and
biogeographic approaches are required to shed light on the contemporary
distribution of coevolutionary processes between bats and the pathogens they
host. Not all areas in which bats, viruses, and human are co-occuring are facing
a risk of spillover towards human populations, and the areas in which this risk
exist may not be facing risks of the same nature and magnitude.

Here, we propose a simple freamework with broad explanatory power that helps
contextualize discoveries like highly divergent nobecoviruses in Madagascar and
the previously-neglected adaptive radiation of sarbecoviruses outside of
southern China and throughout southeast Asia. In doing so, it advances
ecological theory beyond the current state of the art for global maps of
emergence risk. For example, previous studies that have used host richness as
proxy have predicted a high diversity of unsampled bat viruses
[@Olival2017HosVir], bat coronaviruses [@Anthony2017GloPat], and even
specifically betacoronaviruses [@Becker2022OptPre] in both the Amazon and
southeast Asia. While we find that both regions are characterized by highly
divergent host and viral communities, our framework identifies key differences
between the regions. We find that Latin America is a hotspot of both host and
viral distinctiveness, suggesting that this branch of the bat-betacoronavirus
complex may be undergoing independent evolutionary dynamics from the rest of the
global pool, but with limited potential for viral diversification--- a finding
that is supported by previous work indicating a higher rate of codivergence in
Latin America [@Anthony2017GloPat]. In contrast, in southeast Asia, host
richness and viral distinctiveness are high but sharing is low; this suggests a
different type of evolutionary dynamics that could generate high local diversity
of viruses through host switching and viral recombination (see *e.g.*
[@Latinne2020OriCro], as well as the discovery of recombinant viruses that share
genetic material from both the SARS-CoV and SARS-CoV-2 branches of the
Sarbecovirus lineage [@Wu2021ComSur]). Both of these regions are priority areas
for sampling, especially given predictions that they contain many bat hosts of
undiscovered betacoronaviruses [@Becker2022OptPre; @Cohen2022SamStr]. However,
both the evolutionary and ecological aspects of emergence risk are likely higher
in southeast Asia---a fact that will only become more relevant, as bats track
shifting climates and exchange viruses with other species, creating a hotspot of
cross-species transmission unique to the region [@Carlson2022CliCha].

The diversity and diversification potential of bats responds to anthropogenic
factors others than shifting climates [@Alves2018GeoVar]. Land use changes could
significantly decrease bat suitability, notably through effects on diet and
availability of habitats [@Treitler2016EffLoc]. As our results establish that
the diversification of bats betacoronaviruses happens on top of processes
affecting hosts, biogeographic variation in human population density and
anthropogenic disturbances may feed into co-evolutionary dynamics. Increase in
humans-hosts contacts also increase the risk of emergence of novel diseases
[@Johnson2020GloShi], so does the changes in landscape connectivity at
local/regional scales [@Gryseels2017WheVir]. This represents a challenge for
both conservation strategies and disease ecology: some areas can a high
emergence risk and more potential for the acquisition of zoonotic viruses
through bat-human encounters [@Amman2011InvRol]. In particular, the challenge
ahead lies in the need to quantify actual exposure (and risk)  accounting for
several transmission scenarios, including both direct and indirect bat - human
interactions, and feeding back into the provision of ecosystem services by bats.

**Acknowledgements**: We acknowledge that this study was conducted on land
within the traditional unceded territory of the Saint Lawrence Iroquoian,
Anishinabewaki, Mohawk, Huron-Wendat, and Omàmiwininiwak nations. This work was
supported by funding to the Viral Emergence Research Initiative (VERENA)
consortium including NSF BII 2021909 and a grant from Institut de Valorisation
des Données (IVADO). This research was enabled in part by support provided by
Calcul Québec (www.calculquebec.ca) and Compute Canada (www.computecanada.ca).
NF is funded by the NSERC BIOS² CREATE program. TP and NF are funded by the
Courtois Foundation. RLM was supported by Bryce Carmine and Anne Carmine (née
Percival), through the Massey University Foundation. DJB was supported by the
National Institute of General Medical Sciences of the National Institutes of
Health (P20GM134973).

\newpage

# Methods

## Known *Betacoronavirus* hosts

We downloaded the data on bats hosts of *Betacoronavirus* from
`https://www.viralemergence.org/betacov` on Apr. 2022 [@Becker2022OptPre], and
filtered it to "known" hosts (established before the emergence of SARS-CoV-2)
and "novel" hosts (confirmed through sampling and competence assays since the
initial data collection). The original database was assembled by a combination
of data mining and literature surveys, including automated alerts on the "bats"
and "coronavirus" keywords to identify novel empirical evidence of
bats-betacoronaviruses associations; this yielded a total of 126 known hosts, 47
of which were novel hosts.

## Bat occurrences

We downloaded the rangemap of every current bat species that was classified as
an empirically documented host of *Betacoronavirus* from the previous step,
according to recent IUCN data [@IUCN2021IucRed]. The range maps were
subsequently rasterized using the `rasterize` function from `GDAL`
[@RouaultEven2022GdaOgr] at a resolution of approximately 100kmx100km. For every
pixel in the resulting raster where at least one bat host of *Betacoronavirus*
was present, we extract the species pool (list of all known bat hosts),
which was used to calculate the following risk assessment components: bat
phylogenetic diversity, bat compositional uniqueness, and predicted viral
sharing risk.

## Bat phylogenetic diversity

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

## Bat compositional uniqueness

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

For all bat hosts of *Betacoronavirus*, we extracted their predicted viral
sharing network, generated from a previously published generalized additive
mixed model of virus sharing by a tensor function of phylogenetic distance and
geographic range overlap across mammals [@Albery2020PreGlo]. This network stores
pairwise values of viral community similarity. To project viral sharing values
into a single value for every pixel, we averaged the pairwise scores. High
values of the average sharing propensity means that this specific extant bat
assemblage is likely to be proficient at exchanging viruses.

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
[@Gomulkiewicz2000HotSpo], in that this allows multiple independent host-virus
coevolutionary dynamics to take place in the same location. In the colorimetric
space, this correspond to yellow -- because the HSV space is more amenable to
calculations for feature extraction [@Keke2010StuSki], we measured the risk
level by calculating the angular distance of the hue of each pixel to a
reference value of 60 (yellow), and weighted this risk level by the value
component. Specifically, given a pixel with colorimetric coordinates $(h,s,v)$,
its ranged weighted risk value is

$$
v\times\left[1-\frac{\left|\text{atan}\left(\text{cos}(\text{rad}(h)), \text{sin}(\text{rad}(h))\right) - X\right|}{2\pi}\right]\,,
$$

where X is $\text{atan}\left(\text{cos}(\text{rad}(60)),
\text{sin}(\text{rad}(60))\right)$, a constant approximately equal to $0.5235$.

## Viral phylogeography and evolutionary diversification

To next represent phylogeography of betacoronaviruses in bats, we aggregated and
analyzed betacoronavirus sequence data. We used the following query to pull all
*Betacoronavirus* sequence data from the GenBank Nucleotide database except
SARS-CoV-2; ("Betacoronavirus"[Organism] OR betacoronavirus[All Fields]) NOT
("Severe acute respiratory syndrome coronavirus 2"[Organism] OR sars-cov-2[All
Fields]). We added a single representative sequence for SARS-CoV-2 and manually
curated to remove sequences without the RNA-dependent RNA polymerase (RdRp)
sequence or that contained words indicating recombinant or laboratory strains
including “patent”, “mutant”, “GFP”, and “recombinant”. We filtered
over-represented taxa including betacoronavirus 1, hCoV-OC43, Middle East
respiratory syndrome coronavirus, Murine hepatitis virus, and hCoV-HKU1. Curated
betacoronavirus RdRp sequences were then aligned using MAFFT [@Katoh2013MafMul]
v1.4.0 (Algorithm FFT-NS-2, Scoring matrix 200PAM / k=2, gap open penalty 1.53m
offset value 0.123) and a maximum likelihood tree reconstructed in IQ-TREE
[@Nguyen2015IqtFas] v1.6.12 with ModelFinder [@Kalyaanamoorthy2017ModFas]
ultrafast bootstrap approximation [@Hoang2018UfbImp] with a general time
reversible model with empirical base frequencies and the
5-discrete-rate-category FreeRaye model of nucleotide substitution (GTR+F+R5).

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
between the different viruses in the betacoronaviruses tree by using the
cophenetic function in `ape` [@Paradis2019ApeEnv]; subsequently, we take a
principal components analysis of that distance matrix (readily interchangeable
for NMDS in this case) to project the viral tree into an n-dimensional space. We
then take the first two principal components and, as with the evolutionary
distinctiveness analysis, aggregated these to a mean host value and projected
them using a four-color bivariate map.

\newpage

# References
