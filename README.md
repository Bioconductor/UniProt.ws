
# `UniProt.ws`

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("UniProt.ws")
```

## Configuring `UniProt.ws`

The *[UniProt.ws](https://bioconductor.org/packages/3.21/UniProt.ws)*
package provides a `select` interface to the UniProt web service.

``` r
suppressPackageStartupMessages({
    library(UniProt.ws)
})
up <- UniProt.ws(taxId=9606)
```

If you already know about the select interface, you can immediately
learn about the various methods for this object by just looking it’s the
help page.

``` r
help("UniProt.ws")
```

When you load the
*[UniProt.ws](https://bioconductor.org/packages/3.21/UniProt.ws)*
package, it creates a `UniProt.ws` object. If you look at the object you
will see some helpful information about it.

``` r
up
#> UniProt.ws interface object:
#> Taxonomy ID: 9606
#> Species name: Homo sapiens (Human)
#> List species with 'availableUniprotSpecies()'
```

By default, you can see that the `UniProt.ws` object is set to retrieve
records from Homo sapiens. But you can change that of course. In order
to change it, you first need to look up the appropriate taxonomy ID for
the species that you are interested in. Uniprot provides support for
over 20 thousand species, so there are a few to choose from! In order to
make this easier, we have provided the helper function
`availableUniprotSpecies` which will list all the supported species
along with their taxonomy ids. When you call the
`availableUniprotSpecies` function, it’s recommended that you make use
of the pattern argument to limit your queries like this:

``` r
availableUniprotSpecies(pattern="musculus")
#>       kingdom Taxon Node                     Official (scientific) name
#> ANTMS       E     520121                            Anthocoris musculus
#> ANTMU       E     208057                           Anthoscopus musculus
#> APOMU       E     238007                                Apomys musculus
#> BAIMU       E     213557                               Baiomys musculus
#> BALMU       E       9771                          Balaenoptera musculus
#> BLEMU       E     197864                           Blepharisma musculus
#> MOUSE       E      10090                                   Mus musculus
#> MUSMB       E      35531                        Mus musculus bactrianus
#> MUSMC       E      10091                         Mus musculus castaneus
#> MUSMM       E      57486                        Mus musculus molossinus
#> MUSMS       E     186842                     Mus musculus x Mus spretus
#> MUSMX       E     477816 Mus musculus musculus x Mus musculus castaneus
#> POVM1       V    1891730                    Mus musculus polyomavirus 1
```

Once you have learned the taxonomy ID for the species of interest, you
can then change the taxonomy id for the `UniProt.ws` object using
`taxId` setter or by calling the constructor for `UniProt.ws`

``` r
mouseUp <- UniProt.ws(10090)
mouseUp
#> UniProt.ws interface object:
#> Taxonomy ID: 10090
#> Species name: Mus musculus (Mouse)
#> List species with 'availableUniprotSpecies()'
```

As you can see the species is different for the `mouseUp` new object.

## Using `UniProt.ws`

Once you are safisfied that you have an `uniport.ws` that is using the
appropriate organsims, you can make use of the standard set of methods
in a `select` interface. Specifically: `columns`, `keytypes`, `keys` and
`select`.

You will probably notice that there are a large number of columns that
can be retrieved.

``` r
head(keytypes(up))
#> [1] "Allergome"     "ArachnoServer" "Araport"       "BioCyc"       
#> [5] "BioGRID"       "BioMuta"
```

And most (but not all) of these fields can also be used as keytypes.

``` r
head(columns(up))
#> [1] "absorption"              "accession"              
#> [3] "annotation_score"        "cc_activity_regulation" 
#> [5] "cc_allergen"             "cc_alternative_products"
```

If necessary you can also look up the keys of a given type. But please
be warned that the web service is slow at this particular kind of
lookup. So if you really want to do this kind of operation you are
probably going to want to save the result to your R session.

``` r
egs <- keys(up, "GeneID")
```

Finally, you can loop up whatever combinations of columns, keytypes and
keys that you need when using `select`.

*Note*. ‘ENTREZ_GENE’ is now ‘GeneID’

``` r
keys <- c("1","2")
columns <- c("xref_pdb", "xref_hgnc", "sequence")
kt <- "GeneID"
res <- select(up, keys, columns, kt)
res
#>   From  Entry                                                               PDB
#> 1    1 P04217                                                              <NA>
#> 2    1 V9HWD8                                                              <NA>
#> 3    2 P01023 1BV8;2P9R;6TAV;7O7L;7O7M;7O7N;7O7O;7O7P;7O7Q;7O7R;7O7S;7VON;7VOO;
#>      HGNC
#> 1 HGNC:5;
#> 2    <NA>
#> 3 HGNC:7;
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Sequence
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MSMLVVFLLLWGVTWGPVTEAAIFYETQPSLWAESESLLKPLANVTLTCQAHLETPDFQLFKNGVAQEPVHLDSPAIKHQFLLTGDTQGRYRCRSGLSTGWTQLSKLLELTGPKSLPAPWLSMAPVSWITPGLKTTAVCRGVLRGVTFLLRREGDHEFLEVPEAQEDVEATFPVHQPGNYSCSYRTDGEGALSEPSATVTIEELAAPPPPVLMHHGESSQVLHPGNKVTLTCVAPLSGVDFQLRRGEKELLVPRSSTSPDRIFFHLNAVALGDGGHYTCRYRLHDNQNGWSGDSAPVELILSDETLPAPEFSPEPESGRALRLRCLAPLEGARFALVREDRGGRRVHRFQSPAGTEALFELHNISVADSANYSCVYVDLKPPFGGSAPSERLELHVDGPPPRPQLRATWSGAVLAGRDAVLRCEGPIPDVTFELLREGETKAVKTVRTPGAAANLELIFVGPQHAGNYRCRYRSWVPHTFESELSDPVELLVAES
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MSMLVVFLLLWGVTWGPVTEAAIFYETQPSLWAESESLLKPLANVTLTCQAHLETPDFQLFKNGVAQEPVHLDSPAIKHQFLLTGDTQGRYRCRSGLSTGWTQLSKLLELTGPKSLPAPWLSMAPVSWITPGLKTTAVCRGVLRGVTFLLRREGDHEFLEVPEAQEDVEATFPVHQPGNYSCSYRTDGEGALSEPSATVTIEELAAPPPPVLMHHGESSQVLHPGNKVTLTCVAPLSGVDFQLRRGEKELLVPRSSTSPDRIFFHLNAVALGDGGHYTCRYRLHDNQNGWSGDSAPVELILSDETLPAPEFSPEPESGRALRLRCLAPLEGARFALVREDRGGRRVHRFQSPAGTEALFELHNISVADSANYSCVYVDLKPPFGGSAPSERLELHVDGPPPRPQLRATWSGAVLAGRDAVLRCEGPIPDVTFELLREGETKAVKTVRTPGAAANLELIFVGPQHAGNYRCRYRSWVPHTFESELSDPVELLVAES
#> 3 MGKNKLLHPSLVLLLLVLLPTDASVSGKPQYMVLVPSLLHTETTEKGCVLLSYLNETVTVSASLESVRGNRSLFTDLEAENDVLHCVAFAVPKSSSNEEVMFLTVQVKGPTQEFKKRTTVMVKNEDSLVFVQTDKSIYKPGQTVKFRVVSMDENFHPLNELIPLVYIQDPKGNRIAQWQSFQLEGGLKQFSFPLSSEPFQGSYKVVVQKKSGGRTEHPFTVEEFVLPKFEVQVTVPKIITILEEEMNVSVCGLYTYGKPVPGHVTVSICRKYSDASDCHGEDSQAFCEKFSGQLNSHGCFYQQVKTKVFQLKRKEYEMKLHTEAQIQEEGTVVELTGRQSSEITRTITKLSFVKVDSHFRQGIPFFGQVRLVDGKGVPIPNKVIFIRGNEANYYSNATTDEHGLVQFSINTTNVMGTSLTVRVNYKDRSPCYGYQWVSEEHEEAHHTAYLVFSPSKSFVHLEPMSHELPCGHTQTVQAHYILNGGTLLGLKKLSFYYLIMAKGGIVRTGTHGLLVKQEDMKGHFSISIPVKSDIAPVARLLIYAVLPTGDVIGDSAKYDVENCLANKVDLSFSPSQSLPASHAHLRVTAAPQSVCALRAVDQSVLLMKPDAELSASSVYNLLPEKDLTGFPGPLNDQDNEDCINRHNVYINGITYTPVSSTNEKDMYSFLEDMGLKAFTNSKIRKPKMCPQLQQYEMHGPEGLRVGFYESDVMGRGHARLVHVEEPHTETVRKYFPETWIWDLVVVNSAGVAEVGVTVPDTITEWKAGAFCLSEDAGLGISSTASLRAFQPFFVELTMPYSVIRGEAFTLKATVLNYLPKCIRVSVQLEASPAFLAVPVEKEQAPHCICANGRQTVSWAVTPKSLGNVNFTVSAEALESQELCGTEVPSVPEHGRKDTVIKPLLVEPEGLEKETTFNSLLCPSGGEVSEELSLKLPPNVVEESARASVSVLGDILGSAMQNTQNLLQMPYGCGEQNMVLFAPNIYVLDYLNETQQLTPEIKSKAIGYLNTGYQRQLNYKHYDGSYSTFGERYGRNQGNTWLTAFVLKTFAQARAYIFIDEAHITQALIWLSQRQKDNGCFRSSGSLLNNAIKGGVEDEVTLSAYITIALLEIPLTVTHPVVRNALFCLESAWKTAQEGDHGSHVYTKALLAYAFALAGNQDKRKEVLKSLNEEAVKKDNSVHWERPQKPKAPVGHFYEPQAPSAEVEMTSYVLLAYLTAQPAPTSEDLTSATNIVKWITKQQNAQGGFSSTQDTVVALHALSKYGAATFTRTGKAAQVTIQSSGTFSSKFQVDNNNRLLLQQVSLPELPGEYSMKVTGEGCVYLQTSLKYNILPEKEEFPFALGVQTLPQTCDEPKAHTSFQISLSVSYTGSRSASNMAIVDVKMVSGFIPLKPTVKMLERSNHVSRTEVSSNHVLIYLDKVSNQTLSLFFTVLQDVPVRDLKPAIVKVYDYYETDEFAIAEYNAPCSKDLGNA
```

## `sessionInfo()`

``` r
sessionInfo()
#> R Under development (unstable) (2024-11-01 r87285)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] UniProt.ws_2.47.4 BiocStyle_2.35.0 
#> 
#> loaded via a namespace (and not attached):
#>  [1] rappdirs_0.3.3          generics_0.1.3          RSQLite_2.3.9          
#>  [4] hms_1.1.3               digest_0.6.37           magrittr_2.0.3         
#>  [7] evaluate_1.0.1          fastmap_1.2.0           blob_1.2.4             
#> [10] jsonlite_1.8.9          progress_1.2.3          AnnotationDbi_1.69.0   
#> [13] GenomeInfoDb_1.43.2     DBI_1.2.3               BiocManager_1.30.25    
#> [16] httr_1.4.7              purrr_1.0.2             UCSC.utils_1.3.0       
#> [19] Biostrings_2.75.3       codetools_0.2-20        httr2_1.0.7            
#> [22] cli_3.6.3               rlang_1.1.4             crayon_1.5.3           
#> [25] dbplyr_2.5.0            XVector_0.47.1          Biobase_2.67.0         
#> [28] bit64_4.5.2             withr_3.0.2             cachem_1.1.0           
#> [31] yaml_2.3.10             BiocBaseUtils_1.9.0     tools_4.5.0            
#> [34] memoise_2.0.1           dplyr_1.1.4             filelock_1.0.3         
#> [37] GenomeInfoDbData_1.2.13 BiocGenerics_0.53.3     curl_6.0.1             
#> [40] rjsoncons_1.3.1         vctrs_0.6.5             R6_2.5.1               
#> [43] png_0.1-8               stats4_4.5.0            lifecycle_1.0.4        
#> [46] BiocFileCache_2.15.0    KEGGREST_1.47.0         S4Vectors_0.45.2       
#> [49] IRanges_2.41.2          bit_4.5.0.1             pkgconfig_2.0.3        
#> [52] pillar_1.10.0           glue_1.8.0              tidyselect_1.2.1       
#> [55] xfun_0.49               tibble_3.2.1            rstudioapi_0.17.1      
#> [58] knitr_1.49              AnVILBase_1.1.0         htmltools_0.5.8.1      
#> [61] rmarkdown_2.29          compiler_4.5.0          prettyunits_1.2.0
```
