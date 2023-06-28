Changes VERSION 2.40.1
----------------------

USER VISIBLE CHANGES

    - `pageSize` and `n` arguments added to `queryUniProt` to expose
    underlying API request defaults. It is recommended to set the pageSize
    to a large value e.g. 500 for large queries.

BUG FIXES AND MINOR IMPROVEMENTS

    - Fixed issue with pagination with large queries (over 25 results) in
    `queryUniProt` (@jdreyf, #23)

Changes VERSION 2.40.0
----------------------

USER VISIBLE CHANGES

    - Convert vignette from RSweave to web-based RMarkdown

BUG FIXES AND MINOR IMPROVEMENTS

    - Increase fault tolerance of unit tests
    - Add examples to `mapUniProt` documentation

Changes VERSION 2.38.0
----------------------

USER VISIBLE CHANGES

    - `UniProt.ws` uses the https://rest.uniprot.org/ API interface for queries.
    - `mapUniprot` is an exported function to directly map identifiers via
    UniProt and is used by the `select` method.
    - `allToKeys` and `allFromKeys` provide all the available `to` and `from`
    keys for mapping identifiers
    - `returnFields` provides all the possible inputs to the `columns` argument
    in the `select` and `mapUniProt` functions (ids in the `name` column)

Changes VERSION 2.27.0
----------------------

BUG FIX

    o (2.27.1) Fix bug when selecting column GENEID. The mapping mapped both
    GENEID and ENTREZ_GENE to P_ENTREZGENE. When returning columsn used match to
    identify but would only pick up first match which was ENTREZ_GENE entry.


Changes VERSION 1.2.0
--------------------

PKG FEATURES

    o UniProt.ws creates an object that talks directly to the web
    services at UniProt.  As such, it provides access to a tremendous
    library of IDs etc. directly from UniProt.

    o When the package loads there will be acces to a Uniprot.ws
    object, this object has select, keys, cols and keytypes methods
    that behave the way these methods normally do for the other
    annotation packages.  

    o One important difference from other packages is that the user
    must use the species method to set the species to their organism
    of choice (by default it is set for humans).  This is because the
    web resource is too big to return values for all these species at
    the same time.  Uniprot currently has support for over 21,000
    different species.  

    o Please see the manual pages and associated vignette for more
    detailed information about how to use this resource.


