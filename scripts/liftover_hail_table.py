import hail as hl

hl.init(log="/Users/xx20081/git/dcm_pgs/output/hail.log")
hl.default_reference("GRCh38")

rg38 = hl.get_reference("GRCh38")
rg37 = hl.get_reference("GRCh37")

rg38.add_liftover('/Users/xx20081/git/dcm_pgs/data/grch38_to_grch37.over.chain.gz', rg37)

ht = hl.read_table("/Users/xx20081/Downloads/gnomad.v3.1.pca_loadings.ht")
ht.describe()
ht.show(3)


# liftover
ht37 = ht.annotate(locus37 = hl.liftover(ht.locus, "GRCh37"))
ht37 = ht37.filter(hl.is_defined(ht37.locus37))
ht37 = ht37.key_by(locus=ht37.locus37, alleles=ht37.alleles)
ht37 = ht37.drop("locus37")
ht37.describe()


print("GRCh38 rows:", ht.count())
print("GRCh37 rows:", ht37.count())

# save for rsid annotation
out = ht37.select(
    chr = ht37.locus.contig,
    pos = ht37.locus.position,
    ref = ht37.alleles[0],
    alt = ht37.alleles[1]
)
out.export("/Users/xx20081/git/dcm_pgs/data/gnomad_pca_variants_chr_pos_ref_alt.tsv")

# save HAIL table
ht37.write(
    "/Users/xx20081/git/dcm_pgs/data/gnomad.v3.1.pca_loadings_grch37.ht",
    overwrite=True
)
